use core::str;
use std::cmp;
use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::Result;
use askama::Template;
use bio::io::fasta;
use bio_types::genome::AbstractLocus;
use derive_builder::Builder;
use itertools::Itertools;
use regex::Regex;
use rust_htslib::bam::Read as BamRead;
use rust_htslib::{bam, bcf, bcf::Read};

use crate::calling::variants::preprocessing::haplotype_feature_index::HaplotypeFeatureIndex;
use crate::errors;
use crate::utils;
use crate::utils::anonymize::Anonymizer;
use crate::utils::collect_variants::VariantInfo;
use crate::variants::model::{HaplotypeIdentifier, Variant};
use crate::variants::sample;
use crate::{cli, reference};

lazy_static! {
    static ref TESTCASE_RE: Regex =
        Regex::new(r"^(?P<chrom>[^:]+):(?P<pos>\d+)(:(?P<idx>\d+))?$").unwrap();
}

lazy_static! {
    static ref BND_ALLELE: Regex = Regex::new(r"(\d+):(\d+)").unwrap();
}

#[derive(Template)]
#[template(path = "testcase.yml", escape = "none")]
struct TestcaseTemplate {
    samples: HashMap<String, Sample>,
    candidate: String,
    scenario: Option<String>,
    ref_path: String,
    mode: Mode,
    purity: Option<f64>,
}

#[derive(Debug, Clone)]
struct Sample {
    path: String,
    properties: String,
    options: String,
}

#[derive(Debug, Clone, Copy, EnumString, Display)]
pub enum Mode {
    TumorNormal,
    Generic,
}

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Testcase {
    #[builder(setter(into))]
    prefix: PathBuf,
    #[builder(private)]
    chrom_name: Option<Vec<u8>>,
    #[builder(private)]
    pos: Option<u64>,
    #[builder(private)]
    idx: usize,
    #[builder(private)]
    reference_buffer: reference::Buffer,
    candidates: PathBuf,
    #[builder(private)]
    bams: HashMap<String, PathBuf>,
    scenario: Option<PathBuf>,
    #[builder(private)]
    options: HashMap<String, String>,
    #[builder(default = "None")]
    purity: Option<f64>,
    mode: Mode,
    anonymize: bool,
}

pub(crate) fn modify_bnd_alleles(
    alleles: &[&[u8]],
    chromosomal_regions: &HashMap<Vec<u8>, (u64, u64)>,
) -> Result<Vec<Vec<u8>>> {
    let mod_alleles = alleles
        .iter()
        .map(|allele| {
            BND_ALLELE
                .replace(str::from_utf8(allele).unwrap(), |caps: &regex::Captures| {
                    let chrom = caps[1].as_bytes();
                    let pos: u64 = caps[2].parse().unwrap();
                    let ref_start = chromosomal_regions.get(chrom).unwrap().0;
                    let mpos = pos - ref_start;
                    format!("{}:{}", &caps[1], mpos)
                })
                .as_bytes()
                .to_vec()
        })
        .collect();
    Ok(mod_alleles)
}

impl TestcaseBuilder {
    pub(crate) fn reference(self, path: impl AsRef<Path> + std::fmt::Debug) -> Result<Self> {
        Ok(self.reference_buffer(reference::Buffer::from_path(&path, 1)?))
    }

    pub(crate) fn locus(self, locus: &str) -> Result<Self> {
        if locus == "all" {
            Ok(self.chrom_name(None).pos(None).idx(0))
        } else if let Some(captures) = TESTCASE_RE.captures(locus) {
            let chrom_name = captures
                .name("chrom")
                .unwrap()
                .as_str()
                .as_bytes()
                .to_owned();
            let mut pos: u64 = captures.name("pos").unwrap().as_str().parse()?;
            pos -= 1;
            let idx = if let Some(m) = captures.name("idx") {
                let idx: usize = m.as_str().parse()?;
                idx - 1
            } else {
                0
            };

            Ok(self.chrom_name(Some(chrom_name)).pos(Some(pos)).idx(idx))
        } else {
            Err(errors::Error::InvalidLocus.into())
        }
    }

    pub(crate) fn register_sample(
        mut self,
        name: &str,
        bam: impl AsRef<Path>,
        mut options: cli::Varlociraptor,
    ) -> Result<Self> {
        if self.bams.is_none() {
            self = self.bams(HashMap::new());
        }

        self.bams
            .as_mut()
            .unwrap()
            .insert(name.to_owned(), bam.as_ref().to_owned());

        if self.options.is_none() {
            self = self.options(HashMap::new());
        }

        if let cli::Varlociraptor::Preprocess {
            kind:
                cli::PreprocessKind::Variants {
                    ref mut reference,
                    ref mut candidates,
                    ref mut bam,
                    ref mut output,
                    ..
                },
        } = options
        {
            *reference = "?".into();
            *candidates = "?".into();
            *bam = "?".into();
            *output = Some("?".into());
        } else {
            unreachable!();
        }

        self.options
            .as_mut()
            .unwrap()
            .insert(name.to_owned(), serde_json::to_string(&options)?);

        Ok(self)
    }
}

impl Testcase {
    fn candidate_reader(&self) -> Result<bcf::Reader> {
        Ok(bcf::Reader::from_path(&self.candidates)?)
    }

    fn variants(&mut self) -> Result<Vec<bcf::Record>> {
        let mut candidate_reader = self.candidate_reader()?;

        // get variant
        let rid = if let Some(name) = self.chrom_name.as_ref() {
            Some(candidate_reader.header().name2rid(name)?)
        } else {
            None
        };

        let mut found = vec![];
        for res in candidate_reader.records() {
            let rec = res?;
            if let Some(rec_rid) = rec.rid() {
                if let Some(rid) = rid {
                    if rec_rid == rid && rec.pos() as u64 == self.pos.unwrap() {
                        found.push(rec);
                        break;
                    }
                } else {
                    // add all variants
                    found.push(rec);
                }
            }
        }
        if found.is_empty() {
            Err(errors::Error::NoCandidateFound.into())
        } else {
            if rid.is_some() {
                // If not collecting all records, fetch all breakends belonging to the same events.
                let mut breakend_index = None;
                let mut aux_breakends = Vec::new();

                for rec in &mut found {
                    if utils::is_bnd(rec)? {
                        if let Some(event) = HaplotypeIdentifier::from(rec)? {
                            // METHOD: for breakend events, collect all the other breakends.
                            if breakend_index.is_none() {
                                breakend_index =
                                    Some(HaplotypeFeatureIndex::new(&self.candidates)?);
                            }
                            let breakend_index = breakend_index.as_ref().unwrap();
                            let last_idx = breakend_index.last_record_index(&event).unwrap();

                            let mut candidate_reader = self.candidate_reader()?;
                            for (i, res) in candidate_reader.records().enumerate() {
                                let mut other_rec = res?;
                                if let Some(other_event) =
                                    HaplotypeIdentifier::from(&mut other_rec)?
                                {
                                    if event == other_event
                                        && (other_rec.contig() != rec.contig()
                                            || other_rec.pos() != rec.pos())
                                    {
                                        // This is another record of the same event.
                                        aux_breakends.push(other_rec);
                                    }
                                }
                                if i == last_idx {
                                    // Last record of this event processed, stop.
                                    break;
                                }
                            }
                        } else {
                            info!("Skipping collection of mate breakends because EVENT tag is not specified, which is unsupported at the moment.")
                        }
                    }
                }

                found.extend(aux_breakends);
            }

            Ok(found)
        }
    }

    pub(crate) fn collect_start_end(
        &mut self,
        candidate_variant_info: VariantInfo,
        pos: u64,
    ) -> Result<(u64, u64)> {
        let (start, end) = match candidate_variant_info.variant() {
            Variant::Deletion(l) => (pos.saturating_sub(1000), pos + { *l } + 1000),
            Variant::Insertion(ref seq) => {
                (pos.saturating_sub(1000), pos + seq.len() as u64 + 1000)
            }
            Variant::Methylation() => (pos.saturating_sub(100), pos + 2 + 100),
            Variant::Snv(_) => (pos.saturating_sub(100), pos + 1 + 100),
            Variant::Mnv(ref bases) => (pos.saturating_sub(100), pos + bases.len() as u64 + 100),
            Variant::Breakend { .. } => {
                (pos.saturating_sub(1000), pos + 1 + 1000) // TODO collect entire breakend event!
            }
            Variant::Inversion(l) => (pos.saturating_sub(1000), pos + { *l } + 1000),
            Variant::Duplication(l) => (pos.saturating_sub(1000), pos + { *l } + 1000),
            Variant::Replacement { ref ref_allele, .. } => (
                pos.saturating_sub(1000),
                pos + ref_allele.len() as u64 + 1000,
            ),
            Variant::None => (pos.saturating_sub(100), pos + 1 + 100),
        };
        Ok((start, end))
    }

    pub(crate) fn extend_chromosomal_regions(
        &self,
        chromosomal_regions: &HashMap<Vec<u8>, (u64, u64)>,
    ) -> Result<HashMap<Vec<u8>, (u64, u64)>> {
        let mut extended_chromosomal_regions = HashMap::new();
        for path in self.bams.values() {
            let mut bam_reader = bam::IndexedReader::from_path(path)?;

            for (chrom_name, (start, end)) in chromosomal_regions.clone() {
                let mut ref_start = start;
                let mut ref_end = end;
                let tid: u32 = bam_reader.header().tid(&chrom_name).unwrap();
                // first pass, extend reference interval
                bam_reader.fetch((tid, start, end))?;
                for res in bam_reader.records() {
                    let rec = res?;
                    let seq_len = rec.seq().len() as u64;
                    ref_start = cmp::min((rec.pos() as u64).saturating_sub(seq_len), ref_start);
                    ref_end = cmp::max(rec.cigar().end_pos() as u64 + seq_len, ref_end);
                }
                extended_chromosomal_regions.insert(chrom_name.clone(), (ref_start, ref_end));
            }
        }
        Ok(extended_chromosomal_regions)
    }

    pub(crate) fn write(&mut self) -> Result<()> {
        let mut anonymizer = Anonymizer::new();

        fs::create_dir_all(&self.prefix)?;

        let candidate_filename = Path::new("candidates.vcf");
        let mut skips = utils::SimpleCounter::default();

        // get all candidates
        let mut candidates = Vec::new();
        for mut record in self.variants()? {
            let variants =
                utils::collect_variants(&mut record, false, Some(&mut skips), None, None)?;
            for variant in variants {
                candidates.push((variant, record.clone()))
            }
        }

        // get start and end (min/max all start and end pos per chromosome)
        let mut chromosomal_regions: HashMap<Vec<u8>, (u64, u64)> = HashMap::new();
        for (candidate_variant_info, candidate_record) in candidates.clone() {
            let chrom = self
                .candidate_reader()?
                .header()
                .rid2name(candidate_record.rid().unwrap())
                .unwrap()
                .to_owned();
            let rec_pos = candidate_record.pos() as u64;
            let (candidate_start, candidate_end) =
                self.collect_start_end(candidate_variant_info, rec_pos)?;
            chromosomal_regions
                .entry(chrom)
                .and_modify(|e| {
                    e.0 = e.0.min(candidate_start);
                    e.1 = e.1.max(candidate_end);
                })
                .or_insert((candidate_start, candidate_end));
        }

        let extended_chromosomal_regions = self.extend_chromosomal_regions(&chromosomal_regions)?;

        // write bam records
        let mut samples = HashMap::new();
        for (name, path) in &self.bams {
            let properties = sample::estimate_alignment_properties(
                &[path],
                false,
                &mut self.reference_buffer,
                Some(crate::estimation::alignment_properties::NUM_FRAGMENTS),
            )?;

            let filename = Path::new(name).with_extension("bam");

            let mut bam_reader = bam::IndexedReader::from_path(path)?;

            // TODO: create header with just the modified sequence
            // let mut header = bam::header::Header::new();
            // header.push_record(
            //     bam::header::HeaderRecord::new(b"SQ")
            //         .push_tag(b"SN", &str::from_utf8(chrom_name)?)
            //         .push_tag(b"LN", &format!("{}", ref_end - ref_start)),
            // );

            let header = bam::header::Header::from_template(bam_reader.header());

            let mut bam_writer =
                bam::Writer::from_path(self.prefix.join(&filename), &header, bam::Format::Bam)?;
            for (chrom, (start, end)) in chromosomal_regions.clone() {
                let tid: u32 = bam_reader.header().tid(&chrom).unwrap();
                bam_reader.fetch((tid, start, end))?;
                let (ref_start, _) = extended_chromosomal_regions.get(&chrom).unwrap().to_owned();
                for res in bam_reader.records() {
                    let mut rec = res?;
                    // update mapping position to interval
                    rec.set_pos(rec.pos() - ref_start as i64);
                    let mtid = bam_writer.header().tid2name(rec.mtid() as u32);
                    let ref_start_mate = if mtid == b"=" {
                        ref_start
                    } else if let Some(chrom_region) = extended_chromosomal_regions.get(mtid) {
                        chrom_region.0
                    } else {
                        //TODO mate records not being on a candidate chromosome are being ignored by setting offset to 0
                        0
                    };
                    rec.set_mpos(rec.mpos() - ref_start_mate as i64);
                    rec.set_tid(bam_writer.header().tid(&chrom).unwrap() as i32);
                    if rec.remove_aux(b"RG").is_err() {
                        debug!("No RG tag to remove in BAM record.");
                    }
                    if self.anonymize {
                        anonymizer.anonymize_bam_record(&mut rec);
                    }
                    bam_writer.write(&rec)?;
                }
                samples.insert(
                    name.to_owned(),
                    Sample {
                        path: filename.to_str().unwrap().to_owned(),
                        properties: serde_json::to_string(&properties)?,
                        options: self.options.get(name).unwrap().to_owned(),
                    },
                );
            }
        }

        // write candidate
        let mut header = bcf::Header::from_template(self.candidate_reader()?.header());
        if self.anonymize {
            header.remove_generic(b"varlociraptor_preprocess_args");
        }
        let mut candidate_writer = bcf::Writer::from_path(
            self.prefix.join(candidate_filename),
            &header,
            true,
            bcf::Format::Vcf,
        )?;

        // write scenario
        let scenario = if let Some(scenario) = self.scenario.as_ref() {
            fs::copy(scenario, self.prefix.join("scenario.yaml"))?;
            Some("scenario.yaml".to_owned())
        } else {
            None
        };

        let ref_filename = "ref.fa";
        let mut ref_writer = fasta::Writer::to_file(self.prefix.join(ref_filename))?;
        for (_, mut candidate_record) in candidates {
            let chrom = self
                .candidate_reader()?
                .header()
                .rid2name(candidate_record.rid().unwrap())
                .unwrap()
                .to_owned();
            let (ref_start, mut ref_end) =
                extended_chromosomal_regions.get(&chrom).unwrap().to_owned();
            candidate_record.set_pos(candidate_record.pos() - ref_start as i64);
            let mod_alleles =
                modify_bnd_alleles(&candidate_record.alleles(), &extended_chromosomal_regions)?;
            candidate_record.set_alleles(
                &mod_alleles
                    .iter()
                    .map(|inner_vec| inner_vec.as_slice())
                    .collect_vec(),
            )?;
            if let Ok(Some(end)) = candidate_record
                .info(b"END")
                .integer()
                .map(|v| v.map(|v| v[0]))
            {
                candidate_record.push_info_integer(b"END", &[end - ref_start as i32])?;
            }
            if self.anonymize {
                anonymizer.anonymize_bcf_record(&mut candidate_record)?;
            }
            candidate_writer.write(&candidate_record)?;

            // fetch reference
            // limit ref_end
            let ref_name = str::from_utf8(&chrom)?;
            for seq in self.reference_buffer.sequences() {
                if seq.name == ref_name {
                    ref_end = cmp::min(ref_end, seq.len);
                }
            }

            let mut ref_seq = self.reference_buffer.seq(ref_name)?
                [ref_start as usize..ref_end as usize]
                .to_owned();

            if self.anonymize {
                anonymizer.anonymize_seq(&mut ref_seq);
            }

            // write reference

            ref_writer.write(ref_name, None, &ref_seq)?;

            let mut desc = File::create(self.prefix.join("testcase.yaml"))?;
            desc.write_all(
                TestcaseTemplate {
                    samples: samples.clone(),
                    candidate: candidate_filename.to_str().unwrap().to_owned(),
                    ref_path: ref_filename.to_owned(),
                    scenario: scenario.clone(),
                    mode: self.mode,
                    purity: self.purity,
                }
                .render()?
                .as_bytes(),
            )?;
        }

        Ok(())
    }
}
