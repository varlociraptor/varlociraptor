use std::cmp;
use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::str;

use anyhow::Result;
use askama::Template;
use bio::io::fasta;
use bio_types::genome::AbstractLocus;
use derive_builder::Builder;
use regex::Regex;
use rust_htslib::bam::Read as BamRead;
use rust_htslib::{bam, bcf, bcf::Read};

use crate::cli;
use crate::errors;
use crate::utils;
use crate::variants::model::Variant;
use crate::variants::sample;
use crate::variants::types::breakends::BreakendIndex;

lazy_static! {
    static ref TESTCASE_RE: Regex =
        Regex::new(r"^(?P<chrom>[^:]+):(?P<pos>\d+)(:(?P<idx>\d+))?$").unwrap();
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

#[derive(Debug)]
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
    reference_reader: fasta::IndexedReader<File>,
    candidates: PathBuf,
    #[builder(private)]
    bams: HashMap<String, PathBuf>,
    scenario: Option<PathBuf>,
    #[builder(private)]
    options: HashMap<String, String>,
    #[builder(default = "None")]
    purity: Option<f64>,
    mode: Mode,
}

impl TestcaseBuilder {
    pub(crate) fn reference(self, path: impl AsRef<Path>) -> Result<Self> {
        Ok(self.reference_reader(fasta::IndexedReader::from_file(&path)?))
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
        options: &cli::Varlociraptor,
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

        self.options
            .as_mut()
            .unwrap()
            .insert(name.to_owned(), serde_json::to_string(options)?);

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
                        if let Some(event) =
                            utils::info_tag_event(rec)?.map(|event| event.to_owned())
                        {
                            // METHOD: for breakend events, collect all the other breakends.
                            if breakend_index.is_none() {
                                breakend_index = Some(BreakendIndex::new(&self.candidates)?);
                            }
                            let breakend_index = breakend_index.as_ref().unwrap();
                            let last_idx = breakend_index.last_record_index(&event).unwrap();

                            let mut candidate_reader = self.candidate_reader()?;
                            for (i, res) in candidate_reader.records().enumerate() {
                                let mut other_rec = res?;
                                if let Some(other_event) = utils::info_tag_event(&mut other_rec)? {
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

    pub(crate) fn write(&mut self) -> Result<()> {
        fs::create_dir_all(&self.prefix)?;

        let candidate_filename = Path::new("candidates.vcf");

        // get and write candidate
        let mut candidate = None;
        for (i, mut record) in (self.variants()?).into_iter().enumerate() {
            let variants = utils::collect_variants(&mut record)?;
            for variant in variants {
                if i == self.idx {
                    // if no chromosome was specified, we infer the locus from the matching
                    // variant
                    if self.chrom_name.is_none() {
                        self.chrom_name = Some(
                            self.candidate_reader()?
                                .header()
                                .rid2name(record.rid().unwrap())
                                .unwrap()
                                .to_owned(),
                        );
                        self.pos = Some(record.pos() as u64);
                    }

                    candidate = Some((variant, record));

                    break;
                }
            }
        }
        if candidate.is_none() {
            return Err(errors::Error::InvalidIndex.into());
        }
        let candidate = candidate.unwrap();

        let chrom_name = self.chrom_name.as_ref().unwrap();
        let pos = self.pos.unwrap();

        let (start, end) = match candidate {
            (Variant::Deletion(l), _) => (pos.saturating_sub(1000), pos + l as u64 + 1000),
            (Variant::Insertion(ref seq), _) => {
                (pos.saturating_sub(1000), pos + seq.len() as u64 + 1000)
            }
            (Variant::SNV(_), _) => (pos.saturating_sub(100), pos + 1 + 100),
            (Variant::MNV(ref bases), _) => {
                (pos.saturating_sub(100), pos + bases.len() as u64 + 100)
            }
            (Variant::Breakend { .. }, _) => {
                (pos.saturating_sub(1000), pos + 1 + 1000) // TODO collect entire breakend event!
            }
            (Variant::Inversion(l), _) => (pos.saturating_sub(1000), pos + l as u64 + 1000),
            (Variant::Duplication(l), _) => (pos.saturating_sub(1000), pos + l as u64 + 1000),
            (Variant::None, _) => (pos.saturating_sub(100), pos + 1 + 100),
        };

        let mut ref_start = start;
        let mut ref_end = end;
        // first pass, extend reference interval
        for path in self.bams.values() {
            let mut bam_reader = bam::IndexedReader::from_path(path)?;

            let tid = bam_reader.header().tid(chrom_name).unwrap();
            bam_reader.fetch(tid, start, end)?;
            for res in bam_reader.records() {
                let rec = res?;
                let seq_len = rec.seq().len() as u64;
                ref_start = cmp::min((rec.pos() as u64).saturating_sub(seq_len), ref_start);
                ref_end = cmp::max(rec.cigar().end_pos() as u64 + seq_len, ref_end);
            }
        }

        // second pass, write samples
        let mut samples = HashMap::new();
        for (name, path) in &self.bams {
            let properties = sample::estimate_alignment_properties(path, false)?;
            let mut bam_reader = bam::IndexedReader::from_path(path)?;
            let filename = Path::new(name).with_extension("bam");
            let mut bam_writer = bam::Writer::from_path(
                self.prefix.join(&filename),
                &bam::Header::from_template(bam_reader.header()),
                bam::Format::BAM,
            )?;

            let tid = bam_reader.header().tid(chrom_name).unwrap();

            bam_reader.fetch(tid, start, end)?;
            for res in bam_reader.records() {
                let mut rec = res?;
                // update mapping position to interval
                rec.set_pos(rec.pos() - ref_start as i64);
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

        // write candidate
        let mut candidate_writer = bcf::Writer::from_path(
            self.prefix.join(candidate_filename),
            &bcf::Header::from_template(self.candidate_reader()?.header()),
            true,
            bcf::Format::BCF,
        )?;
        let (_, mut candidate_record) = candidate;
        candidate_record.set_pos(candidate_record.pos() - ref_start as i64);
        if let Ok(Some(end)) = candidate_record
            .info(b"END")
            .integer()
            .map(|v| v.map(|v| v[0]))
        {
            candidate_record.push_info_integer(b"END", &[end - ref_start as i32])?;
        }
        candidate_writer.write(&candidate_record)?;

        // write scenario
        let scenario = if let Some(scenario) = self.scenario.as_ref() {
            fs::copy(scenario, self.prefix.join("scenario.yaml"))?;
            Some("scenario.yaml".to_owned())
        } else {
            None
        };

        // fetch reference
        let ref_name = str::from_utf8(&chrom_name)?;
        self.reference_reader
            .fetch(ref_name, ref_start as u64, ref_end as u64)?;
        let mut ref_seq = Vec::new();
        self.reference_reader.read(&mut ref_seq)?;

        // write reference
        let ref_filename = "ref.fa";
        let mut ref_writer = fasta::Writer::to_file(self.prefix.join(ref_filename))?;
        ref_writer.write(ref_name, None, &ref_seq)?;

        let mut desc = File::create(self.prefix.join("testcase.yaml"))?;
        desc.write_all(
            TestcaseTemplate {
                samples,
                candidate: candidate_filename.to_str().unwrap().to_owned(),
                ref_path: ref_filename.to_owned(),
                scenario,
                mode: self.mode,
                purity: self.purity,
            }
            .render()?
            .as_bytes(),
        )?;

        Ok(())
    }
}
