// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::HashMap;
use std::path::PathBuf;
use std::rc::Rc;
use std::sync::{Arc, Mutex, RwLock};

use anyhow::{bail, Context, Result};
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractLocus};
use progress_logger::ProgressLogger;
use rust_htslib::bam::{self, Read as BAMRead};
use rust_htslib::bcf::{self, Read as BCFRead};

use crate::calling::variants::obs_processing::observation_processor::ObservationProcessor;
use crate::calling::variants::{Call, CallBuilder, VariantBuilder};
use crate::cli;
use crate::errors;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::reference;
use crate::utils::aux_info::AuxInfoCollector;
use crate::utils::collect_variants::VariantInfo;
use crate::utils::variant_buffer::{VariantBuffer, Variants};
use crate::variants;
use crate::variants::evidence::observations::depth_observation::DepthObservation;
use crate::variants::evidence::observations::observation::{AltLocus, ReadPosition, Strand};
use crate::variants::evidence::observations::pileup::Pileup;
use crate::variants::evidence::observations::read_observation::ReadObservationBuilder;
use crate::variants::evidence::realignment::{self, Realignable};
use crate::variants::model::bias::Artifacts;
use crate::variants::model::{self, HaplotypeIdentifier};
use crate::variants::sample::Sample;
use crate::variants::sample::SampleBuilder;
use crate::variants::types::haplotype_block::HaplotypeBlock;
use crate::variants::types::{breakends::Breakend, Loci};

use crate::calling::variants::obs_processing::haplotype_feature_index::HaplotypeFeatureIndex;
use crate::calling::variants::obs_processing::OBSERVATION_FORMAT_VERSION;

#[derive(TypedBuilder)]
pub(crate) struct Preprocessor<R: realignment::Realigner + Clone + 'static> {
    alignment_properties: AlignmentProperties,
    max_depth: usize,
    reference_buffer: Arc<reference::Buffer>,
    realigner: R,
    inbcf: PathBuf,
    outbcf: Option<PathBuf>,
    inbam: PathBuf,
    aux_info_fields: Vec<Vec<u8>>,
    min_bam_refetch_distance: u64,
    options: cli::Varlociraptor,
    haplotype_feature_index: HaplotypeFeatureIndex,
    #[builder(default)]
    breakend_group_builders: RwLock<
        HashMap<Vec<u8>, Mutex<Option<variants::types::breakends::BreakendGroupBuilder<R>>>>,
    >,
    #[builder(default)]
    breakend_groups: RwLock<HashMap<Vec<u8>, Mutex<variants::types::breakends::BreakendGroup<R>>>>,
    #[builder(default)]
    haplotype_blocks: RwLock<
        HashMap<HaplotypeIdentifier, Mutex<variants::types::haplotype_block::HaplotypeBlock>>,
    >,
    log_each_record: bool,
    raw_observation_output: Option<PathBuf>,
    report_fragment_ids: bool,
    adjust_prob_mapping: bool,
    atomic_candidate_variants: bool,
    variant_heterozygosity_field: Option<Vec<u8>>,
    max_number_cn: usize,
    variant_somatic_effective_mutation_rate_field: Option<Vec<u8>>,
}

impl<R: realignment::Realigner + Clone + std::marker::Send + std::marker::Sync> ObservationProcessor
    for Preprocessor<R>
{
    type Realigner = R;

    fn writer(&self, aux_info_collector: &AuxInfoCollector) -> Result<bcf::Writer> {
        let mut header = bcf::Header::new();

        // register tags
        header.push_record(
            b"##INFO=<ID=SVLEN,Number=.,Type=Integer,\
              Description=\"Difference in length between REF and ALT alleles\">",
        );
        header.push_record(
            b"##INFO=<ID=END,Number=1,Type=Integer,\
              Description=\"End position of structural variant (inclusive, 1-based).\">",
        );
        header.push_record(
            b"##INFO=<ID=SVTYPE,Number=1,Type=String,\
              Description=\"Type of structural variant\">",
        );
        header.push_record(
            b"##INFO=<ID=EVENT,Number=1,Type=String,\
              Description=\"ID of event associated to breakend\">",
        );
        header.push_record(
            b"##INFO=<ID=MATEID,Number=.,Type=String,\
              Description=\"ID of mate breakends\">",
        );
        header.push_record(
            b"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">"
        );
        header.push_record(
            b"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">"
        );
        header.push_record(
            b"##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">"
        );
        header.push_record(
            b"##INFO=<ID=HETEROZYGOSITY,Number=A,Type=Float,Description=\"PHRED scaled expected heterozygosity of this particular variant (equivalent to population allele frequency)\">"
        );
        header.push_record(
            b"##INFO=<ID=SOMATIC_EFFECTIVE_MUTATION_RATE,Number=A,Type=Float,Description=\"PHRED scaled expected somatic effective mutation rate of this particular variant (see Williams et al. Nature Genetics 2016)\">"
        );

        // register sequences
        for sequence in self.reference_buffer.sequences() {
            header.push_record(
                format!("##contig=<ID={},length={}>", sequence.name, sequence.len).as_bytes(),
            );
        }

        // store observations
        for name in &vec![
            "FRAGMENT_ID",
            "PROB_MAPPING",
            "PROB_ALT",
            "PROB_REF",
            "PROB_DEPTH",
            "PROB_MISSED_ALLELE",
            "PROB_SAMPLE_ALT",
            "PROB_DOUBLE_OVERLAP",
            "PROB_HIT_BASE",
            "STRAND",
            "READ_ORIENTATION",
            "READ_POSITION",
            "SOFTCLIPPED",
            "ALT_INDEL_OPERATIONS",
            "PAIRED",
            "PROB_HOMOPOLYMER_ARTIFACT_OBSERVABLE",
            "PROB_HOMOPOLYMER_VARIANT_OBSERVABLE",
            "HOMOPOLYMER_INDEL_LEN",
            "IS_MAX_MAPQ",
            "ALT_LOCUS",
            "THIRD_ALLELE_EVIDENCE",
        ] {
            header.push_record(
                format!("##INFO=<ID={},Number=.,Type=Integer,Description=\"Varlociraptor observations (binary encoded, meant for internal use only).\">", name).as_bytes()
            );
        }

        // store options
        header.push_record(
            format!(
                "##varlociraptor_preprocess_args={}",
                serde_json::to_string(&self.options)?
            )
            .as_bytes(),
        );

        // store observation format version
        header.push_record(
            format!(
                "##varlociraptor_observation_format_version={}",
                OBSERVATION_FORMAT_VERSION
            )
            .as_bytes(),
        );

        aux_info_collector.write_header_info(&mut header);

        Ok(if let Some(ref path) = self.outbcf {
            bcf::Writer::from_path(path, &header, false, bcf::Format::Bcf)
                .context(format!("Unable to write BCF to {}.", path.display()))?
        } else {
            bcf::Writer::from_stdout(&header, false, bcf::Format::Bcf)
                .context("Unable to write BCF to STDOUT.")?
        })
    }

    fn process(&mut self) -> Result<()> {
        let mut bcf_reader = bcf::Reader::from_path(&self.inbcf)?;
        bcf_reader.set_threads(1)?;

        let progress_logger = ProgressLogger::builder()
            .with_items_name("records")
            .with_frequency(std::time::Duration::from_secs(20))
            .start();

        let aux_info_collector = AuxInfoCollector::new(&self.aux_info_fields, &bcf_reader)?;

        let mut bcf_writer = self.writer(&aux_info_collector)?;
        bcf_writer.set_threads(1)?;

        let mut variant_buffer = VariantBuffer::new(
            bcf_reader,
            progress_logger,
            self.log_each_record,
            aux_info_collector,
            self.variant_heterozygosity_field.clone(),
            self.variant_somatic_effective_mutation_rate_field.clone(),
        );

        let mut bam_reader =
            bam::IndexedReader::from_path(&self.inbam).context("Unable to read BAM/CRAM file.")?;
        bam_reader.set_threads(1)?;
        bam_reader
            .set_reference(
                self.reference_buffer.reference_path().expect(
                    "bug: reference buffer seemingly has not been created from reference file",
                ),
            )
            .context("Unable to read reference FASTA")?;

        let mut sample = SampleBuilder::default()
            .max_depth(self.max_depth)
            .report_fragment_ids(self.report_fragment_ids)
            .adjust_prob_mapping(self.adjust_prob_mapping)
            .alignments(
                &self.inbam,
                bam_reader,
                self.alignment_properties.clone(),
                self.min_bam_refetch_distance,
            )
            .build()
            .unwrap();

        while let Some(variants) = variant_buffer.next()? {
            let calls = self.process_variant(variants, &mut sample)?;
            for call in calls {
                call.write_preprocessed_record(&mut bcf_writer)?;
            }
        }

        Ok(())
    }

    pub(crate) fn compute_cn(&mut self) -> Result<()> {
        let mut bcf_reader = bcf::Reader::from_path(&self.inbcf)?;
        bcf_reader.set_threads(1)?;

        let progress_logger = ProgressLogger::builder()
            .with_items_name("records")
            .with_frequency(std::time::Duration::from_secs(20))
            .start();

        let aux_info_collector = AuxInfoCollector::new(&self.aux_info_fields, &bcf_reader)?;
        // dbg!(&aux_info_collector);
        // let mut bcf_writer = self.writer(&aux_info_collector)?;

        // Build writer
        let mut header = bcf::Header::from_template(bcf_reader.header());
        header.push_record(
            b"##INFO=<ID=SVLEN,Number=.,Type=Integer,\
              Description=\"Difference in length between REF and ALT alleles\">",
        );
        header.push_record(
            b"##FORMAT=<ID=CN,Number=.,Type=Integer,\
              Description=\"Copy Number\">",
        );
        header.push_record(
            b"##FORMAT=<ID=AF,Number=.,Type=Float,\
              Description=\"Allele Frequency\">",
        );

        let mut bcf_writer = if let Some(ref path) = self.outbcf {
            bcf::Writer::from_path(path, &header, false, bcf::Format::Bcf)
                .context(format!("Unable to write BCF to {}.", path.display()))?
        } else {
            bcf::Writer::from_stdout(&header, false, bcf::Format::Bcf)
                .context("Unable to write BCF to STDOUT.")?
        };
        bcf_writer.set_threads(1)?;

        let mut variant_buffer = VariantBuffer::new(
            bcf_reader,
            progress_logger,
            self.log_each_record,
            aux_info_collector,
            self.variant_heterozygosity_field.clone(),
            self.variant_somatic_effective_mutation_rate_field.clone(),
        );

        let mut bam_reader =
            bam::IndexedReader::from_path(&self.inbam).context("Unable to read BAM/CRAM file.")?;
        bam_reader.set_threads(1)?;
        bam_reader
            .set_reference(
                self.reference_buffer.reference_path().expect(
                    "bug: reference buffer seemingly has not been created from reference file",
                ),
            )
            .context("Unable to read reference FASTA")?;

        let mut sample = SampleBuilder::default()
            .max_depth(self.max_depth)
            .report_fragment_ids(self.report_fragment_ids)
            .adjust_prob_mapping(self.adjust_prob_mapping)
            .alignments(
                &self.inbam,
                bam_reader,
                self.alignment_properties.clone(),
                self.min_bam_refetch_distance,
            )
            .build()
            .unwrap();

        while let Some(variants) = variant_buffer.next()? {
            dbg!(&variants.variant_of_interest());
            let calls = self.process_variant(variants, &mut sample)?;
            for call in calls {
                call.write_final_record(&mut bcf_writer)?;
                // dbg!(&call);
            }
        }

        Ok(())
    }

    fn write_observations(&self, pileup: &Pileup, variants: &Variants) -> Result<()> {
        if let Some(prefix) = &self.raw_observation_output {
            let path = prefix.join(format!(
                "{}_{}_{}.tsv",
                variants.locus().contig(),
                variants.locus().pos(),
                variants.variant_of_interest().variant()
            ));
            if let Some(parent) = path.parent() {
                std::fs::create_dir_all(parent).with_context(|| {
                    format!(
                        "Failed to create directories for writing observations to {}",
                        parent.to_string_lossy()
                    )
                })?;
            }
            let mut wrt = csv::WriterBuilder::new()
                .delimiter(b'\t')
                .from_path(&path)
                .with_context(|| {
                    format!("failed to write observations to {}", path.to_string_lossy())
                })?;

            for obs in pileup.read_observations().iter() {
                wrt.serialize(obs)?;
            }
            for obs in pileup.depth_observations().iter() {
                wrt.serialize(obs)?;
            }
        }

        Ok(())
    }

    fn process_variant(&self, variants: Variants, sample: &mut Sample) -> Result<Vec<Call>> {
        let call_builder = |chrom, start, id| {
            let mut builder = CallBuilder::default();
            builder
                .chrom(chrom)
                .pos(start)
                .id({
                    if id == b"." {
                        None
                    } else {
                        Some(id)
                    }
                })
                .heterozygosity(variants.variant_of_interest().heterozygosity())
                .somatic_effective_mutation_rate(
                    variants
                        .variant_of_interest()
                        .somatic_effective_mutation_rate(),
                );
            builder
        };

        match variants.variant_of_interest() {
            VariantInfo {
                variant,
                haplotype: None,
                ..
            } => {
                let mut call = call_builder(
                    variants.locus().contig().as_bytes().to_owned(),
                    variants.locus().pos(),
                    variants.record_info().id().clone(),
                )
                .aux_info(variants.record_info().aux_info().clone())
                .build()
                .unwrap();

                let chrom_seq = self.reference_buffer.seq(variants.locus().contig())?;
                let pileup = self.process_pileup(&variants, sample)?.unwrap(); // only breakends can lead to None, and they are handled below

                self.write_observations(&pileup, &variants)?;
                match variant {
                    model::Variant::Cnv(_, af) => {
                        // add variant information

                        let sample_info = SampleInfoBuilder::default()
                            .allelefreq_estimate(*af)
                            .pileup(Rc::new(pileup.clone()))
                            .artifacts(Artifacts::none())
                            .vaf_dist(None)
                            .build()
                            .unwrap();

                        // let mut dummy_event_probs: HashMap<String, LogProb> = HashMap::new();
                        // dummy_event_probs.insert("high".to_string(), LogProb::from(0.0)); // 0.0 als Dummy

                        call.variant = Some(
                            VariantBuilder::default()
                                .variant(
                                    variant,
                                    &None,
                                    variants.locus().pos() as usize,
                                    Some(chrom_seq.as_ref()),
                                )
                                .pileup(Some(Rc::new(pileup)))
                                .sample_info(vec![Some(sample_info)])
                                // .event_probs(Some(dummy_event_probs))
                                .build()
                                .unwrap(),
                        );
                    }
                    _ => {
                        // add variant information
                        call.variant = Some(
                            VariantBuilder::default()
                                .variant(
                                    variant,
                                    &None,
                                    variants.locus().pos() as usize,
                                    Some(chrom_seq.as_ref()),
                                )
                                .pileup(Some(Rc::new(pileup)))
                                .build()
                                .unwrap(),
                        );
                    }
                }

                Ok(vec![call])
            }
            VariantInfo {
                haplotype: Some(haplotype),
                ..
            } => {
                let mut calls = Vec::new();
                // process breakend
                match variants.variant_of_interest().variant() {
                    model::Variant::Breakend { .. } => {
                        match haplotype {
                            HaplotypeIdentifier::Event(event) => {
                                if let Some(pileup) = self.process_pileup(&variants, sample)? {
                                    self.write_observations(&pileup, &variants)?;

                                    let pileup = Rc::new(pileup);
                                    for breakend in self
                                        .breakend_groups
                                        .read()
                                        .unwrap()
                                        .get(event)
                                        .as_ref()
                                        .unwrap()
                                        .lock()
                                        .unwrap()
                                        .breakends()
                                    {
                                        let mut call_builder = call_builder(
                                            breakend.locus().contig().as_bytes().to_owned(),
                                            breakend.locus().pos(),
                                            breakend.id().to_owned(),
                                        );
                                        call_builder.mateid(breakend.mateid().to_owned());
                                        if let Some(ref aux_info) = breakend.aux_info() {
                                            call_builder.aux_info(aux_info.clone());
                                        }
                                        let mut call = call_builder.build().unwrap();

                                        // add variant information
                                        call.variant = Some(
                                            VariantBuilder::default()
                                                .variant(
                                                    &breakend.to_variant(),
                                                    &Some(haplotype.to_owned()),
                                                    breakend.locus().pos() as usize,
                                                    None,
                                                )
                                                .precision(breakend.precision().to_owned())
                                                .pileup(Some(Rc::clone(&pileup)))
                                                .build()
                                                .unwrap(),
                                        );
                                        calls.push(call);
                                    }
                                    // As all records a written, the breakend group can be discarded.
                                    self.breakend_groups.write().unwrap().remove(event);
                                }
                            }
                        }
                    }
                    _ => {
                        // process haplotype block
                        // if this is the last variant in the block
                        if let Some(pileup) = self.process_pileup(&variants, sample)? {
                            self.write_observations(&pileup, &variants)?;

                            let pileup = Rc::new(pileup);
                            {
                                let haplotype_blocks = self.haplotype_blocks.read().unwrap();
                                let haplotype_block = haplotype_blocks
                                    .get(haplotype)
                                    .as_ref()
                                    .unwrap()
                                    .lock()
                                    .unwrap();

                                fn to_call<L: Loci>(
                                    loci: &L,
                                    variant: model::Variant,
                                    reference_buffer: Arc<reference::Buffer>,
                                    haplotype: &HaplotypeIdentifier,
                                    pileup: Rc<Pileup>,
                                ) -> Result<Call> {
                                    if let Some(contig) = loci.contig() {
                                        let mut call = CallBuilder::default()
                                            .chrom(loci.contig().unwrap().as_bytes().to_owned())
                                            .pos(loci.first_pos())
                                            .build()
                                            .unwrap();

                                        let chrom_seq = reference_buffer.seq(contig)?;

                                        // add variant information
                                        call.variant = Some(
                                            VariantBuilder::default()
                                                .variant(
                                                    &variant,
                                                    &Some(haplotype.to_owned()),
                                                    loci.first_pos() as usize,
                                                    Some(chrom_seq.as_ref()),
                                                )
                                                .pileup(Some(Rc::clone(&pileup)))
                                                .build()
                                                .unwrap(),
                                        );
                                        Ok(call)
                                    } else {
                                        unreachable!("bug: unexpected multi-contig variant (must be breakend, which is unsupported for now) in haplotype block");
                                    }
                                }

                                for variant in haplotype_block.variants() {
                                    calls.push(to_call(
                                        variant.loci(),
                                        variant.to_variant_representation(),
                                        Arc::clone(&self.reference_buffer),
                                        haplotype,
                                        Rc::clone(&pileup),
                                    )?);
                                }
                            }
                            // As all records a written, the haplotype block can be discarded.
                            self.haplotype_blocks.write().unwrap().remove(haplotype);
                        }
                    }
                }
                Ok(calls)
            }
        }
    }

    fn process_pileup(&self, variants: &Variants, sample: &mut Sample) -> Result<Option<Pileup>> {
        let interval = |len: u64, offset: u64| {
            genome::Interval::new(
                variants.locus().contig().to_owned(),
                variants.locus().pos() - offset..variants.locus().pos() - offset + len,
            )
        };

        let ref_base = || {
            self.reference_buffer
                .seq(variants.locus().contig())?
                .get(variants.locus().pos() as usize)
                .cloned()
                .ok_or_else(|| -> anyhow::Error {
                    errors::invalid_bcf_record(
                        variants.locus().contig(),
                        variants.locus().pos() as i64,
                        "position larger than reference length",
                    )
                    .into()
                })
        };

        let parse_snv = |alt| -> Result<variants::types::Snv<R>> {
            let locus = variants.locus().clone();
            Ok(variants::types::Snv::new(
                locus,
                ref_base()?,
                alt,
                self.realigner.clone(),
                !self.atomic_candidate_variants,
            ))
        };

        let parse_mnv = |alt: &Vec<u8>| -> Result<variants::types::Mnv<R>> {
            Ok(variants::types::Mnv::new(
                variants.locus().clone(),
                self.reference_buffer
                    .seq(variants.locus().contig())?
                    .get(
                        variants.locus().pos() as usize
                            ..variants.locus().pos() as usize + alt.len(),
                    )
                    .ok_or_else(|| -> anyhow::Error {
                        errors::invalid_bcf_record(
                            variants.locus().contig(),
                            variants.locus().pos() as i64,
                            "MNV exceeds reference length",
                        )
                        .into()
                    })?
                    .to_owned(),
                alt.to_owned(),
                self.realigner.clone(),
                !self.atomic_candidate_variants,
            ))
        };

        let parse_none = || -> Result<variants::types::None> {
            Ok(variants::types::None::new(
                variants.locus().clone(),
                ref_base()?,
            ))
        };

        let parse_deletion =
            |len| variants::types::Deletion::new(interval(len, 0), self.realigner.clone());

        let cnv_alt =
            |len| variants::types::Deletion::new(interval(len, 1), self.realigner.clone());

        let parse_insertion = |seq: &Vec<u8>| {
            variants::types::Insertion::new(
                variants.locus().clone(),
                seq.to_owned(),
                self.realigner.clone(),
            )
        };

        let parse_inversion = |len| -> Result<variants::types::Inversion<R>> {
            Ok(variants::types::Inversion::new(
                interval(len, 0),
                self.realigner.clone(),
                self.reference_buffer
                    .seq(variants.locus().contig())?
                    .as_ref(),
            ))
        };

        let parse_cnv = |len, af| -> Result<variants::types::Cnv<R>> {
            Ok(variants::types::Cnv::new(
                interval(len, 0),
                self.realigner.clone(),
                self.reference_buffer
                    .seq(variants.locus().contig())?
                    .as_ref(),
                af,
            ))
        };

        let parse_duplication = |len| -> Result<variants::types::Duplication<R>> {
            Ok(variants::types::Duplication::new(
                interval(len, 0),
                self.realigner.clone(),
                self.reference_buffer
                    .seq(variants.locus().contig())?
                    .as_ref(),
            ))
        };

        let parse_replacement = |ref_allele: &Vec<u8>,
                                 alt_allele: &Vec<u8>|
         -> Result<variants::types::Replacement<R>> {
            variants::types::Replacement::new(
                interval(ref_allele.len() as u64, 0),
                alt_allele.to_owned(),
                self.realigner.clone(),
            )
        };

        let mut alt_variants = variants
            .alt_variants()
            .filter(|variant_info| {
                !variant_info.variant().is_breakend() && !variant_info.variant().is_none()
            })
            .map(|variant_info| -> Result<Box<dyn Realignable>> {
                Ok(match variant_info.variant() {
                    model::Variant::Snv(alt) => Box::new(parse_snv(*alt)?),
                    model::Variant::Mnv(alt) => Box::new(parse_mnv(alt)?),
                    model::Variant::Deletion(l) => Box::new(parse_deletion(*l)?),
                    model::Variant::Insertion(seq) => Box::new(parse_insertion(seq)?),
                    model::Variant::Inversion(len) => Box::new(parse_inversion(*len)?),
                    model::Variant::Cnv(len, af) => Box::new(parse_cnv(*len, *af)?),
                    model::Variant::Duplication(len) => Box::new(parse_duplication(*len)?),
                    model::Variant::Replacement {
                        ref_allele,
                        alt_allele,
                    } => Box::new(parse_replacement(ref_allele, alt_allele)?),
                    model::Variant::Breakend { .. } => {
                        unreachable!();
                    }
                    model::Variant::None => {
                        unreachable!();
                    }
                })
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Some(
            if let Some(ref haplotype) = variants.variant_of_interest().haplotype() {
                match variants.variant_of_interest().variant() {
                    model::Variant::Breakend {
                        ref_allele,
                        spec,
                        precision,
                    } => {
                        let HaplotypeIdentifier::Event(event) = haplotype;
                        {
                            if !self
                                .breakend_group_builders
                                .read()
                                .unwrap()
                                .contains_key(event)
                            {
                                let mut builder =
                                    variants::types::breakends::BreakendGroupBuilder::new();
                                builder.realigner(self.realigner.clone());
                                self.breakend_group_builders
                                    .write()
                                    .unwrap()
                                    .insert(event.to_owned(), Mutex::new(Some(builder)));
                            }
                        }
                        let group_builders = self.breakend_group_builders.read().unwrap();

                        let mut group = group_builders.get(event).unwrap().lock().unwrap();

                        if let Some(group) = group.as_mut() {
                            if let Some(breakend) = Breakend::new(
                                variants.locus().clone(),
                                ref_allele,
                                spec,
                                variants.record_info().id(),
                                variants.record_info().mateid().clone(),
                                precision.clone(),
                                variants.record_info().aux_info().clone(),
                            )? {
                                group.push_breakend(breakend);

                                // If this is the last breakend in the group, process!
                                if self
                                    .haplotype_feature_index
                                    .last_record_index(haplotype)
                                    .unwrap()
                                    == variants.record_info().index()
                                {
                                    // METHOD: last record of the breakend event. Hence, we can extract observations.
                                    if let Some(group) = group.build() {
                                        let breakend_group = Mutex::new(group);
                                        self.breakend_groups
                                            .write()
                                            .unwrap()
                                            .insert(event.to_owned(), breakend_group);
                                        sample.extract_read_observations(
                                            &*self
                                                .breakend_groups
                                                .read()
                                                .unwrap()
                                                .get(event)
                                                .unwrap()
                                                .lock()
                                                .unwrap(),
                                            &Vec::new(), // Do not consider alt variants in case of breakends
                                        )?
                                    } else {
                                        // skipped with message in the builder
                                        return Ok(None);
                                    }
                                } else {
                                    // Not yet the last one, skip.
                                    return Ok(None);
                                }
                            } else {
                                // Breakend type not supported, remove breakend group.
                                self.breakend_group_builders
                                    .write()
                                    .unwrap()
                                    .insert(event.to_owned(), Mutex::new(None));
                                return Ok(None);
                            }
                        } else {
                            // Breakend group has been removed before because one breakend was invalid.
                            return Ok(None);
                        }
                    }
                    variant => {
                        // other variants are grouped into haplotype blocks
                        if !self
                            .haplotype_blocks
                            .read()
                            .unwrap()
                            .contains_key(haplotype)
                        {
                            // add new haplotype block
                            self.haplotype_blocks
                                .write()
                                .unwrap()
                                .insert(haplotype.to_owned(), Mutex::new(HaplotypeBlock::new()));
                        }
                        let haplotype_blocks = self.haplotype_blocks.read().unwrap();

                        let mut haplotype_block =
                            haplotype_blocks.get(haplotype).unwrap().lock().unwrap();

                        match variant {
                            model::Variant::Snv(alt) => {
                                haplotype_block.push_variant(Box::new(parse_snv(*alt)?))
                            }
                            model::Variant::Mnv(alt) => {
                                haplotype_block.push_variant(Box::new(parse_mnv(alt)?))
                            }
                            model::Variant::None => {
                                haplotype_block.push_variant(Box::new(parse_none()?))
                            }
                            model::Variant::Deletion(l) => {
                                haplotype_block.push_variant(Box::new(parse_deletion(*l)?))
                            }
                            model::Variant::Insertion(seq) => {
                                haplotype_block.push_variant(Box::new(parse_insertion(seq)?))
                            }
                            model::Variant::Inversion(len) => {
                                haplotype_block.push_variant(Box::new(parse_inversion(*len)?))
                            }
                            model::Variant::Cnv(len, af) => {
                                haplotype_block.push_variant(Box::new(parse_cnv(*len, *af)?))
                            }
                            model::Variant::Duplication(len) => {
                                haplotype_block.push_variant(Box::new(parse_duplication(*len)?))
                            }
                            model::Variant::Replacement {
                                ref_allele,
                                alt_allele,
                            } => haplotype_block
                                .push_variant(Box::new(parse_replacement(ref_allele, alt_allele)?)),
                            model::Variant::Breakend { .. } => {
                                bail!(errors::Error::HaplotypeBlockWithBreakend);
                            }
                        }

                        // If this is the last variant in the group, process!
                        if self
                            .haplotype_feature_index
                            .last_record_index(haplotype)
                            .unwrap()
                            == variants.record_info().index()
                        {
                            sample.extract_read_observations(&*haplotype_block, &Vec::new())?
                        } else {
                            return Ok(None);
                        }
                    }
                }
            } else {
                // single variants
                match variants.variant_of_interest().variant() {
                    model::Variant::Snv(alt) => {
                        sample.extract_read_observations(&parse_snv(*alt)?, &Vec::new())?
                    }
                    model::Variant::Mnv(alt) => {
                        sample.extract_read_observations(&parse_mnv(alt)?, &alt_variants)?
                    }
                    model::Variant::None => {
                        sample.extract_read_observations(&parse_none()?, &alt_variants)?
                    }
                    model::Variant::Deletion(l) => {
                        sample.extract_read_observations(&parse_deletion(*l)?, &alt_variants)?
                    }
                    model::Variant::Insertion(seq) => {
                        sample.extract_read_observations(&parse_insertion(seq)?, &alt_variants)?
                    }
                    model::Variant::Inversion(len) => {
                        sample.extract_read_observations(&parse_inversion(*len)?, &alt_variants)?
                    }
                    model::Variant::Cnv(len, af) => {
                        // Add deletion of the cn as an alternative variant.
                        // alt_variants.push(Box::new(cnv_alt(*len, *af)?));
                        sample.extract_depth_observations(
                            &parse_cnv(*len, *af)?,
                            &alt_variants,
                            self.max_number_cn,
                        )?
                        // sample.reset_record_buffer();

                        // // Determine if CNV should be treated as a deletion based on depth probabilities.
                        // if let Some((max_idx, _)) =
                        //     pileup.depth_observations().first().and_then(|obs| {
                        //         obs.cnv_probs()
                        //             .iter()
                        //             .enumerate()
                        //             .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                        //     })
                        // {
                        //     let ploidy = 2; // TODO: make dynamic
                        //     if max_idx < ploidy {
                        //         pileup = sample.extract_read_observations(
                        //             &parse_deletion(*len)?,
                        //             &alt_variants,
                        //         )?;
                        //     }
                        // }

                        // pileup
                    }
                    model::Variant::Duplication(len) => sample
                        .extract_read_observations(&parse_duplication(*len)?, &alt_variants)?,
                    model::Variant::Replacement {
                        ref_allele,
                        alt_allele,
                    } => sample.extract_read_observations(
                        &parse_replacement(ref_allele, alt_allele)?,
                        &alt_variants,
                    )?,
                    model::Variant::Breakend { .. } => unimplemented!(
                        "bug: breakends without haplotype events should be ignored for now"
                    ),
                }
            },
        ))
    }
}
