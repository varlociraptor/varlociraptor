// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::HashMap;
use std::fmt::Debug;
use std::fs;
use std::path::{Path, PathBuf};
use std::str;
use std::sync::{Arc, Mutex, RwLock};

use anyhow::{Context, Result};
use bio::io::fasta;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractLocus};
use bv::BitVec;
use byteorder::{ByteOrder, LittleEndian};
use crossbeam::channel::{Receiver, Sender};
use derive_builder::Builder;
use itertools::Itertools;
use rust_htslib::bam;
use rust_htslib::bcf::{self, Read};

use crate::calling::variants::{chrom, Call, CallBuilder, VariantBuilder};
use crate::cli;
use crate::errors;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::reference;
use crate::utils;
use crate::utils::worker_pool;
use crate::utils::worker_pool::Orderable;
use crate::utils::MiniLogProb;
use crate::variants;
use crate::variants::evidence::observation::{Observation, ObservationBuilder};
use crate::variants::evidence::realignment;
use crate::variants::model;
use crate::variants::sample::Sample;
use crate::variants::sample::{ProtocolStrandedness, SampleBuilder};
use crate::variants::types::breakends::{Breakend, BreakendIndex};

#[derive(Builder)]
#[builder(pattern = "owned")]
pub(crate) struct ObservationProcessor {
    threads: usize,
    buffer_capacity: usize,
    alignment_properties: AlignmentProperties,
    max_depth: usize,
    protocol_strandedness: ProtocolStrandedness,
    #[builder(private)]
    reference_buffer: Arc<reference::Buffer>,
    #[builder(private)]
    realigner: realignment::Realigner,
    inbcf: PathBuf,
    outbcf: Option<PathBuf>,
    inbam: PathBuf,
    options: cli::Varlociraptor,
    breakend_index: BreakendIndex,
    #[builder(default)]
    breakend_group_builders:
        RwLock<HashMap<Vec<u8>, Mutex<Option<variants::types::breakends::BreakendGroupBuilder>>>>,
    #[builder(default)]
    breakend_groups: RwLock<HashMap<Vec<u8>, Mutex<variants::types::breakends::BreakendGroup>>>,
}

impl ObservationProcessorBuilder {
    pub(crate) fn reference(self, reader: fasta::IndexedReader<fs::File>) -> Self {
        self.reference_buffer(Arc::new(reference::Buffer::new(reader, 3)))
    }

    pub(crate) fn realignment(
        self,
        gap_params: realignment::pairhmm::GapParams,
        window: u64,
    ) -> Self {
        let ref_buffer = Arc::clone(
            self.reference_buffer
                .as_ref()
                .expect("You need to set reference before setting realignment parameters"),
        );

        self.realigner(realignment::Realigner::new(ref_buffer, gap_params, window))
    }
}

impl ObservationProcessor {
    fn writer(&self) -> Result<bcf::Writer> {
        let mut header = bcf::Header::new();

        // register tags
        header.push_record(
            b"##INFO=<ID=SVLEN,Number=A,Type=Integer,\
              Description=\"Difference in length between REF and ALT alleles\">",
        );
        header.push_record(
            b"##INFO=<ID=END,Number=A,Type=Integer,\
              Description=\"End position of structural variant (inclusive, 1-based).\">",
        );
        header.push_record(
            b"##INFO=<ID=SVTYPE,Number=A,Type=String,\
              Description=\"Structural variant type\">",
        );
        header.push_record(
            b"##INFO=<ID=EVENT,Number=A,Type=String,\
              Description=\"ID of event associated to breakend\">",
        );
        header.push_record(
            b"##INFO=<ID=MATEID,Number=1,Type=String,\
              Description=\"ID of mate breakend\">",
        );

        // register sequences
        for sequence in self.reference_buffer.sequences() {
            header.push_record(
                format!("##contig=<ID={},length={}>", sequence.name, sequence.len).as_bytes(),
            );
        }

        // store observations
        for name in &vec![
            "PROB_MAPPING",
            "PROB_ALT",
            "PROB_REF",
            "PROB_MISSED_ALLELE",
            "PROB_SAMPLE_ALT",
            "PROB_DOUBLE_OVERLAP",
            "PROB_ANY_STRAND",
            "FORWARD_STRAND",
            "REVERSE_STRAND",
        ] {
            header.push_record(
                format!("##INFO=<ID={},Number=.,Type=Integer,Description=\"Varlociraptor observations (binary encoded, meant internal use only).\"", name).as_bytes()
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

        Ok(if let Some(ref path) = self.outbcf {
            bcf::Writer::from_path(path, &header, false, bcf::Format::BCF)
                .context(format!("Unable to write BCF to {}.", path.display()))?
        } else {
            bcf::Writer::from_stdout(&header, false, bcf::Format::BCF)
                .context("Unable to write BCF to STDOUT.")?
        })
    }

    pub(crate) fn process(&mut self) -> Result<()> {
        let mut workers = Vec::new();

        for _ in 0..self.threads {
            let worker = |receiver: Receiver<WorkItem>, sender: Sender<Box<Calls>>| -> Result<()> {
                let bam_reader = bam::IndexedReader::from_path(&self.inbam)
                    .context("Unable to read BAM/CRAM file.")?;

                let mut sample = SampleBuilder::default()
                    .max_depth(self.max_depth)
                    .protocol_strandedness(self.protocol_strandedness)
                    .alignments(bam_reader, self.alignment_properties)
                    .build()
                    .unwrap();
                for work_item in receiver {
                    let calls = self.process_record(work_item, &mut sample)?;
                    sender.send(calls).unwrap();
                }
                Ok(())
            };
            workers.push(worker);
        }

        let postprocessor = |receiver: Receiver<Box<Calls>>| -> Result<()> {
            let mut bcf_writer = self.writer()?;
            for calls in receiver {
                for call in calls.iter() {
                    call.write_preprocessed_record(&mut bcf_writer)?;

                    if calls.index() % 100 == 0 {
                        info!("{} records processed.", calls.index());
                    }
                }
            }

            Ok(())
        };

        let preprocessor = |sender: Sender<WorkItem>,
                            buffer_guard: Arc<utils::worker_pool::BufferGuard>|
         -> Result<()> {
            let mut bcf_reader = bcf::Reader::from_path(&self.inbcf)?;

            let mut i = 0;
            loop {
                let mut record = bcf_reader.empty_record();
                if !bcf_reader.read(&mut record)? {
                    return Ok(());
                }

                // process record
                let work_item = WorkItem {
                    start: record.pos() as u64,
                    chrom: String::from_utf8(chrom(&bcf_reader, &record).to_owned()).unwrap(),
                    variants: utils::collect_variants(&mut record)?,
                    record_id: record.id(),
                    record_mateid: utils::info_tag_mateid(&mut record)
                        .map_or(None, |mateid| mateid.map(|mateid| mateid.to_owned())),
                    record_index: i,
                };

                sender.send(work_item).unwrap();

                i += 1;

                buffer_guard.wait_for_free();
            }
        };

        worker_pool(
            preprocessor,
            workers.into_iter(),
            postprocessor,
            self.buffer_capacity,
        )
    }

    fn process_record(&self, work_item: WorkItem, sample: &mut Sample) -> Result<Box<Calls>> {
        if work_item.variants.is_empty() {
            return Ok(Box::new(Calls::new(work_item.record_index, vec![])));
        }

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
                .variants(Vec::new());
            builder
        };

        if work_item
            .variants
            .iter()
            .all(|variant| !variant.is_breakend())
        {
            let mut call = call_builder(
                work_item.chrom.as_bytes().to_owned(),
                work_item.start,
                work_item.record_id.clone(),
            )
            .build()
            .unwrap();

            for variant in work_item
                .variants
                .iter()
                .filter(|variant| !variant.is_breakend())
            {
                let chrom_seq = self.reference_buffer.seq(&work_item.chrom)?;
                let pileup = self.process_variant(&variant, &work_item, sample)?.unwrap(); // only breakends can lead to None, and they are handled below

                // add variant information
                call.variants.push(
                    VariantBuilder::default()
                        .variant(&variant, work_item.start as usize, Some(chrom_seq.as_ref()))
                        .observations(Some(pileup))
                        .build()
                        .unwrap(),
                );
            }

            Ok(Box::new(Calls::new(work_item.record_index, vec![call])))
        } else {
            let mut calls = Vec::new();
            for variant in work_item.variants.iter() {
                if let model::Variant::Breakend { event, .. } = variant {
                    if let Some(pileup) = self.process_variant(variant, &work_item, sample)? {
                        let mut pileup = Some(pileup);
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
                            let mut call = call_builder(
                                breakend.locus().contig().as_bytes().to_owned(),
                                breakend.locus().pos(),
                                breakend.id().to_owned(),
                            )
                            .mateid(breakend.mateid().to_owned())
                            .build()
                            .unwrap();

                            // add variant information
                            call.variants.push(
                                VariantBuilder::default()
                                    .variant(
                                        &breakend.to_variant(event),
                                        breakend.locus().pos() as usize,
                                        None,
                                    )
                                    .observations(pileup)
                                    .build()
                                    .unwrap(),
                            );
                            calls.push(call);
                            // Subsequent calls should not contain the pileup again (saving space).
                            pileup = None;
                        }
                        // As all records a written, the breakend group can be discarded.
                        self.breakend_groups.write().unwrap().remove(event);
                    }
                }
            }
            Ok(Box::new(Calls::new(work_item.record_index, calls)))
        }
    }

    fn process_variant(
        &self,
        variant: &model::Variant,
        work_item: &WorkItem,
        sample: &mut Sample,
    ) -> Result<Option<Vec<Observation>>> {
        let locus = || genome::Locus::new(work_item.chrom.clone(), work_item.start);
        let interval = |len: u64| {
            genome::Interval::new(
                work_item.chrom.clone(),
                work_item.start..work_item.start + len,
            )
        };
        let start = work_item.start as usize;

        Ok(Some(match variant {
            model::Variant::SNV(alt) => sample.extract_observations(&variants::types::SNV::new(
                locus(),
                self.reference_buffer.seq(&work_item.chrom)?[start],
                *alt,
            ))?,
            model::Variant::MNV(alt) => sample.extract_observations(&variants::types::MNV::new(
                locus(),
                self.reference_buffer.seq(&work_item.chrom)?[start..start + alt.len()].to_owned(),
                alt.to_owned(),
            ))?,
            model::Variant::None => sample.extract_observations(&variants::types::None::new(
                locus(),
                self.reference_buffer.seq(&work_item.chrom)?[start],
            ))?,
            model::Variant::Deletion(l) => sample.extract_observations(
                &variants::types::Deletion::new(interval(*l), self.realigner.clone()),
            )?,
            model::Variant::Insertion(seq) => sample.extract_observations(
                &variants::types::Insertion::new(locus(), seq.to_owned(), self.realigner.clone()),
            )?,
            model::Variant::Inversion(len) => {
                sample.extract_observations(&variants::types::Inversion::new(
                    interval(*len),
                    self.realigner.clone(),
                    self.reference_buffer.seq(&work_item.chrom)?.as_ref(),
                ))?
            }
            model::Variant::Duplication(len) => {
                sample.extract_observations(&variants::types::Duplication::new(
                    interval(*len),
                    self.realigner.clone(),
                    self.reference_buffer.seq(&work_item.chrom)?.as_ref(),
                ))?
            }
            model::Variant::Replacement {
                ref_allele,
                alt_allele,
            } => sample.extract_observations(&variants::types::Replacement::new(
                interval(ref_allele.len() as u64),
                alt_allele.to_owned(),
                self.realigner.clone(),
                self.reference_buffer.seq(&work_item.chrom)?.as_ref(),
            ))?,
            model::Variant::Breakend {
                ref_allele,
                spec,
                event,
            } => {
                {
                    if !self
                        .breakend_group_builders
                        .read()
                        .unwrap()
                        .contains_key(event)
                    {
                        let mut builder =
                            variants::types::breakends::BreakendGroupBuilder::default();
                        builder.set_realigner(self.realigner.clone());
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
                        locus(),
                        ref_allele,
                        spec,
                        &work_item.record_id,
                        work_item.record_mateid.clone(),
                    )? {
                        group.push_breakend(breakend);

                        if self.breakend_index.last_record_index(event).unwrap()
                            == work_item.record_index
                        {
                            // METHOD: last record of the breakend event. Hence, we can extract observations.
                            let breakend_group = Mutex::new(group.build().unwrap());
                            self.breakend_groups
                                .write()
                                .unwrap()
                                .insert(event.to_owned(), breakend_group);
                            sample.extract_observations(
                                &*self
                                    .breakend_groups
                                    .read()
                                    .unwrap()
                                    .get(event)
                                    .unwrap()
                                    .lock()
                                    .unwrap(),
                            )?
                        } else {
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
        }))
    }
}

pub(crate) static OBSERVATION_FORMAT_VERSION: &str = "2";

/// Read observations from BCF record.
pub(crate) fn read_observations(record: &mut bcf::Record) -> Result<Vec<Observation>> {
    fn read_values<T>(record: &mut bcf::Record, tag: &[u8]) -> Result<T>
    where
        T: serde::de::DeserializeOwned + Debug,
    {
        let raw_values =
            record
                .info(tag)
                .integer()?
                .ok_or_else(|| errors::Error::InvalidBCFRecord {
                    msg: "No varlociraptor observations found in record.".to_owned(),
                })?;

        // decode from i32 to u16 to u8
        let mut values_u8 = Vec::new();
        for v in raw_values {
            let mut buf = [0; 2];
            LittleEndian::write_u16(&mut buf, *v as u16);
            values_u8.extend(&buf);
        }

        // deserialize
        let values = bincode::deserialize(&values_u8)?;

        Ok(values)
    }

    let prob_mapping: Vec<MiniLogProb> = read_values(record, b"PROB_MAPPING")?;
    let prob_ref: Vec<MiniLogProb> = read_values(record, b"PROB_REF")?;
    let prob_alt: Vec<MiniLogProb> = read_values(record, b"PROB_ALT")?;
    let prob_missed_allele: Vec<MiniLogProb> = read_values(record, b"PROB_MISSED_ALLELE")?;
    let prob_sample_alt: Vec<MiniLogProb> = read_values(record, b"PROB_SAMPLE_ALT")?;
    let prob_double_overlap: Vec<MiniLogProb> = read_values(record, b"PROB_DOUBLE_OVERLAP")?;
    let prob_any_strand: Vec<MiniLogProb> = read_values(record, b"PROB_ANY_STRAND")?;
    let forward_strand: BitVec<u8> = read_values(record, b"FORWARD_STRAND")?;
    let reverse_strand: BitVec<u8> = read_values(record, b"REVERSE_STRAND")?;

    let obs = (0..prob_mapping.len())
        .map(|i| {
            ObservationBuilder::default()
                .prob_mapping_mismapping(prob_mapping[i].to_logprob())
                .prob_alt(prob_alt[i].to_logprob())
                .prob_ref(prob_ref[i].to_logprob())
                .prob_missed_allele(prob_missed_allele[i].to_logprob())
                .prob_sample_alt(prob_sample_alt[i].to_logprob())
                .prob_overlap(prob_double_overlap[i].to_logprob())
                .prob_any_strand(prob_any_strand[i].to_logprob())
                .forward_strand(forward_strand[i as u64])
                .reverse_strand(reverse_strand[i as u64])
                .build()
                .unwrap()
        })
        .collect_vec();

    Ok(obs)
}

pub(crate) fn write_observations(
    observations: &[Observation],
    record: &mut bcf::Record,
) -> Result<()> {
    let vec = || Vec::with_capacity(observations.len());
    let mut prob_mapping = vec();
    let mut prob_ref = vec();
    let mut prob_alt = vec();
    let mut prob_missed_allele = vec();
    let mut prob_sample_alt = vec();
    let mut prob_double_overlap = vec();
    let mut prob_any_strand = vec();
    let mut forward_strand: BitVec<u8> = BitVec::with_capacity(observations.len() as u64);
    let mut reverse_strand: BitVec<u8> = BitVec::with_capacity(observations.len() as u64);
    let encode_logprob = |prob: LogProb| utils::MiniLogProb::new(prob);

    for obs in observations {
        prob_mapping.push(encode_logprob(obs.prob_mapping));
        prob_ref.push(encode_logprob(obs.prob_ref));
        prob_alt.push(encode_logprob(obs.prob_alt));
        prob_missed_allele.push(encode_logprob(obs.prob_missed_allele));
        prob_sample_alt.push(encode_logprob(obs.prob_sample_alt));
        prob_double_overlap.push(encode_logprob(obs.prob_double_overlap));
        prob_any_strand.push(encode_logprob(obs.prob_any_strand));
        forward_strand.push(obs.forward_strand);
        reverse_strand.push(obs.reverse_strand);
    }

    fn push_values<T>(record: &mut bcf::Record, tag: &[u8], values: &T) -> Result<()>
    where
        T: serde::Serialize + Debug,
    {
        // serialize
        let mut encoded_values = bincode::serialize(values)?;

        // add padding zero if length is odd
        if encoded_values.len() % 2 > 0 {
            encoded_values.push(0);
        }

        // Encode as i32 (must first encode as u16, because the maximum i32 is used internally by BCF to indicate vector end)
        // This should not cause much wasted space, because similar (empty) bytes will be compressed away.
        let values_i32 = (0..encoded_values.len())
            .step_by(2)
            .map(|i| LittleEndian::read_u16(&encoded_values[i..i + 2]) as i32)
            .collect_vec();
        // write to record
        record.push_info_integer(tag, &values_i32)?;

        Ok(())
    }

    push_values(record, b"PROB_MAPPING", &prob_mapping)?;
    push_values(record, b"PROB_REF", &prob_ref)?;
    push_values(record, b"PROB_ALT", &prob_alt)?;
    push_values(record, b"PROB_MISSED_ALLELE", &prob_missed_allele)?;
    push_values(record, b"PROB_SAMPLE_ALT", &prob_sample_alt)?;
    push_values(record, b"PROB_DOUBLE_OVERLAP", &prob_double_overlap)?;
    push_values(record, b"PROB_ANY_STRAND", &prob_any_strand)?;
    push_values(record, b"FORWARD_STRAND", &forward_strand)?;
    push_values(record, b"REVERSE_STRAND", &reverse_strand)?;

    Ok(())
}

pub(crate) fn remove_observation_header_entries(header: &mut bcf::Header) {
    header.remove_info(b"PROB_MAPPING");
    header.remove_info(b"PROB_REF");
    header.remove_info(b"PROB_ALT");
    header.remove_info(b"PROB_MISSED_ALLELE");
    header.remove_info(b"PROB_SAMPLE_ALT");
    header.remove_info(b"PROB_DOUBLE_OVERLAP");
    header.remove_info(b"PROB_ANY_STRAND");
    header.remove_info(b"FORWARD_STRAND");
    header.remove_info(b"REVERSE_STRAND");
}

pub(crate) fn read_preprocess_options<P: AsRef<Path>>(bcfpath: P) -> Result<cli::Varlociraptor> {
    let reader = bcf::Reader::from_path(&bcfpath)?;
    for rec in reader.header().header_records() {
        if let bcf::header::HeaderRecord::Generic { ref key, ref value } = rec {
            if key == "varlociraptor_preprocess_args" {
                return Ok(serde_json::from_str(value)?);
            }
        }
    }
    Err(errors::Error::InvalidObservations {
        path: bcfpath.as_ref().to_owned(),
    }
    .into())
}

struct WorkItem {
    start: u64,
    chrom: String,
    variants: Vec<model::Variant>,
    record_id: Vec<u8>,
    record_mateid: Option<Vec<u8>>,
    record_index: usize,
}

#[derive(Derefable, new, Debug)]
struct Calls {
    index: usize,
    #[deref]
    inner: Vec<Call>,
}

impl utils::worker_pool::Orderable for Calls {
    fn index(&self) -> usize {
        self.index
    }
}
