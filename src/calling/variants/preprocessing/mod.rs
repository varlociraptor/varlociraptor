// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::fmt::Debug;
use std::fs;
use std::path::Path;
use std::str;
use std::sync::Arc;

use anyhow::Result;
use bincode;
use bio::io::fasta;
use bio::stats::LogProb;
use bio_types::genome;
use bv::BitVec;
use byteorder::{ByteOrder, LittleEndian};
use derive_builder::Builder;
use itertools::Itertools;
use rust_htslib::bcf::{self, Read};
use serde_json;

use crate::calling::variants::{chrom, Call, CallBuilder, VariantBuilder};
use crate::cli;
use crate::errors;
use crate::reference;
use crate::utils;
use crate::utils::MiniLogProb;
use crate::variants;
use crate::variants::evidence::observation::{Observation, ObservationBuilder};
use crate::variants::evidence::realignment;
use crate::variants::model;
use crate::variants::sample::Sample;

#[derive(Builder)]
#[builder(pattern = "owned")]
pub(crate) struct ObservationProcessor {
    sample: Sample,
    #[builder(private)]
    reference_buffer: reference::Buffer,
    #[builder(private)]
    bcf_reader: bcf::Reader,
    #[builder(private)]
    bcf_writer: bcf::Writer,
    omit_snvs: bool,
    omit_indels: bool,
    gap_params: realignment::pairhmm::GapParams,
    realignment_window: u64,
}

impl ObservationProcessorBuilder {
    pub(crate) fn reference(self, reader: fasta::IndexedReader<fs::File>) -> Result<Self> {
        Ok(self.reference_buffer(reference::Buffer::new(reader)))
    }

    pub(crate) fn inbcf<P: AsRef<Path>>(self, path: Option<P>) -> Result<Self> {
        Ok(self.bcf_reader(if let Some(path) = path {
            bcf::Reader::from_path(path)?
        } else {
            bcf::Reader::from_stdin()?
        }))
    }

    pub(crate) fn outbcf<P: AsRef<Path>>(
        self,
        path: Option<P>,
        options: &cli::Varlociraptor,
    ) -> Result<Self> {
        let mut header = bcf::Header::new();

        // register SVLEN
        header.push_record(
            b"##INFO=<ID=SVLEN,Number=A,Type=Integer,\
              Description=\"Difference in length between REF and ALT alleles\">",
        );

        // register sequences
        for sequence in self
            .reference_buffer
            .as_ref()
            .expect(".reference() has to be called before .outbcf()")
            .reader
            .index
            .sequences()
        {
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
                serde_json::to_string(options)?
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

        let writer = if let Some(path) = path {
            bcf::Writer::from_path(path, &header, false, bcf::Format::BCF)?
        } else {
            bcf::Writer::from_stdout(&header, false, bcf::Format::BCF)?
        };
        Ok(self.bcf_writer(writer))
    }
}

impl ObservationProcessor {
    pub(crate) fn process(&mut self) -> Result<()> {
        let mut i = 0;
        loop {
            let mut record = self.bcf_reader.empty_record();
            if !self.bcf_reader.read(&mut record)? {
                return Ok(());
            }

            i += 1;

            let call = self.process_record(&mut record)?;
            if let Some(call) = call {
                call.write_preprocessed_record(&mut self.bcf_writer)?;
            }
            if i % 100 == 0 {
                info!("{} records processed.", i);
            }
        }
    }

    fn process_record(&mut self, record: &mut bcf::Record) -> Result<Option<Call>> {
        let start = record.pos() as u64;
        let chrom = chrom(&self.bcf_reader, &record).to_owned();
        let variants = utils::collect_variants(record, self.omit_snvs, self.omit_indels, None)?;

        if variants.is_empty() {
            return Ok(None);
        }

        let mut call = CallBuilder::default()
            .chrom(chrom.to_owned())
            .pos(start as u64)
            .id({
                let id = record.id();
                if id == b"." {
                    None
                } else {
                    Some(id)
                }
            })
            .variants(Vec::new())
            .build()
            .unwrap();

        for variant in variants.into_iter() {
            let chrom_seq = self.reference_buffer.seq(&chrom)?;
            let pileup =
                self.extract_observations(start, &chrom, Arc::clone(&chrom_seq), &variant)?;

            let start = start as usize;

            // add variant information
            call.variants.push(
                VariantBuilder::default()
                    .variant(&variant, start, chrom_seq.as_ref())
                    .observations(Some(pileup))
                    .build()
                    .unwrap(),
            );
        }

        Ok(Some(call))
    }

    fn extract_observations(
        &mut self,
        start: u64,
        chrom: &[u8],
        chrom_seq: Arc<Vec<u8>>,
        variant: &model::Variant,
    ) -> Result<Vec<Observation>> {
        let chrom = str::from_utf8(chrom).unwrap().to_owned();
        let locus = || genome::Locus::new(chrom.clone(), start);
        let interval = |len: u64| genome::Interval::new(chrom.clone(), start..start + len);
        let start = start as usize;
        let realignment_window = self.realignment_window;
        let gap_params = self.gap_params.clone();

        let realigner =
            || realignment::Realigner::new(Arc::clone(&chrom_seq), gap_params, realignment_window);

        match variant {
            model::Variant::SNV(alt) => self
                .sample
                .extract_observations(&variants::types::SNV::new(locus(), chrom_seq[start], *alt)),
            model::Variant::MNV(alt) => {
                self.sample.extract_observations(&variants::types::MNV::new(
                    locus(),
                    chrom_seq[start..start + alt.len()].to_owned(),
                    alt.to_owned(),
                ))
            }
            model::Variant::None => self
                .sample
                .extract_observations(&variants::types::None::new(locus(), chrom_seq[start])),
            model::Variant::Deletion(l) => self
                .sample
                .extract_observations(&variants::types::Deletion::new(interval(*l), realigner())),
            model::Variant::Insertion(seq) => {
                self.sample
                    .extract_observations(&variants::types::Insertion::new(
                        locus(),
                        seq.to_owned(),
                        realigner(),
                    ))
            }
        }
    }
}

pub(crate) static OBSERVATION_FORMAT_VERSION: &'static str = "2";

/// Read observations from BCF record.
pub(crate) fn read_observations<'a>(record: &'a mut bcf::Record) -> Result<Vec<Observation>> {
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

pub(crate) fn write_observations(observations: &[Observation], record: &mut bcf::Record) -> Result<()> {
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
        match rec {
            bcf::header::HeaderRecord::Generic { ref key, ref value } => {
                if key == "varlociraptor_preprocess_args" {
                    return Ok(serde_json::from_str(value)?);
                }
            }
            _ => (),
        }
    }
    Err(errors::Error::InvalidObservations {
        path: bcfpath.as_ref().to_owned(),
    })?
}
