// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::error::Error;
use std::fs;
use std::path::Path;

use bio::io::fasta;
use bio::stats::LogProb;
use derive_builder::Builder;
use itertools::Itertools;
use rust_htslib::bcf::{self, Read};

use crate::errors;
use crate::model::evidence::{Observation, observation::ObservationBuilder};
use crate::model::sample::Sample;
use crate::utils;
use crate::calling::variants::{chrom, Call, CallBuilder, VariantBuilder};

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct ObservationProcessor {
    sample: Sample,
    #[builder(private)]
    reference_buffer: utils::ReferenceBuffer,
    #[builder(private)]
    bcf_reader: bcf::Reader,
    #[builder(private)]
    bcf_writer: bcf::Writer,
    omit_snvs: bool,
    omit_indels: bool,
    max_indel_len: u32,
}

impl ObservationProcessorBuilder {
    pub fn reference(self, reader: fasta::IndexedReader<fs::File>) -> Result<Self, Box<dyn Error>> {
        Ok(self.reference_buffer(utils::ReferenceBuffer::new(reader)))
    }

    pub fn inbcf<P: AsRef<Path>>(self, path: Option<P>) -> Result<Self, Box<dyn Error>> {
        Ok(self.bcf_reader(if let Some(path) = path {
            bcf::Reader::from_path(path)?
        } else {
            bcf::Reader::from_stdin()?
        }))
    }

    pub fn outbcf<P: AsRef<Path>>(self, path: Option<P>) -> Result<Self, Box<dyn Error>> {
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

        let writer = if let Some(path) = path {
            bcf::Writer::from_path(path, &header, false, bcf::Format::BCF)?
        } else {
            bcf::Writer::from_stdout(&header, false, bcf::Format::BCF)?
        };
        Ok(self.bcf_writer(writer))
    }
}

impl ObservationProcessor {
    pub fn process(&mut self) -> Result<(), Box<dyn Error>> {
        let mut i = 0;
        let mut record = self.bcf_reader.empty_record();
        loop {
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

    fn process_record(
        &mut self,
        record: &mut bcf::Record,
    ) -> Result<Option<Call>, Box<dyn Error>> {
        let start = record.pos();
        let chrom = chrom(&self.bcf_reader, &record);
        let variants = utils::collect_variants(
            record,
            self.omit_snvs,
            self.omit_indels,
            Some(0..self.max_indel_len + 1),
        )?;

        if variants.is_empty() || variants.iter().all(|v| v.is_none()) {
            return Ok(None);
        }

        let mut call = CallBuilder::default()
            .chrom(chrom.to_owned())
            .pos(start)
            .id({
                let id = record.id();
                if id == b"." {
                    None
                } else {
                    Some(id)
                }
            })
            .variants(Vec::new())
            .build()?;

        for variant in variants.into_iter() {
            if let Some(variant) = variant {
                let chrom_seq = self.reference_buffer.seq(&chrom)?;
                let pileup = self.sample.extract_observations(start, &variant, chrom, chrom_seq)?;

                let start = start as usize;

                // add variant information
                call.variants.push(
                    VariantBuilder::default()
                        .variant(&variant, start, chrom_seq)
                        .observations(Some(pileup))
                        .build()?
                );
            }
        }

        Ok(Some(call))
    }
}



/// Read observations from BCF record.
pub fn read_observations<'a>(record: &'a mut bcf::Record) -> Result<Vec<Observation>, Box<dyn Error>> {
    fn read_float_entry<'a>(record: &'a mut bcf::Record, entry_name: &'a [u8]) -> Result<&'a [f32], Box<dyn Error>> {
        Ok(record.info(entry_name).float()?.ok_or_else(|| errors::Error::InvalidBCFRecord { msg: "No varlociraptor observations found in record.".to_owned() } )?)
    };

    fn read_string_entry<'a>(record: &'a mut bcf::Record, entry_name: &'a [u8]) ->  Result<Vec<&'a [u8]>, Box<dyn Error>> {
        Ok(record.info(entry_name).string()?.ok_or_else(|| errors::Error::InvalidBCFRecord { msg: "No varlociraptor observations found in record.".to_owned() } )?)
    };

    let prob_mapping = read_float_entry(record, b"PROB_MAPPING")?.to_owned();
    let prob_ref = read_float_entry(record, b"PROB_REF")?.to_owned();
    let prob_alt = read_float_entry(record, b"PROB_ALT")?.to_owned();
    let prob_missed_allele = read_float_entry(record, b"PROB_MISSED_ALLELE")?.to_owned();
    let prob_sample_alt = read_float_entry(record, b"PROB_SAMPLE_ALT")?.to_owned();
    let prob_double_overlap = read_float_entry(record, b"PROB_DOUBLE_OVERLAP")?.to_owned();
    let prob_any_strand = read_float_entry(record, b"PROB_ANY_STRAND")?.to_owned();
    let forward_strand = read_string_entry(record, b"FORWARD_STRAND")?[0].to_owned();
    let reverse_strand = read_string_entry(record, b"REVERSE_STRAND")?[0].to_owned();

    let decode_bool = |value| if value == b'1' { true } else { false };
    
    Ok((0..prob_mapping.len()).map(|i| {
        ObservationBuilder::default()
            .prob_mapping_mismapping(LogProb(prob_mapping[i] as f64))
            .prob_alt(LogProb(prob_alt[i] as f64))
            .prob_ref(LogProb(prob_ref[i] as f64))
            .prob_missed_allele(LogProb(prob_missed_allele[i] as f64))
            .prob_sample_alt(LogProb(prob_sample_alt[i] as f64))
            .prob_overlap(LogProb(prob_double_overlap[i] as f64))
            .prob_any_strand(LogProb(prob_any_strand[i] as f64))
            .forward_strand(decode_bool(forward_strand[i]))
            .reverse_strand(decode_bool(reverse_strand[i]))
            .build()
            .unwrap()
    }).collect_vec())
}

pub fn write_observations(observations: &[Observation], record: &mut bcf::Record) -> Result<(), Box<dyn Error>> {
    let vec = || Vec::with_capacity(observations.len());
    let mut prob_mapping = vec();
    let mut prob_ref = vec();
    let mut prob_alt = vec();
    let mut prob_missed_allele = vec();
    let mut prob_sample_alt = vec();
    let mut prob_double_overlap = vec();
    let mut prob_any_strand = vec();
    let mut forward_strand = Vec::with_capacity(observations.len());
    let mut reverse_strand = Vec::with_capacity(observations.len());
    let encode_bool = |value| if value { b'1' } else { b'0' };

    // TODO use bincode to compress as f16 and store in i32 blocks.
    for obs in observations {
        prob_mapping.push(*obs.prob_mapping as f32);
        prob_ref.push(*obs.prob_ref as f32);
        prob_alt.push(*obs.prob_alt as f32);
        prob_missed_allele.push(*obs.prob_missed_allele as f32);
        prob_sample_alt.push(*obs.prob_sample_alt as f32);
        prob_double_overlap.push(*obs.prob_double_overlap as f32);
        prob_any_strand.push(*obs.prob_any_strand as f32);
        forward_strand.push(encode_bool(obs.forward_strand));
        reverse_strand.push(encode_bool(obs.reverse_strand));
    }

    record.push_info_float(b"PROB_MAPPING", &prob_mapping)?;
    record.push_info_float(b"PROB_REF", &prob_ref)?;
    record.push_info_float(b"PROB_ALT", &prob_alt)?;
    record.push_info_float(b"PROB_MISSED_ALLELE", &prob_missed_allele)?;
    record.push_info_float(b"PROB_SAMPLE_ALT", &prob_sample_alt)?;
    record.push_info_float(b"PROB_DOUBLE_OVERLAP", &prob_double_overlap)?;
    record.push_info_float(b"PROB_ANY_STRAND", &prob_any_strand)?;
    record.push_info_string(b"FORWARD_STRAND", &[&forward_strand])?;
    record.push_info_string(b"REVERSE_STRAND", &[&reverse_strand])?;

    Ok(())
}

pub fn remove_observation_header_entries(header: &mut bcf::Header) {
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