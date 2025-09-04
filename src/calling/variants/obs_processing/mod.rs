// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::fmt::Debug;
use std::path::{Path, PathBuf};
use std::str;

use crate::cli;
use crate::errors;
use crate::utils;
use crate::utils::MiniLogProb;
use crate::variants::evidence::observations::pileup::Pileup;
use crate::variants::evidence::observations::read_observation::{
    AltLocus, ReadObservationBuilder, ReadPosition, Strand,
};
use anyhow::{bail, Context, Result};
use bio_types::genome::{self, AbstractLocus};
use bio_types::sequence::SequenceReadPairOrientation;
use bv::BitVec;
use byteorder::{ByteOrder, LittleEndian};
use itertools::Itertools;
use rust_htslib::bam::{self, Read as BAMRead};
use rust_htslib::bcf::{self, Read as BCFRead};

pub(crate) mod haplotype_feature_index;
pub(crate) mod observation_processor;
pub(crate) mod preprocessing;

pub(crate) static OBSERVATION_FORMAT_VERSION: &str = "15";

pub struct Observations {
    pub pileup: Pileup,
    pub is_homopolymer_indel: bool,
}

/// Read observations from BCF record.
pub fn read_observations(record: &mut bcf::Record) -> Result<Observations> {
    fn read_values<T>(record: &mut bcf::Record, tag: &[u8], allow_missing: bool) -> Result<T>
    where
        T: serde::de::DeserializeOwned + Debug + Default,
    {
        let info = record.info(tag).integer()?;
        if allow_missing && info.is_none() {
            return Ok(T::default());
        }
        let raw_values = info.ok_or_else(|| {
            errors::invalid_bcf_record(
                record.contig(),
                record.pos(),
                "No varlociraptor observations found in record.",
            )
        })?;

        // decode from i32 to u16 to u8
        let mut values_u8 = Vec::new();
        for v in raw_values.iter() {
            let mut buf = [0; 2];
            LittleEndian::write_u16(&mut buf, *v as u16);
            values_u8.extend(&buf);
        }

        // deserialize
        let values = bincode::deserialize(&values_u8)?;

        Ok(values)
    }

    let ids: Vec<Option<u64>> = read_values(record, b"FRAGMENT_ID", false)?;
    let prob_mapping: Vec<MiniLogProb> = read_values(record, b"PROB_MAPPING", false)?;
    let prob_ref: Vec<MiniLogProb> = read_values(record, b"PROB_REF", false)?;
    let prob_alt: Vec<MiniLogProb> = read_values(record, b"PROB_ALT", false)?;
    let prob_missed_allele: Vec<MiniLogProb> = read_values(record, b"PROB_MISSED_ALLELE", false)?;
    let prob_sample_alt: Vec<MiniLogProb> = read_values(record, b"PROB_SAMPLE_ALT", false)?;
    let prob_double_overlap: Vec<MiniLogProb> = read_values(record, b"PROB_DOUBLE_OVERLAP", false)?;
    let prob_hit_base: Vec<MiniLogProb> = read_values(record, b"PROB_HIT_BASE", false)?;
    let strand: Vec<Strand> = read_values(record, b"STRAND", false)?;
    let read_orientation: Vec<SequenceReadPairOrientation> =
        read_values(record, b"READ_ORIENTATION", false)?;
    let read_position: Vec<ReadPosition> = read_values(record, b"READ_POSITION", false)?;
    let softclipped: BitVec<u8> = read_values(record, b"SOFTCLIPPED", false)?;
    let paired: BitVec<u8> = read_values(record, b"PAIRED", false)?;
    let prob_observable_at_homopolymer_artifact: Vec<Option<MiniLogProb>> =
        read_values(record, b"PROB_HOMOPOLYMER_ARTIFACT_OBSERVABLE", true)?;
    let prob_observable_at_homopolymer_variant: Vec<Option<MiniLogProb>> =
        read_values(record, b"PROB_HOMOPOLYMER_VARIANT_OBSERVABLE", true)?;
    let homopolymer_indel_len: Vec<Option<i8>> =
        read_values(record, b"HOMOPOLYMER_INDEL_LEN", true)?;
    let is_homopolymer_indel = !prob_observable_at_homopolymer_artifact.is_empty();
    let is_max_mapq: BitVec<u8> = read_values(record, b"IS_MAX_MAPQ", false)?;
    let alt_locus: Vec<AltLocus> = read_values(record, b"ALT_LOCUS", false)?;
    let third_allele_evidence: Vec<Option<u32>> =
        read_values(record, b"THIRD_ALLELE_EVIDENCE", false)?;

    let read_obs = (0..prob_mapping.len())
        .map(|i| {
            let mut obs = ReadObservationBuilder::default();
            obs.name(None) // we do not pass the read names to the calling process
                .fragment_id(ids[i])
                .prob_mapping_mismapping(prob_mapping[i].to_logprob())
                .prob_alt(prob_alt[i].to_logprob())
                .prob_ref(prob_ref[i].to_logprob())
                .prob_missed_allele(prob_missed_allele[i].to_logprob())
                .prob_sample_alt(prob_sample_alt[i].to_logprob())
                .prob_overlap(prob_double_overlap[i].to_logprob())
                .prob_hit_base(prob_hit_base[i].to_logprob())
                .strand(strand[i])
                .read_orientation(read_orientation[i])
                .read_position(read_position[i])
                .softclipped(softclipped[i as u64])
                .paired(paired[i as u64])
                .is_max_mapq(is_max_mapq[i as u64])
                .alt_locus(alt_locus[i])
                .third_allele_evidence(third_allele_evidence[i]);

            if is_homopolymer_indel {
                obs.homopolymer_indel_len(homopolymer_indel_len[i])
                    .prob_observable_at_homopolymer_artifact(
                        prob_observable_at_homopolymer_artifact[i].map(|prob| prob.to_logprob()),
                    )
                    .prob_observable_at_homopolymer_variant(
                        prob_observable_at_homopolymer_variant[i].map(|prob| prob.to_logprob()),
                    );
            } else {
                obs.homopolymer_indel_len(None)
                    .prob_observable_at_homopolymer_artifact(None)
                    .prob_observable_at_homopolymer_variant(None);
            }
            obs.build().unwrap()
        })
        .collect_vec();

    let depth_obs = Vec::new(); // TODO: read depth observations!

    Ok(Observations {
        pileup: Pileup::new(read_obs, depth_obs),
        is_homopolymer_indel,
    })
}

pub(crate) fn write_observations(pileup: &Pileup, record: &mut bcf::Record) -> Result<()> {
    // TODO: write depth observations
    let read_observations = pileup.read_observations();

    let vec = || Vec::with_capacity(read_observations.len());
    let mut ids = Vec::with_capacity(read_observations.len());
    let mut prob_mapping = vec();
    let mut prob_ref = vec();
    let mut prob_alt = vec();
    let mut prob_missed_allele = vec();
    let mut prob_sample_alt = vec();
    let mut prob_double_overlap = vec();
    let mut strand = Vec::with_capacity(read_observations.len());
    let mut read_orientation = Vec::with_capacity(read_observations.len());
    let mut softclipped: BitVec<u8> = BitVec::with_capacity(read_observations.len() as u64);
    let mut paired: BitVec<u8> = BitVec::with_capacity(read_observations.len() as u64);
    let mut read_position = Vec::with_capacity(read_observations.len());
    let mut prob_hit_base = vec();
    let mut prob_observable_at_homopolymer_artifact: Vec<Option<MiniLogProb>> =
        Vec::with_capacity(read_observations.len());
    let mut prob_observable_at_homopolymer_variant: Vec<Option<MiniLogProb>> =
        Vec::with_capacity(read_observations.len());
    let mut homopolymer_indel_len: Vec<Option<i8>> = Vec::with_capacity(read_observations.len());
    let mut is_max_mapq: BitVec<u8> = BitVec::with_capacity(read_observations.len() as u64);
    let mut alt_locus = Vec::with_capacity(read_observations.len());
    let mut third_allele_evidence = Vec::with_capacity(read_observations.len());

    let encode_logprob = utils::MiniLogProb::new;
    for obs in read_observations {
        ids.push(obs.fragment_id);
        prob_mapping.push(encode_logprob(obs.prob_mapping()));
        prob_ref.push(encode_logprob(obs.prob_ref()));
        prob_alt.push(encode_logprob(obs.prob_alt()));
        prob_missed_allele.push(encode_logprob(obs.prob_missed_allele));
        prob_sample_alt.push(encode_logprob(obs.prob_sample_alt));
        prob_double_overlap.push(encode_logprob(obs.prob_double_overlap));
        prob_hit_base.push(encode_logprob(obs.prob_hit_base));
        strand.push(obs.strand);
        read_orientation.push(obs.read_orientation);
        softclipped.push(obs.softclipped);
        paired.push(obs.paired);
        read_position.push(obs.read_position);
        is_max_mapq.push(obs.is_max_mapq);
        alt_locus.push(obs.alt_locus);
        third_allele_evidence.push(obs.third_allele_evidence);

        prob_observable_at_homopolymer_artifact.push(
            obs.prob_observable_at_homopolymer_artifact
                .map(encode_logprob),
        );
        prob_observable_at_homopolymer_variant.push(
            obs.prob_observable_at_homopolymer_variant
                .map(encode_logprob),
        );
        homopolymer_indel_len.push(obs.homopolymer_indel_len);
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

    push_values(record, b"FRAGMENT_ID", &ids)?;
    push_values(record, b"PROB_MAPPING", &prob_mapping)?;
    push_values(record, b"PROB_REF", &prob_ref)?;
    push_values(record, b"PROB_ALT", &prob_alt)?;
    push_values(record, b"PROB_MISSED_ALLELE", &prob_missed_allele)?;
    push_values(record, b"PROB_SAMPLE_ALT", &prob_sample_alt)?;
    push_values(record, b"PROB_DOUBLE_OVERLAP", &prob_double_overlap)?;
    push_values(record, b"STRAND", &strand)?;
    push_values(record, b"READ_ORIENTATION", &read_orientation)?;
    push_values(record, b"SOFTCLIPPED", &softclipped)?;
    push_values(record, b"PAIRED", &paired)?;
    push_values(record, b"READ_POSITION", &read_position)?;
    push_values(record, b"PROB_HIT_BASE", &prob_hit_base)?;
    push_values(record, b"IS_MAX_MAPQ", &is_max_mapq)?;
    push_values(record, b"ALT_LOCUS", &alt_locus)?;
    push_values(record, b"THIRD_ALLELE_EVIDENCE", &third_allele_evidence)?;

    if prob_observable_at_homopolymer_artifact
        .iter()
        .any(|prob| prob.is_some())
    {
        // only record values if there is any homopolymer error observation
        push_values(
            record,
            b"PROB_HOMOPOLYMER_ARTIFACT_OBSERVABLE",
            &prob_observable_at_homopolymer_artifact,
        )?;
        push_values(
            record,
            b"PROB_HOMOPOLYMER_VARIANT_OBSERVABLE",
            &prob_observable_at_homopolymer_variant,
        )?;
        push_values(record, b"HOMOPOLYMER_INDEL_LEN", &homopolymer_indel_len)?;
    }

    Ok(())
}

pub(crate) fn remove_observation_header_entries(header: &mut bcf::Header) {
    header.remove_info(b"FRAGMENT_ID");
    header.remove_info(b"PROB_MAPPING");
    header.remove_info(b"PROB_REF");
    header.remove_info(b"PROB_ALT");
    header.remove_info(b"PROB_MISSED_ALLELE");
    header.remove_info(b"PROB_SAMPLE_ALT");
    header.remove_info(b"PROB_DOUBLE_OVERLAP");
    header.remove_info(b"STRAND");
    header.remove_info(b"READ_ORIENTATION");
    header.remove_info(b"SOFTCLIPPED");
    header.remove_info(b"PAIRED");
    header.remove_info(b"PROB_HIT_BASE");
    header.remove_info(b"READ_POSITION");
    header.remove_info(b"PROB_HOMOPOLYMER_ARTIFACT_OBSERVABLE");
    header.remove_info(b"PROB_HOMOPOLYMER_VARIANT_OBSERVABLE");
    header.remove_info(b"HOMOPOLYMER_INDEL_LEN");
    header.remove_info(b"IS_MAX_MAPQ");
    header.remove_info(b"ALT_LOCUS");
    header.remove_info(b"THIRD_ALLELE_EVIDENCE");
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
