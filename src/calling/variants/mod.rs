// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

pub mod calling;
pub mod preprocessing;

use std::collections::HashMap;
use std::str;
use std::u8;

use anyhow::Result;
use bio::stats::{LogProb, PHREDProb};
use derive_builder::Builder;
use itertools::join;
use itertools::Itertools;
use rust_htslib::bcf::{self, record::Numeric, Read};
use vec_map::VecMap;

use crate::calling::variants::preprocessing::write_observations;
use crate::model;
use crate::model::evidence::observation::expected_depth;
use crate::model::evidence::Observation;
use crate::model::{AlleleFreq, StrandBias};
use crate::utils;

pub use crate::calling::variants::calling::{Caller, CallerBuilder};
pub use crate::calling::variants::preprocessing::{
    ObservationProcessor, ObservationProcessorBuilder,
};

#[derive(Default, Clone, Debug, Builder)]
pub struct Call {
    chrom: Vec<u8>,
    pos: u32,
    #[builder(default = "None")]
    id: Option<Vec<u8>>,
    #[builder(default = "Vec::new()")]
    variants: Vec<Variant>,
}

impl Call {
    pub fn write_preprocessed_record(&self, bcf_writer: &mut bcf::Writer) -> Result<()> {
        let rid = bcf_writer.header().name2rid(&self.chrom)?;
        for variant in self.variants.iter() {
            let mut record = bcf_writer.empty_record();
            record.set_rid(Some(rid));
            record.set_pos(self.pos as i64);
            // set ID if present
            if let Some(ref id) = self.id {
                record.set_id(id)?;
            }

            // set alleles
            record.set_alleles(&[&variant.ref_allele[..], &variant.alt_allele[..]])?;
            if let Some(svlen) = variant.svlen {
                record.push_info_integer(b"SVLEN", &[svlen])?;
            }

            // set qual
            record.set_qual(f32::missing());

            // add raw observations
            write_observations(variant.observations.as_ref().unwrap(), &mut record)?;

            bcf_writer.write(&record)?;
        }

        Ok(())
    }

    pub fn write_final_record(&self, bcf_writer: &mut bcf::Writer) -> Result<()> {
        let rid = bcf_writer.header().name2rid(&self.chrom)?;
        for (first_grouper, group) in self
            .variants
            .iter()
            .group_by(|variant| BCFGrouper(&variant))
            .into_iter()
        {
            let ref_allele = &first_grouper.0.ref_allele;
            let mut record = bcf_writer.empty_record();
            record.set_rid(Some(rid));
            record.set_pos(self.pos as i64);
            // set ID if present
            if let Some(ref id) = self.id {
                record.set_id(id)?;
            }

            let mut event_probs = HashMap::new();
            let mut allelefreq_estimates = VecMap::new();
            let mut observations = VecMap::new();
            let mut obs_counts = VecMap::new();
            let mut strand_bias = VecMap::new();
            let mut alleles = Vec::new();
            let mut svlens = Vec::new();
            alleles.push(&ref_allele[..]);

            // collect per group information
            for variant in group {
                alleles.push(&variant.alt_allele[..]);

                for (event, prob) in variant
                    .event_probs
                    .as_ref()
                    .expect("bug: event probs must be set")
                {
                    event_probs
                        .entry(event)
                        .or_insert_with(|| Vec::new())
                        .push(*prob);
                }

                for (i, sample_info) in variant
                    .sample_info
                    .as_ref()
                    .expect("bug: sample_info must be set")
                    .iter()
                    .enumerate()
                {
                    strand_bias.entry(i).or_insert_with(Vec::new).push(
                        match sample_info.strand_bias {
                            StrandBias::None => '.',
                            StrandBias::Forward => '+',
                            StrandBias::Reverse => '-',
                        },
                    );

                    allelefreq_estimates
                        .entry(i)
                        .or_insert_with(|| Vec::new())
                        .push(*sample_info.allelefreq_estimate as f32);

                    obs_counts
                        .entry(i)
                        .or_insert_with(|| Vec::new())
                        .push(expected_depth(&sample_info.observations) as i32);

                    observations.entry(i).or_insert_with(|| Vec::new()).push({
                        utils::generalized_cigar(
                            sample_info.observations.iter().map(|obs| {
                                let score = utils::evidence_kass_raftery_to_letter(
                                    obs.bayes_factor_alt().evidence_kass_raftery(),
                                );
                                format!(
                                    "{}{}",
                                    if obs.prob_mapping < LogProb(0.95_f64.ln()) {
                                        score.to_ascii_lowercase()
                                    } else {
                                        score.to_ascii_uppercase()
                                    },
                                    match (obs.forward_strand, obs.reverse_strand) {
                                        (true, true) => '*',
                                        (false, true) => '-',
                                        (true, false) => '+',
                                        _ => panic!("bug: unknown strandedness"),
                                    }
                                )
                            }),
                            false,
                        )
                    })
                }

                if let Some(svlen) = variant.svlen {
                    svlens.push(svlen);
                } else {
                    svlens.push(i32::missing());
                }
            }

            // set alleles
            record.set_alleles(&alleles)?;
            record.push_info_integer(b"SVLEN", &svlens)?;

            // set qual
            record.set_qual(f32::missing());

            // set event probabilities
            for (event, probs) in event_probs {
                let probs = probs
                    .iter()
                    .map(|p| PHREDProb::from(*p).abs() as f32)
                    .collect_vec();
                record.push_info_float(event_tag_name(event).as_bytes(), &probs)?;
            }

            // set sample info
            let dp = obs_counts.values().flatten().cloned().collect_vec();
            record.push_format_integer(b"DP", &dp)?;

            let afs = allelefreq_estimates
                .values()
                .flatten()
                .cloned()
                .collect_vec();
            record.push_format_float(b"AF", &afs)?;

            let obs = if observations.values().any(|obs| !obs.is_empty()) {
                observations
                    .values()
                    .map(|allele_obs| {
                        join(
                            allele_obs
                                .iter()
                                .map(|o| if o.is_empty() { "." } else { o }),
                            ",",
                        )
                        .into_bytes()
                    })
                    .collect_vec()
            } else {
                vec![b".".to_vec()]
            };
            record.push_format_string(b"OBS", &obs)?;

            let sb = strand_bias
                .values()
                .map(|sb| join(sb.iter(), ",").into_bytes())
                .collect_vec();

            record.push_format_string(b"SB", &sb)?;

            bcf_writer.write(&record)?;
        }

        Ok(())
    }
}

#[derive(Default, Clone, Debug, Builder)]
pub struct Variant {
    #[builder(private)]
    ref_allele: Vec<u8>,
    #[builder(private)]
    alt_allele: Vec<u8>,
    #[builder(private, default = "None")]
    svlen: Option<i32>,
    #[builder(private, default = "None")]
    event_probs: Option<HashMap<String, LogProb>>,
    #[builder(default = "None")]
    observations: Option<Vec<Observation>>,
    #[builder(default = "None")]
    sample_info: Option<Vec<SampleInfo>>,
}

impl VariantBuilder {
    /// Build the variant from a single-allele bcf record.
    pub fn record(&mut self, record: &mut bcf::Record) -> Result<&mut Self> {
        let alleles = record.alleles();
        Ok(self
            .ref_allele(alleles[0].to_owned())
            .alt_allele(alleles[1].to_owned())
            // TODO remove abs once https://github.com/samtools/bcftools/issues/874 is finally properly fixed
            .svlen(record.info(b"SVLEN").integer()?.map(|v| v[0].abs())))
    }

    pub fn variant(
        &mut self,
        variant: &model::Variant,
        start: usize,
        chrom_seq: &[u8],
    ) -> &mut Self {
        match variant {
            model::Variant::Deletion(l) => {
                let l = l.clone();
                // TODO make negative again once https://github.com/samtools/bcftools/issues/874 is finally properly fixed
                let svlen = l as i32;
                if l <= 50 {
                    self.ref_allele(chrom_seq[start..start + 1 + l as usize].to_ascii_uppercase())
                        .alt_allele(chrom_seq[start..start + 1].to_ascii_uppercase())
                        .svlen(Some(svlen))
                } else {
                    self.ref_allele(chrom_seq[start..start + 1].to_ascii_uppercase())
                        .alt_allele(b"<DEL>".to_ascii_uppercase())
                        .svlen(Some(svlen))
                }
            }
            model::Variant::Insertion(ref seq) => {
                let svlen = seq.len() as i32;
                let ref_allele = vec![chrom_seq[start]];
                let mut alt_allele = ref_allele.clone();
                alt_allele.extend(seq);

                self.ref_allele(ref_allele.to_ascii_uppercase())
                    .alt_allele(alt_allele.to_ascii_uppercase())
                    .svlen(Some(svlen))
            }
            model::Variant::SNV(base) => self
                .ref_allele(chrom_seq[start..start + 1].to_ascii_uppercase())
                .alt_allele(vec![*base].to_ascii_uppercase()),
            model::Variant::MNV(bases) => self
                .ref_allele(chrom_seq[start..start + bases.len()].to_ascii_uppercase())
                .alt_allele(bases.to_ascii_uppercase()),
            model::Variant::None => self
                .ref_allele(chrom_seq[start..start + 1].to_ascii_uppercase())
                .alt_allele(b"<REF>".to_ascii_uppercase()),
        }
    }
}

#[derive(Default, Clone, Debug, Builder)]
pub struct SampleInfo {
    allelefreq_estimate: AlleleFreq,
    #[builder(default = "Vec::new()")]
    observations: Vec<Observation>,
    strand_bias: StrandBias,
}

/// Wrapper for comparing alleles for compatibility in BCF files.
/// PartialEq::eq() returns true for all alleles that can occur in the same BCF record.
pub struct BCFGrouper<'a>(pub &'a Variant);

impl<'a> PartialEq for BCFGrouper<'a> {
    fn eq(&self, _other: &BCFGrouper) -> bool {
        // Currently, we want all variants to be in a separate record.
        // This might change again in the future.
        false

        // let s = self.0;
        // let o = other.0;
        // // Ensure that all compatible alleles have the same ref.
        // // Disallow two <DEL> alleles in the same record (because e.g. htsjdk fails then, many others likely as well).
        // s.ref_allele.eq(&o.ref_allele) && !(&s.alt_allele == b"<DEL>" && o.alt_allele == b"<DEL>")
    }
}

fn chrom<'a>(inbcf: &'a bcf::Reader, record: &bcf::Record) -> &'a [u8] {
    inbcf.header().rid2name(record.rid().unwrap()).unwrap()
}

pub fn event_tag_name(event: &str) -> String {
    format!("PROB_{}", event.to_ascii_uppercase())
}
