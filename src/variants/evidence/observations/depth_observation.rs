// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use anyhow::Result;
use bio::stats::LogProb;

// use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use statrs::distribution::{Discrete, Poisson};

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::utils::homopolymers::HomopolymerErrorModel;
use crate::variants::evidence::observations::observation::Evidence;
use crate::variants::sample;
use crate::variants::types::DepthVariant;

use crate::variants::evidence::realignment::Realignable;

use super::fragment_id_factory::FragmentIdFactory;

#[derive(Debug, Clone, Getters)]
pub struct DepthObservation {
    #[getset(get = "pub")]
    pub cnv_probs: Vec<LogProb>,
}

impl DepthObservation {
    pub fn new(cnv_probs: Vec<LogProb>) -> Self {
        Self { cnv_probs }
    }

    pub fn compute_cnv_probs(
        cnv_positions_depth: Vec<f64>,
        avg_depth: u32,
        max_number_cn: usize,
    ) -> Self {
        // TODO How to know the ploidy?
        let ploidy = 2.0;

        let cnv_probs = (0..=max_number_cn)
            .map(|cn| {
                let lambda = (cn as f64 / ploidy) * avg_depth as f64;

                let log_likelihood = if lambda == 0.0 {
                    if cnv_positions_depth.iter().all(|&d| d == 0.0) {
                        LogProb(0.0)
                    } else {
                        LogProb(f64::NEG_INFINITY)
                    }
                } else {
                    let poisson = Poisson::new(lambda).unwrap();
                    let cnv_prob: LogProb = cnv_positions_depth
                        .iter()
                        .map(|&d| LogProb::from(poisson.pmf(d as u64)))
                        .sum();
                    cnv_prob
                };

                log_likelihood
            })
            .collect();
        Self { cnv_probs }
    }
}

/// Something that can be converted into observations.
pub(crate) trait DepthObservable: DepthVariant {
    fn extract_observations(
        &self,
        buffer: &mut sample::RecordBuffer,
        alignment_properties: &mut AlignmentProperties,
        max_depth: usize,
        alt_variants: &[Box<dyn Realignable>],
        observation_id_factory: &mut Option<&mut FragmentIdFactory>,
        max_number_cn: usize,
    ) -> Result<Vec<DepthObservation>>;

    /// Convert MAPQ (from read mapper) to LogProb for the event that the read maps
    /// correctly.
    fn prob_mapping(&self, evidence: &Evidence) -> LogProb;

    /// Return the minimum MAPQ of all records involved in the given evidence.
    fn min_mapq(&self, evidence: &Evidence) -> u8;

    /// Calculate an observation from the given evidence.
    fn evidence_to_observation(
        &self,
        evidence: &[Evidence],
        alignment_properties: &mut AlignmentProperties,
        homopolymer_error_model: &Option<HomopolymerErrorModel>,
        alt_variants: &[Box<dyn Realignable>],
        observation_id_factory: &mut Option<&mut FragmentIdFactory>,
        max_number_cn: usize,
    ) -> Result<Option<DepthObservation>> {
        // let id = observation_id_factory
        //     .as_mut()
        //     .map(|factory| factory.register(evidence));

        Ok(
            match self.allele_support(evidence, alignment_properties, alt_variants)? {
                // METHOD: for precise variants,
                // only consider allele support if it comes either from forward or reverse strand.
                // Unstranded observations (e.g. only insert size), are too unreliable, or do not contain
                // any information (e.g. no overlap).
                Some(cnv_read_depths) => {
                    let mut obs = DepthObservation::compute_cnv_probs(
                        cnv_read_depths,
                        alignment_properties.avg_depth,
                        max_number_cn,
                    );
                    // if allele_support.strand() != Strand::None || self.is_imprecise() =>
                    // {
                    // let alt_indel_len = allele_support.homopolymer_indel_len().unwrap_or(0);

                    // let mut obs = DepthObservationBuilder::default();
                    // obs.name(Some(evidence.id().to_string().to_owned()))
                    //     .fragment_id(id)
                    //     .prob_mapping_mismapping(self.prob_mapping(evidence))
                    //     .prob_alt(allele_support.prob_alt_allele())
                    //     .prob_ref(allele_support.prob_ref_allele())
                    //     .prob_sample_alt(self.prob_sample_alt(evidence, alignment_properties))
                    //     .prob_missed_allele(allele_support.prob_missed_allele())
                    //     .prob_overlap(if allele_support.strand() == Strand::Both {
                    //         LogProb::ln_one()
                    //     } else {
                    //         LogProb::ln_zero()
                    //     })
                    //     .strand(allele_support.strand())
                    //     .read_orientation(evidence.read_orientation()?)
                    //     .softclipped(evidence.softclipped())
                    //     .read_position(allele_support.read_position())
                    //     .paired(evidence.is_paired())
                    //     .prob_hit_base(LogProb::ln_one() - LogProb((evidence.len() as f64).ln()))
                    //     .is_max_mapq(self.min_mapq(evidence) == alignment_properties.max_mapq)
                    //     .alt_locus(evidence.alt_loci())
                    //     .third_allele_evidence(allele_support.third_allele_evidence().map(|d| *d));

                    // if let Some(homopolymer_error_model) = homopolymer_error_model {
                    //     let ref_indel_len =
                    //         alt_indel_len + homopolymer_error_model.variant_homopolymer_indel_len();

                    //     obs.homopolymer_indel_len(Some(ref_indel_len));

                    //     if ref_indel_len == 0 || alt_indel_len == 0 {
                    //         // no homopolymer indel in read compared to reference
                    //         obs.prob_observable_at_homopolymer_artifact(None)
                    //             .prob_observable_at_homopolymer_variant(None);
                    //     } else {
                    //         obs.prob_observable_at_homopolymer_variant(Some(
                    //             homopolymer_error_model.prob_observable(alt_indel_len),
                    //         ))
                    //         .prob_observable_at_homopolymer_artifact(Some(
                    //             homopolymer_error_model.prob_observable(ref_indel_len),
                    //         ));
                    //     }
                    // } else {
                    //     obs.homopolymer_indel_len(None)
                    //         .prob_observable_at_homopolymer_artifact(None)
                    //         .prob_observable_at_homopolymer_variant(None);
                    // }

                    // Some(obs.build().unwrap())
                    Some(obs)
                }
                _ => None,
            },
        )
    }
}
