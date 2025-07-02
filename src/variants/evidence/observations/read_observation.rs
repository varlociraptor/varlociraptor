// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::char;
use std::hash::{Hash, Hasher};
use std::ops;
use std::rc::Rc;
use std::str;

use anyhow::Result;
use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::{LogProb, PHREDProb};
use bio_types::genome::{self, AbstractLocus};
use bio_types::sequence::SequenceReadPairOrientation;
use counter::Counter;

use rust_htslib::bam;

use serde::Serialize;
// use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::bayesian::BayesFactor;
use itertools::Itertools;

use super::pileup::Pileup;
use crate::errors::{self, Error};
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::utils::homopolymers::HomopolymerErrorModel;
use crate::utils::{self, PROB_05};
use crate::utils::{bayes_factor_to_letter, PROB_095};
use crate::variants::evidence::observations::observation::{
    calc_major_feature, locus_to_bucket, AltLocus, Evidence, ExactAltLoci, MaxBayesFactor,
    ReadPosition, Strand,
};
use crate::variants::sample;
use crate::variants::types::ReadVariant;

use crate::variants::evidence::realignment::Realignable;

use super::fragment_id_factory::FragmentIdFactory;

const INVALID_XA_FORMAT_MSG: &str = "XA tag of bam records in unexpected format. Expecting string (type Z) in bwa format (chr,pos,CIGAR,NM;).";

/// Calculate expected value of sequencing depth, considering mapping quality.
pub(crate) fn expected_depth(obs: &[ProcessedReadObservation]) -> u32 {
    LogProb::ln_sum_exp(&obs.iter().map(|o| o.prob_mapping).collect_vec())
        .exp()
        .round() as u32
}

pub type ProcessedReadObservation = ReadObservation<ReadPosition, AltLocus>;

impl<P: Clone, A: Clone> ReadObservationBuilder<P, A> {
    pub fn prob_mapping_mismapping(&mut self, prob_mapping: LogProb) -> &mut Self {
        self.prob_mapping(prob_mapping)
            .prob_mismapping(prob_mapping.ln_one_minus_exp())
    }

    pub fn prob_overlap(&mut self, prob_double_overlap: LogProb) -> &mut Self {
        self.prob_double_overlap(prob_double_overlap)
            .prob_single_overlap(prob_double_overlap.ln_one_minus_exp())
    }
}

impl ReadObservation<Option<u32>, ExactAltLoci> {
    pub(crate) fn process(
        &self,
        major_read_position: Option<u32>,
        major_alt_locus: &Option<genome::Locus>,
        alignment_properties: &AlignmentProperties,
    ) -> ReadObservation<ReadPosition, AltLocus> {
        ReadObservation {
            name: self.name.clone(),
            fragment_id: self.fragment_id,
            prob_mapping: self.prob_mapping,
            prob_mismapping: self.prob_mismapping,
            prob_mapping_adj: self.prob_mapping_adj,
            prob_mismapping_adj: self.prob_mismapping_adj,
            prob_alt: self.prob_alt,
            prob_ref: self.prob_ref,
            prob_alt_adj: self.prob_alt_adj,
            prob_ref_adj: self.prob_ref_adj,
            prob_missed_allele: self.prob_missed_allele,
            prob_sample_alt: self.prob_sample_alt,
            prob_double_overlap: self.prob_double_overlap,
            prob_single_overlap: self.prob_single_overlap,
            prob_hit_base: self.prob_hit_base,
            strand: self.strand,
            read_orientation: self.read_orientation,
            softclipped: self.softclipped,
            paired: self.paired,
            read_position: self.read_position.map_or(ReadPosition::Some, |pos| {
                if let Some(major_pos) = major_read_position {
                    if pos == major_pos {
                        ReadPosition::Major
                    } else {
                        ReadPosition::Some
                    }
                } else {
                    ReadPosition::Some
                }
            }),
            prob_observable_at_homopolymer_artifact: self.prob_observable_at_homopolymer_artifact,
            prob_observable_at_homopolymer_variant: self.prob_observable_at_homopolymer_variant,
            homopolymer_indel_len: self.homopolymer_indel_len,
            is_max_mapq: self.is_max_mapq,
            alt_locus: if let Some(major_alt_locus) = major_alt_locus {
                if self.alt_locus.iter().any(|alt_locus| {
                    locus_to_bucket(alt_locus, alignment_properties) == *major_alt_locus
                }) {
                    AltLocus::Major
                } else if self.alt_locus.is_empty() {
                    AltLocus::None
                } else {
                    AltLocus::Some
                }
            } else {
                AltLocus::None
            },
            third_allele_evidence: self.third_allele_evidence,
        }
    }
}

impl<P: Clone, A: Clone> ReadObservation<P, A> {
    // Represent evidence as bayes factor, returning either the one for ref against alt
    // or alt against ref, depending on which is higher
    pub(crate) fn max_bayes_factor(&self) -> MaxBayesFactor {
        let bf_alt = self.bayes_factor_alt();
        let bf_ref = self.bayes_factor_ref();
        if bf_alt > bf_ref {
            MaxBayesFactor::Alt(bf_alt)
        } else if bf_ref > bf_alt {
            MaxBayesFactor::Ref(bf_ref)
        } else {
            MaxBayesFactor::Equal
        }
    }

    pub fn bayes_factor_alt(&self) -> BayesFactor {
        BayesFactor::new(self.prob_alt, self.prob_ref)
    }

    pub fn bayes_factor_ref(&self) -> BayesFactor {
        BayesFactor::new(self.prob_ref, self.prob_alt)
    }

    pub fn prob_mapping_orig(&self) -> LogProb {
        self.prob_mapping
    }

    pub fn prob_mapping(&self) -> LogProb {
        self.prob_mapping_adj.unwrap_or(self.prob_mapping)
    }

    pub fn prob_mismapping(&self) -> LogProb {
        self.prob_mismapping_adj.unwrap_or(self.prob_mismapping)
    }

    pub fn prob_alt(&self) -> LogProb {
        self.prob_alt_adj.unwrap_or(self.prob_alt)
    }

    pub fn prob_ref(&self) -> LogProb {
        self.prob_ref_adj.unwrap_or(self.prob_ref)
    }

    pub fn prob_alt_orig(&self) -> LogProb {
        self.prob_alt
    }

    pub fn prob_ref_orig(&self) -> LogProb {
        self.prob_ref
    }

    pub fn is_uniquely_mapping(&self) -> bool {
        self.prob_mapping() >= *PROB_095
    }

    pub fn is_strong_alt_support(&self) -> bool {
        BayesFactor::new(self.prob_alt, self.prob_ref).evidence_kass_raftery()
            >= KassRaftery::Strong
    }

    pub fn is_strong_ref_support(&self) -> bool {
        BayesFactor::new(self.prob_ref, self.prob_alt).evidence_kass_raftery()
            >= KassRaftery::Strong
    }

    pub fn is_ref_support(&self) -> bool {
        self.prob_ref > self.prob_alt
    }

    pub fn is_positive_ref_support(&self) -> bool {
        BayesFactor::new(self.prob_ref, self.prob_alt).evidence_kass_raftery()
            >= KassRaftery::Positive
    }

    pub fn has_homopolymer_error(&self) -> bool {
        self.homopolymer_indel_len
            .map(|indel_len| indel_len != 0)
            .unwrap_or(false)
    }
}

impl ProcessedReadObservation {
    pub(crate) fn adjust_prob_mapping(
        pileup: &mut [Self],
        alignment_properties: &AlignmentProperties,
    ) {
        if !pileup.is_empty() {
            // METHOD: adjust MAPQ to get rid of stochastically inflated ones.
            // The idea here is to estimate the probability that this locus is
            // actually filled with reads from a homolog (e.g. induced by another variant somewhere else), say ph.
            // We do this by considering all non-max MAPQs as indicative of this.
            // Conservatively, we recalibrate those to a probability of 0.5.
            // Then, adjusted MAPQs are calculated as ph * 0.5 + 1-ph * max_mapq.
            // Technically, we just compute the mean here, which yields the same result.

            let max_prob_mapping =
                LogProb::from(PHREDProb(alignment_properties.max_mapq as f64)).ln_one_minus_exp();

            let probs = pileup
                .iter()
                .map(|obs| {
                    if relative_eq!(*obs.prob_mapping_orig(), *max_prob_mapping) {
                        obs.prob_mapping_orig()
                    } else {
                        *PROB_05
                    }
                })
                .collect_vec();

            let mut prob_sum = LogProb::ln_sum_exp(&probs);

            let calc_average = |prob_sum: LogProb, n| LogProb(*prob_sum - (n as f64).ln());
            let mut average = calc_average(prob_sum, pileup.len());
            if pileup.len() < 20 {
                // METHOD: for low depths, this method does not reliably work because it can be that by accident the
                // low MAPQ reads are not in the pileup. In order to correct for this sampling issue,
                // we add one pseudo low MAPQ observation. The higher the depth becomes, the less this observation
                // plays a role.
                prob_sum = prob_sum.ln_add_exp(*PROB_05);

                average = calc_average(prob_sum, pileup.len() + 1);
            }

            for obs in pileup {
                obs.prob_mapping_adj = Some(average);
                obs.prob_mismapping_adj = Some(average.ln_one_minus_exp());
            }
        }
    }
}

/// Adjusts probabilities for singleton evidence (variants supported by only one read
/// across all samples).
///
/// # Arguments
/// * `pileups` - Mutable slice of pileups to analyze and potentially adjust
///
/// # Returns
/// `true` if an adjustment was made (singleton evidence was found and adjusted)
///
/// # Side effects
/// When singleton evidence is found, both `prob_alt_adj` and `prob_ref_adj` of the observation
/// are set to 0.5 (PROB_05) to reflect the increased uncertainty that the variant could be
/// an artifact like a PCR error.
pub(crate) fn adjust_singleton_evidence(pileups: &mut [Pileup]) -> bool {
    let mut alt_observations: Vec<&mut ReadObservation<_, _>> = pileups
        .iter_mut()
        .flat_map(|pileup| pileup.read_observations_mut().iter_mut())
        .filter(|obs| obs.prob_alt > obs.prob_ref)
        .collect();
    if alt_observations.len() == 1 {
        let obs = alt_observations.first_mut().unwrap();
        obs.prob_alt_adj = Some(*PROB_05);
        obs.prob_ref_adj = Some(*PROB_05);
        true
    } else {
        false
    }
}

pub(crate) fn major_read_position(
    pileup: &[ReadObservation<Option<u32>, ExactAltLoci>],
) -> Option<u32> {
    calc_major_feature(pileup.iter().filter_map(|obs| {
        if obs.prob_alt > obs.prob_ref {
            obs.read_position
        } else {
            None
        }
    }))
}

pub(crate) fn major_alt_locus(
    pileup: &[ReadObservation<Option<u32>, ExactAltLoci>],
    alignment_properties: &AlignmentProperties,
) -> Option<genome::Locus> {
    calc_major_feature(
        pileup
            .iter()
            // TODO filter for alt obs only?
            //.filter(|obs| obs.prob_alt > obs.prob_ref)
            .flat_map(|obs| {
                obs.alt_locus
                    .iter()
                    .map(|locus| locus_to_bucket(locus, alignment_properties))
            }),
    )
}

/// Something that can be converted into observations.
pub(crate) trait ReadObservable: ReadVariant {
    fn extract_observations(
        &self,
        buffer: &mut sample::RecordBuffer,
        alignment_properties: &mut AlignmentProperties,
        max_depth: usize,
        alt_variants: &[Box<dyn Realignable>],
        observation_id_factory: &mut Option<&mut FragmentIdFactory>,
    ) -> Result<Vec<ReadObservation>>;

    /// Convert MAPQ (from read mapper) to LogProb for the event that the read maps
    /// correctly.
    fn prob_mapping(&self, evidence: &Evidence) -> LogProb;

    /// Return the minimum MAPQ of all records involved in the given evidence.
    fn min_mapq(&self, evidence: &Evidence) -> u8;

    /// Calculate an observation from the given evidence.
    fn evidence_to_observation(
        &self,
        evidence: &Evidence,
        alignment_properties: &mut AlignmentProperties,
        homopolymer_error_model: &Option<HomopolymerErrorModel>,
        alt_variants: &[Box<dyn Realignable>],
        observation_id_factory: &mut Option<&mut FragmentIdFactory>,
    ) -> Result<Option<ReadObservation>> {
        let id = observation_id_factory
            .as_mut()
            .map(|factory| factory.register(evidence));

        Ok(
            match self.allele_support(evidence, alignment_properties, alt_variants)? {
                // METHOD: for precise variants,
                // only consider allele support if it comes either from forward or reverse strand.
                // Unstranded observations (e.g. only insert size), are too unreliable, or do not contain
                // any information (e.g. no overlap).
                Some(allele_support)
                    if allele_support.strand() != Strand::None || self.is_imprecise() =>
                {
                    let alt_indel_len = allele_support.homopolymer_indel_len().unwrap_or(0);

                    let mut obs = ReadObservationBuilder::default();
                    obs.name(Some(evidence.id().to_string().to_owned()))
                        .fragment_id(id)
                        .prob_mapping_mismapping(self.prob_mapping(evidence))
                        .prob_alt(allele_support.prob_alt_allele())
                        .prob_ref(allele_support.prob_ref_allele())
                        .prob_sample_alt(self.prob_sample_alt(evidence, alignment_properties))
                        .prob_missed_allele(allele_support.prob_missed_allele())
                        .prob_overlap(if allele_support.strand() == Strand::Both {
                            LogProb::ln_one()
                        } else {
                            LogProb::ln_zero()
                        })
                        .strand(allele_support.strand())
                        .read_orientation(evidence.read_orientation()?)
                        .softclipped(evidence.softclipped())
                        .read_position(allele_support.read_position())
                        .paired(evidence.is_paired())
                        .prob_hit_base(LogProb::ln_one() - LogProb((evidence.len() as f64).ln()))
                        .is_max_mapq(self.min_mapq(evidence) == alignment_properties.max_mapq)
                        .alt_locus(evidence.alt_loci())
                        .third_allele_evidence(allele_support.third_allele_evidence().map(|d| *d));

                    if let Some(homopolymer_error_model) = homopolymer_error_model {
                        let ref_indel_len =
                            alt_indel_len + homopolymer_error_model.variant_homopolymer_indel_len();

                        obs.homopolymer_indel_len(Some(ref_indel_len));

                        if ref_indel_len == 0 || alt_indel_len == 0 {
                            // no homopolymer indel in read compared to reference
                            obs.prob_observable_at_homopolymer_artifact(None)
                                .prob_observable_at_homopolymer_variant(None);
                        } else {
                            obs.prob_observable_at_homopolymer_variant(Some(
                                homopolymer_error_model.prob_observable(alt_indel_len),
                            ))
                            .prob_observable_at_homopolymer_artifact(Some(
                                homopolymer_error_model.prob_observable(ref_indel_len),
                            ));
                        }
                    } else {
                        obs.homopolymer_indel_len(None)
                            .prob_observable_at_homopolymer_artifact(None)
                            .prob_observable_at_homopolymer_variant(None);
                    }

                    Some(obs.build().unwrap())
                }
                _ => None,
            },
        )
    }
}

/// An observation for or against a variant.
#[derive(Clone, Debug, Builder, Default, Serialize)]
pub struct ReadObservation<P = Option<u32>, A = ExactAltLoci>
where
    P: Clone,
    A: Clone,
{
    pub name: Option<String>,
    pub fragment_id: Option<u64>,
    /// Posterior probability that the read/read-pair has been mapped correctly (1 - MAPQ).
    prob_mapping: LogProb,
    /// Posterior probability that the read/read-pair has been mapped incorrectly (MAPQ).
    prob_mismapping: LogProb,
    /// Posterior probability that the read/read-pair has been mapped correctly (1 - MAPQ), adjusted form.
    #[builder(private, default = "None")]
    prob_mapping_adj: Option<LogProb>,
    /// Posterior probability that the read/read-pair has been mapped incorrectly (MAPQ), adjusted form.
    #[builder(private, default = "None")]
    prob_mismapping_adj: Option<LogProb>,
    /// Probability that the read/read-pair comes from the alternative allele.
    prob_alt: LogProb,
    /// Probability that the read/read-pair comes from the reference allele.
    prob_ref: LogProb,
    #[builder(private, default = "None")]
    prob_alt_adj: Option<LogProb>,
    #[builder(private, default = "None")]
    prob_ref_adj: Option<LogProb>,
    /// Probability that the read/read-pair comes from an unknown allele at an unknown true
    /// locus (in case it is mismapped). This should usually be set as the product of the maxima
    /// of prob_ref and prob_alt per read.
    pub prob_missed_allele: LogProb,
    /// Probability to sample the alt allele
    pub prob_sample_alt: LogProb,
    /// Probability to overlap with both strands
    #[builder(private)]
    pub prob_double_overlap: LogProb,
    /// Probability to overlap with one strand only (1-prob_double_overlap)
    #[builder(private)]
    pub prob_single_overlap: LogProb,
    pub prob_hit_base: LogProb,
    /// Strand evidence this observation relies on
    pub strand: Strand,
    /// Read orientation support this observation relies on
    pub read_orientation: SequenceReadPairOrientation,
    /// True if obervation contains s
    pub softclipped: bool,
    pub paired: bool,
    /// Read position of the variant in the read (for SNV and MNV)
    pub read_position: P,
    /// Probability to make this observation at a homopolymer artifact
    pub prob_observable_at_homopolymer_artifact: Option<LogProb>,
    pub prob_observable_at_homopolymer_variant: Option<LogProb>,
    /// Homopolymer indel length (None if there is no homopolymer indel compared to reference)
    pub homopolymer_indel_len: Option<i8>,
    pub is_max_mapq: bool,
    pub alt_locus: A,
    /// Edit distance of the read against the alt allele. Only recorded if it is higher than
    /// the expected number of sequencing errors of each type.
    pub third_allele_evidence: Option<u32>,
}
