use std::cmp;

use bio::stats::probs::LogProb;

use itertools::Itertools;
use strum::IntoEnumIterator;

use crate::variants::evidence::observations::{
    pileup::Pileup, read_observation::ProcessedReadObservation,
};

pub(crate) mod alt_locus_bias;
pub(crate) mod homopolymer_error;
pub(crate) mod read_orientation_bias;
pub(crate) mod read_position_bias;
pub(crate) mod softclip_bias;
pub(crate) mod strand_bias;

pub(crate) use alt_locus_bias::AltLocusBias;
pub(crate) use homopolymer_error::HomopolymerError;
pub(crate) use read_orientation_bias::ReadOrientationBias;
pub(crate) use read_position_bias::ReadPositionBias;
pub(crate) use softclip_bias::SoftclipBias;
pub(crate) use strand_bias::StrandBias;

pub(crate) trait Bias: Default + cmp::PartialEq + std::fmt::Debug {
    fn prob_alt(&self, observation: &ProcessedReadObservation) -> LogProb;

    fn prob_ref(&self, observation: &ProcessedReadObservation) -> LogProb {
        self.prob_any(observation)
    }

    fn prob_any(&self, observation: &ProcessedReadObservation) -> LogProb;

    fn is_artifact(&self) -> bool;

    fn is_possible(&self, pileups: &[Pileup]) -> bool {
        if !self.is_artifact() {
            return true;
        }

        pileups.iter().any(|pileup| {
            pileup
                .read_observations()
                .iter()
                .any(|observation| self.prob_alt(observation) != LogProb::ln_zero())
        })
    }

    fn is_informative(&self, _pileups: &[Pileup]) -> bool {
        true
    }

    fn is_bias_evidence(&self, observation: &ProcessedReadObservation) -> bool {
        self.prob_alt(observation) != LogProb::ln_zero()
    }

    fn min_strong_evidence_ratio(&self) -> f64 {
        0.66666
    }

    fn is_likely(&self, pileups: &[Pileup]) -> bool {
        if !self.is_artifact() {
            true
        } else {
            pileups.iter().any(|pileup| {
                let strong_all = pileup
                    .read_observations()
                    .iter()
                    .filter(|obs| obs.is_uniquely_mapping() && obs.is_strong_alt_support())
                    .count();
                if strong_all >= 10 {
                    let strong_bias_evidence = pileup
                        .read_observations()
                        .iter()
                        .filter(|obs| {
                            obs.is_uniquely_mapping()
                                && obs.is_strong_alt_support()
                                && self.is_bias_evidence(obs)
                        })
                        .count();
                    // METHOD: there is bias evidence if we have at least two third of the strong observations supporting the bias
                    let ratio = strong_bias_evidence as f64 / strong_all as f64;
                    ratio >= self.min_strong_evidence_ratio()
                } else if pileup
                    .read_observations()
                    .iter()
                    .all(|obs| obs.is_ref_support())
                {
                    // METHOD: if all obs are towards REF allele, there is no need to consider biases.
                    // The variant will anyway be called as absent.
                    // This can safe a lot of time and also avoids unexpected reporting
                    // of artifacts in ambiguous cases.
                    false
                } else if pileup.read_observations().is_empty() {
                    // METHOD: no reads, no need to consider for biases, hence, skip with false
                    false
                } else {
                    // METHOD: not enough reads, rather consider all biases to be sure
                    true
                }
            })
        }
    }

    /// Learn parameters needed for estimation on current pileup.
    fn learn_parameters(&mut self, _pileups: &[Pileup]) {
        // METHOD: by default, there is nothing to learn, however, a bias can use this to
        // infer some parameters over which we would otherwise need to integrate (which would hamper
        // performance too much).
    }
}

#[derive(Builder, CopyGetters, Getters, Debug, Clone, Eq, PartialEq, Hash)]
pub(crate) struct Artifacts {
    #[getset(get = "pub(crate)")]
    strand_bias: StrandBias,
    #[getset(get = "pub(crate)")]
    read_orientation_bias: ReadOrientationBias,
    #[getset(get = "pub(crate)")]
    read_position_bias: ReadPositionBias,
    #[getset(get = "pub(crate)")]
    softclip_bias: SoftclipBias,
    #[getset(get = "pub(crate)")]
    homopolymer_error: HomopolymerError,
    #[getset(get = "pub(crate)")]
    alt_locus_bias: AltLocusBias,
}

impl Artifacts {
    pub(crate) fn all_artifact_combinations(
        consider_read_orientation_bias: bool,
        consider_strand_bias: bool,
        consider_read_position_bias: bool,
        consider_softclip_bias: bool,
        consider_homopolymer_error: bool,
        consider_alt_locus_bias: bool,
    ) -> Box<dyn Iterator<Item = Self>> {
        if !consider_strand_bias
            && !consider_read_orientation_bias
            && !consider_read_position_bias
            && !consider_softclip_bias
            && !consider_homopolymer_error
            && !consider_alt_locus_bias
        {
            return Box::new(std::iter::empty());
        }

        let strand_biases = if consider_strand_bias {
            StrandBias::iter().collect_vec()
        } else {
            vec![StrandBias::default()]
        };
        let read_position_biases = if consider_read_position_bias {
            ReadPositionBias::iter().collect_vec()
        } else {
            vec![ReadPositionBias::None]
        };
        let read_orientation_biases = if consider_read_orientation_bias {
            ReadOrientationBias::iter().collect_vec()
        } else {
            vec![ReadOrientationBias::None]
        };
        let softclip_biases = if consider_softclip_bias {
            SoftclipBias::iter().collect_vec()
        } else {
            vec![SoftclipBias::None]
        };
        let homopolymer_error = if consider_homopolymer_error {
            HomopolymerError::values()
        } else {
            vec![HomopolymerError::default()]
        };
        let alt_locus_bias = if consider_alt_locus_bias {
            AltLocusBias::iter().collect_vec()
        } else {
            vec![AltLocusBias::default()]
        };

        Box::new(
            strand_biases
                .into_iter()
                .cartesian_product(read_orientation_biases.into_iter())
                .cartesian_product(read_position_biases.into_iter())
                .cartesian_product(softclip_biases.into_iter())
                .cartesian_product(homopolymer_error.into_iter())
                .cartesian_product(alt_locus_bias.into_iter())
                .filter_map(|(((((sb, rob), rpb), scb), dib), alb)| {
                    if [
                        sb.is_artifact(),
                        rob.is_artifact(),
                        rpb.is_artifact(),
                        scb.is_artifact(),
                        dib.is_artifact(),
                        alb.is_artifact(),
                    ]
                    .iter()
                    .map(|artifact| if *artifact { 1 } else { 0 })
                    .sum::<usize>()
                        == 1
                    {
                        Some(
                            ArtifactsBuilder::default()
                                .strand_bias(sb)
                                .read_orientation_bias(rob)
                                .read_position_bias(rpb)
                                .softclip_bias(scb)
                                .homopolymer_error(dib)
                                .alt_locus_bias(alb)
                                .build()
                                .unwrap(),
                        )
                    } else {
                        None
                    }
                }),
        )
    }

    pub(crate) fn none() -> Self {
        ArtifactsBuilder::default()
            .strand_bias(StrandBias::default())
            .read_orientation_bias(ReadOrientationBias::None)
            .read_position_bias(ReadPositionBias::None)
            .softclip_bias(SoftclipBias::None)
            .homopolymer_error(HomopolymerError::default())
            .alt_locus_bias(AltLocusBias::None)
            .build()
            .unwrap()
    }

    pub(crate) fn is_possible(&self, pileups: &[Pileup]) -> bool {
        self.strand_bias.is_possible(pileups)
            && self.read_orientation_bias.is_possible(pileups)
            && self.read_position_bias.is_possible(pileups)
            && self.softclip_bias.is_possible(pileups)
            && self.homopolymer_error.is_possible(pileups)
            && self.alt_locus_bias.is_possible(pileups)
    }

    pub(crate) fn is_informative(&self, pileups: &[Pileup]) -> bool {
        self.strand_bias.is_informative(pileups)
            && self.read_orientation_bias.is_informative(pileups)
            && self.read_position_bias.is_informative(pileups)
            && self.softclip_bias.is_informative(pileups)
            && self.homopolymer_error.is_informative(pileups)
            && self.alt_locus_bias.is_informative(pileups)
    }

    pub(crate) fn is_likely(&self, pileups: &[Pileup]) -> bool {
        self.strand_bias.is_likely(pileups)
            && self.read_orientation_bias.is_likely(pileups)
            && self.read_position_bias.is_likely(pileups)
            && self.softclip_bias.is_likely(pileups)
            && self.homopolymer_error.is_likely(pileups)
            && self.alt_locus_bias.is_likely(pileups)
    }

    pub(crate) fn prob_alt(&self, observation: &ProcessedReadObservation) -> LogProb {
        self.strand_bias.prob_alt(observation)
            + self.read_orientation_bias.prob_alt(observation)
            + self.read_position_bias.prob_alt(observation)
            + self.softclip_bias.prob_alt(observation)
            + self.homopolymer_error.prob_alt(observation)
            + self.alt_locus_bias.prob_alt(observation)
    }

    pub(crate) fn prob_ref(&self, observation: &ProcessedReadObservation) -> LogProb {
        self.strand_bias.prob_ref(observation)
            + self.read_orientation_bias.prob_ref(observation)
            + self.read_position_bias.prob_ref(observation)
            + self.softclip_bias.prob_ref(observation)
            + self.homopolymer_error.prob_ref(observation)
            + self.alt_locus_bias.prob_ref(observation)
    }

    pub(crate) fn prob_any(&self, observation: &ProcessedReadObservation) -> LogProb {
        self.strand_bias.prob_any(observation)
            + self.read_orientation_bias.prob_any(observation)
            + self.read_position_bias.prob_any(observation)
            + self.softclip_bias.prob_any(observation)
            + self.homopolymer_error.prob_any(observation)
            + self.alt_locus_bias.prob_any(observation)
    }

    pub(crate) fn is_artifact(&self) -> bool {
        self.strand_bias.is_artifact()
            || self.read_orientation_bias.is_artifact()
            || self.read_position_bias.is_artifact()
            || self.softclip_bias.is_artifact()
            || self.homopolymer_error.is_artifact()
            || self.alt_locus_bias.is_artifact()
    }

    pub(crate) fn learn_parameters(&mut self, pileups: &[Pileup]) {
        self.homopolymer_error.learn_parameters(pileups);
        self.strand_bias.learn_parameters(pileups);
    }
}
