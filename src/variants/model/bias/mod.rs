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
    fn artifact_values() -> Vec<Self>;

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

#[derive(TypedBuilder, CopyGetters, Getters, Debug, Clone, Eq, PartialEq, Hash)]
pub struct Artifacts {
    #[builder(default)]
    strand_bias: Option<StrandBias>,
    #[builder(default)]
    read_orientation_bias: Option<ReadOrientationBias>,
    #[builder(default)]
    read_position_bias: Option<ReadPositionBias>,
    #[builder(default)]
    softclip_bias: Option<SoftclipBias>,
    #[builder(default)]
    homopolymer_error:Option<HomopolymerError>,
    #[builder(default)]
    alt_locus_bias: Option<AltLocusBias>,
}

impl Artifacts {
    pub(crate) fn all_artifact_combinations(
        consider_read_orientation_bias: bool,
        consider_strand_bias: bool,
        consider_read_position_bias: bool,
        consider_softclip_bias: bool,
        consider_homopolymer_error: bool,
        consider_alt_locus_bias: bool,
    ) -> Vec<Self> {
        let mut artifacts: Vec<Artifacts> = Vec::new();

        if consider_strand_bias {
            for bias in StrandBias::artifact_values() {
                artifacts.push(Artifacts::builder().strand_bias(Some(bias)).build());
            }
        }

        if consider_read_position_bias {
            for bias in ReadPositionBias::artifact_values() {
                artifacts.push(Artifacts::builder().read_position_bias(Some(bias)).build());
            }
        }

        if consider_read_orientation_bias {
            for bias in ReadOrientationBias::artifact_values() {
                artifacts.push(Artifacts::builder().read_orientation_bias(Some(bias)).build());
            }
        }

        if consider_softclip_bias {
            for bias in SoftclipBias::artifact_values() {
                artifacts.push(Artifacts::builder().softclip_bias(Some(bias)).build());
            }
        }

        if consider_homopolymer_error {
            for bias in HomopolymerError::artifact_values() {
                artifacts.push(Artifacts::builder().homopolymer_error(Some(bias)).build());
            }
        }

        if consider_alt_locus_bias {
            for bias in AltLocusBias::artifact_values() {
                artifacts.push(Artifacts::builder().alt_locus_bias(Some(bias)).build());
            }
        }
        artifacts
    }

    

    pub(crate) fn none(
        consider_read_orientation_bias: bool,
        consider_strand_bias: bool,
        consider_read_position_bias: bool,
        consider_softclip_bias: bool,
        consider_homopolymer_error: bool,
        consider_alt_locus_bias: bool,
    ) -> Self {
        Artifacts::builder()
            .strand_bias(if consider_strand_bias { Some(StrandBias::default()) } else { None })
            .read_orientation_bias(if consider_read_orientation_bias {
                Some(ReadOrientationBias::None)
            } else {
                None
            })
            .read_position_bias(if consider_read_position_bias {
                Some(ReadPositionBias::None)
            } else {
                None
            })
            .softclip_bias(if consider_softclip_bias {
                Some(SoftclipBias::None)
            } else {
                None
            })
            .homopolymer_error(if consider_homopolymer_error {
                Some(HomopolymerError::None)
            } else {
                None
            })
            .alt_locus_bias(if consider_alt_locus_bias {
                Some(AltLocusBias::None)
            } else {
                None
            })
            .build()
    }

    pub(crate) fn is_possible(&self, pileups: &[Pileup]) -> bool {
        self.strand_bias.map_or(true, |bias| bias.is_possible(pileups))
            && self.read_orientation_bias.map_or(true, |bias| bias.is_possible(pileups))
            && self.read_position_bias.map_or(true, |bias| bias.is_possible(pileups))
            && self.softclip_bias.map_or(true, |bias| bias.is_possible(pileups))
            && self.homopolymer_error.map_or(true, |bias| bias.is_possible(pileups))
            && self.alt_locus_bias.map_or(true, |bias| bias.is_possible(pileups))
    }

    pub(crate) fn is_informative(&self, pileups: &[Pileup]) -> bool {
        self.strand_bias.map_or(true, |bias| bias.is_informative(pileups))
            && self.read_orientation_bias.map_or(true, |bias| bias.is_informative(pileups))
            && self.read_position_bias.map_or(true, |bias| bias.is_informative(pileups))
            && self.softclip_bias.map_or(true, |bias| bias.is_informative(pileups))
            && self.homopolymer_error.map_or(true, |bias| bias.is_informative(pileups))
            && self.alt_locus_bias.map_or(true, |bias| bias.is_informative(pileups))
    }

    pub(crate) fn is_likely(&self, pileups: &[Pileup]) -> bool {
        self.strand_bias.map_or(true, |bias| bias.is_likely(pileups))
            && self.read_orientation_bias.map_or(true, |bias| bias.is_likely(pileups))
            && self.read_position_bias.map_or(true, |bias| bias.is_likely(pileups))
            && self.softclip_bias.map_or(true, |bias| bias.is_likely(pileups))
            && self.homopolymer_error.map_or(true, |bias| bias.is_likely(pileups))
            && self.alt_locus_bias.map_or(true, |bias| bias.is_likely(pileups))
    }

    pub(crate) fn prob_alt(&self, observation: &ProcessedReadObservation) -> LogProb {
        self.strand_bias.map_or(LogProb::ln_one(), |bias| bias.prob_alt(observation))
            + self.read_orientation_bias.map_or(LogProb::ln_one(), |bias| bias.prob_alt(observation))
            + self.read_position_bias.map_or(LogProb::ln_one(), |bias| bias.prob_alt(observation))
            + self.softclip_bias.map_or(LogProb::ln_one(), |bias| bias.prob_alt(observation))
            + self.homopolymer_error.map_or(LogProb::ln_one(), |bias| bias.prob_alt(observation))
            + self.alt_locus_bias.map_or(LogProb::ln_one(), |bias| bias.prob_alt(observation))
    }

    pub(crate) fn prob_ref(&self, observation: &ProcessedReadObservation) -> LogProb {
        self.strand_bias.map_or(LogProb::ln_one(), |bias| bias.prob_ref(observation))
            + self.read_orientation_bias.map_or(LogProb::ln_one(), |bias| bias.prob_ref(observation))
            + self.read_position_bias.map_or(LogProb::ln_one(), |bias| bias.prob_ref(observation))
            + self.softclip_bias.map_or(LogProb::ln_one(), |bias| bias.prob_ref(observation))
            + self.homopolymer_error.map_or(LogProb::ln_one(), |bias| bias.prob_ref(observation))
            + self.alt_locus_bias.map_or(LogProb::ln_one(), |bias| bias.prob_ref(observation))
    }

    pub(crate) fn prob_any(&self, observation: &ProcessedReadObservation) -> LogProb {
        self.strand_bias.map_or(LogProb::ln_one(), |bias| bias.prob_any(observation))
            + self.read_orientation_bias.map_or(LogProb::ln_one(), |bias| bias.prob_any(observation))
            + self.read_position_bias.map_or(LogProb::ln_one(), |bias| bias.prob_any(observation))
            + self.softclip_bias.map_or(LogProb::ln_one(), |bias| bias.prob_any(observation))
            + self.homopolymer_error.map_or(LogProb::ln_one(), |bias| bias.prob_any(observation))
            + self.alt_locus_bias.map_or(LogProb::ln_one(), |bias| bias.prob_any(observation))
    }

    pub(crate) fn is_artifact(&self) -> bool {
        self.strand_bias.map_or(false, |bias| bias.is_artifact())
            || self.read_orientation_bias.map_or(false, |bias| bias.is_artifact())
            || self.read_position_bias.map_or(false, |bias| bias.is_artifact())
            || self.softclip_bias.map_or(false, |bias| bias.is_artifact())
            || self.homopolymer_error.map_or(false, |bias| bias.is_artifact())
            || self.alt_locus_bias.map_or(false, |bias| bias.is_artifact())
    }

    pub(crate) fn learn_parameters(&mut self, pileups: &[Pileup]) {
        if let Some(ref mut bias) = self.homopolymer_error {
            bias.learn_parameters(pileups);
        }
        if let Some(ref mut bias) = self.strand_bias {
            bias.learn_parameters(pileups);
        }
    }

    pub(crate) fn strand_bias(&self) -> StrandBias {
        self.strand_bias.unwrap_or_default()
    }

    pub(crate) fn read_orientation_bias(&self) -> ReadOrientationBias {
        self.read_orientation_bias.unwrap_or_default()
    }

    pub(crate) fn read_position_bias(&self) -> ReadPositionBias {
        self.read_position_bias.unwrap_or_default()
    }

    pub(crate) fn softclip_bias(&self) -> SoftclipBias {
        self.softclip_bias.unwrap_or_default()
    }

    pub(crate) fn homopolymer_error(&self) -> HomopolymerError {
        self.homopolymer_error.unwrap_or_default()
    }

    pub(crate) fn alt_locus_bias(&self) -> AltLocusBias {
        self.alt_locus_bias.unwrap_or_default()
    }
}
