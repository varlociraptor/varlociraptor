use bio::stats::probs::LogProb;
use itertools::Itertools;

use crate::variants::evidence::observations::pileup::Pileup;
use crate::variants::evidence::observations::read_observation::{
    ProcessedReadObservation, ReadPosition,
};
use crate::variants::model::bias::Bias;

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord, EnumIter, Hash)]
pub(crate) enum ReadPositionBias {
    None,
    Some,
}

impl Default for ReadPositionBias {
    fn default() -> Self {
        ReadPositionBias::None
    }
}

impl Bias for ReadPositionBias {
    fn artifact_values() -> Vec<Self> {
        vec![ReadPositionBias::Some]
    }

    fn prob_alt(&self, observation: &ProcessedReadObservation) -> LogProb {
        match (self, observation.read_position) {
            (ReadPositionBias::None, ReadPosition::Major) => {
                observation.prob_hit_base
            }
            (ReadPositionBias::None, ReadPosition::Some) => {
                observation.prob_hit_base.ln_one_minus_exp()
            }
            (ReadPositionBias::Some, ReadPosition::Major) => LogProb::ln_one(), // bias
            (ReadPositionBias::Some, ReadPosition::Some) => LogProb::ln_zero(), // no bias
        }
    }

    fn prob_any(&self, observation: &ProcessedReadObservation) -> LogProb {
        // METHOD: It is important that the probability of the read position bias is
        // is the same for both alt and ref reads in case the bias is None.
        // Otherwise, the model can be drawn to wrong AF estimates.
        match observation.read_position {
            ReadPosition::Major => {
                observation.prob_hit_base
            }
            ReadPosition::Some => {
                observation.prob_hit_base.ln_one_minus_exp()
            }
        }
    }

    fn is_artifact(&self) -> bool {
        *self != ReadPositionBias::None
    }

    fn is_informative(&self, pileups: &[Pileup]) -> bool {
        // METHOD: if all reads overlap the variant at the major pos,
        // we cannot estimate a read position bias and None is the only informative one.
        !self.is_artifact() || Self::has_valid_major_rate(pileups)
    }
}

impl ReadPositionBias {
    fn has_valid_major_rate(pileups: &[Pileup]) -> bool {
        let strong_all = LogProb::ln_sum_exp(
            &pileups
                .iter()
                .flat_map(|pileup| {
                    pileup.read_observations().iter().filter_map(|obs| {
                        if obs.is_strong_ref_support() {
                            Some(obs.prob_mapping())
                        } else {
                            None
                        }
                    })
                })
                .collect_vec(),
        )
        .exp();
        let strong_major = LogProb::ln_sum_exp(
            &pileups
                .iter()
                .flat_map(|pileup| {
                    pileup.read_observations().iter().filter_map(|obs| {
                        if obs.is_strong_ref_support() && obs.read_position == ReadPosition::Major {
                            Some(obs.prob_mapping())
                        } else {
                            None
                        }
                    })
                })
                .collect_vec(),
        )
        .exp();
        let any_major = pileups.iter().any(|pileup| {
            pileup
                .read_observations()
                .iter()
                .any(|obs| obs.read_position == ReadPosition::Major)
        });

        if strong_all >= 10.0 {
            let major_fraction = strong_major / strong_all;
            // METHOD: if the major read position is rare in the ref supporting reads,
            // we consider read position bias. Otherwise, reads are somehow
            // all starting at the same position, e.g. caused by amplicon sequencing.
            major_fraction < 0.1
        } else {
            // METHOD: not enough data to say anything, let the model decide
            any_major
        }
    }
}
