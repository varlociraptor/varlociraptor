use std::hash::Hash;

use bio::stats::probs::LogProb;

use crate::utils::PROB_05;
use crate::variants::evidence::observation::{Observation, ReadPosition};
use crate::variants::model::bias::Bias;

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub(crate) enum HomopolymerError {
    None,
    Some,
}

impl HomopolymerError {
    pub(crate) fn values() -> Vec<Self> {
        vec![HomopolymerError::None, HomopolymerError::Some]
    }
}

impl Default for HomopolymerError {
    fn default() -> Self {
        HomopolymerError::None
    }
}

impl Bias for HomopolymerError {
    fn prob_alt(&self, observation: &Observation<ReadPosition>) -> LogProb {
        match (self, observation.homopolymer_indel_len) {
            (HomopolymerError::None, Some(len)) => {
                if len != 0 {
                    LogProb::ln_zero()
                } else {
                    LogProb::ln_one()
                }
            }
            (HomopolymerError::Some, Some(len)) => {
                if len != 0 {
                    observation.prob_homopolymer_error.unwrap()
                } else {
                    observation
                        .prob_homopolymer_error
                        .unwrap()
                        .ln_one_minus_exp()
                }
            }
            (HomopolymerError::None, None) => *PROB_05,
            (HomopolymerError::Some, None) => *PROB_05,
        }
    }

    fn prob_ref(&self, observation: &Observation<ReadPosition>) -> LogProb {
        self.prob_alt(observation)
    }

    fn prob_any(&self, _observation: &Observation<ReadPosition>) -> LogProb {
        LogProb::ln_one() // TODO check this
    }

    fn is_artifact(&self) -> bool {
        *self != HomopolymerError::None
    }

    fn is_informative(&self, pileups: &[Vec<Observation<ReadPosition>>]) -> bool {
        if !self.is_artifact() {
            return true;
        }
        // METHOD: this bias is only relevant if there is at least one recorded indel operation (indel operations are only recorded for some variants).
        pileups
            .iter()
            .any(|pileup| pileup.iter().any(|obs| self.is_bias_evidence(obs)))
    }

    fn is_possible(&self, pileups: &[Vec<Observation<ReadPosition>>]) -> bool {
        self.is_informative(pileups)
    }

    fn is_bias_evidence(&self, observation: &Observation<ReadPosition>) -> bool {
        observation.homopolymer_indel_len.unwrap_or(0) != 0
    }

    fn min_strong_evidence_ratio(&self) -> f64 {
        0.01 // TODO possibly tune, this is very permissive
    }
}
