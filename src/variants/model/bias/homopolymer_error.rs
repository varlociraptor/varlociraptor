use std::hash::Hash;

use bio::stats::probs::LogProb;

use crate::variants::evidence::observations::pileup::Pileup;
use crate::variants::evidence::observations::read_observation::ProcessedReadObservation;
use crate::variants::model::bias::Bias;

#[derive(Clone, Debug, Hash, PartialEq, Eq, Default)]
pub(crate) enum HomopolymerError {
    #[default]
    None,
    Some,
}

impl HomopolymerError {
    pub(crate) fn values() -> Vec<Self> {
        vec![HomopolymerError::None, HomopolymerError::Some]
    }
}

impl Bias for HomopolymerError {
    fn prob_alt(&self, observation: &ProcessedReadObservation) -> LogProb {
        match self {
            HomopolymerError::Some => observation
                .prob_observable_at_homopolymer_artifact
                .unwrap_or(LogProb::ln_one()),
            HomopolymerError::None => observation
                .prob_observable_at_homopolymer_variant
                .unwrap_or(LogProb::ln_one()),
        }
    }

    fn prob_ref(&self, observation: &ProcessedReadObservation) -> LogProb {
        self.prob_alt(observation)
    }

    fn prob_any(&self, _observation: &ProcessedReadObservation) -> LogProb {
        LogProb::ln_one() // TODO check this
    }

    fn is_artifact(&self) -> bool {
        *self != HomopolymerError::None
    }

    fn is_informative(&self, pileups: &[Pileup]) -> bool {
        if !self.is_artifact() {
            return true;
        }
        // METHOD: we require all alt supporting samples to have at least one
        // homopolymer indel in both directions relative to the ref allele.
        // Otherwise, we can assume that it is rather not a homopolymer error
        // because it seems to rather support an indel in one particular direction.
        pileups.iter().all(|pileup| {
            let has_homopolymer_indel = |ins: bool| {
                pileup.read_observations().iter().any(|obs| {
                    let indel = obs.homopolymer_indel_len.unwrap_or(0);
                    if ins {
                        indel > 0
                    } else {
                        indel < 0
                    }
                })
            };

            !pileup
                .read_observations()
                .iter()
                .any(|obs| obs.is_strong_alt_support())
                || (has_homopolymer_indel(true) && has_homopolymer_indel(false))
        })
    }

    fn is_possible(&self, pileups: &[Pileup]) -> bool {
        self.is_informative(pileups)
    }

    fn is_likely(&self, pileups: &[Pileup]) -> bool {
        self.is_informative(pileups)
    }

    fn is_bias_evidence(&self, observation: &ProcessedReadObservation) -> bool {
        observation.homopolymer_indel_len.unwrap_or(0) != 0
    }

    fn min_strong_evidence_ratio(&self) -> f64 {
        0.01 // TODO possibly tune, this is very permissive
    }
}
