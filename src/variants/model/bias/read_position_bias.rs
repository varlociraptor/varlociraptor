use bio::stats::probs::LogProb;
use itertools::Itertools;

use crate::variants::evidence::observations::observation::ReadPosition;
use crate::variants::evidence::observations::pileup::Pileup;
use crate::variants::evidence::observations::read_observation::ProcessedReadObservation;
use crate::variants::model::bias::Bias;

#[derive(Copy, Clone, Default, PartialOrd, PartialEq, Eq, Debug, Ord, EnumIter, Hash)]
pub(crate) enum ReadPositionBias {
    #[default]
    None,
    Some,
}

impl Bias for ReadPositionBias {
    fn prob_alt(&self, observation: &ProcessedReadObservation) -> LogProb {
        match (self, observation.read_position) {
            (ReadPositionBias::None, ReadPosition::Major) => observation.prob_hit_base,
            (ReadPositionBias::None, ReadPosition::Some) => {
                Self::one_minus_prob_hit_base(observation)
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
            ReadPosition::Major => observation.prob_hit_base,
            ReadPosition::Some => Self::one_minus_prob_hit_base(observation),
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
    fn one_minus_prob_hit_base(observation: &ProcessedReadObservation) -> LogProb {
        if observation.prob_hit_base != LogProb::ln_one() {
            observation.prob_hit_base.ln_one_minus_exp()
        } else {
            // METHOD: read has length 1. Hence, prob_hit_base is 1.0.
            // We cannot distinguish between the major and minor position.
            // There is only one possible position, and the probabilit for
            // that is always 1.0.
            LogProb::ln_one()
        }
    }

    fn has_valid_major_rate(pileups: &[Pileup]) -> bool {
        // TODO what happens when we combine amplicon and WGS here?
        // The former might give a hint for a bias because of the way amplicons
        // are designed, while the latter might not have a bias.
        pileups.iter().any(|pileup| {
            let expected_all = LogProb::ln_sum_exp(
                &pileup
                    .read_observations()
                    .iter()
                    .filter_map(|obs| {
                        if obs.is_strong_ref_support() {
                            Some(obs.prob_mapping())
                        } else {
                            None
                        }
                    })
                    .collect_vec(),
            )
            .exp();

            if expected_all > 10.0 {
                let expected_major = LogProb::ln_sum_exp(
                    &pileup
                        .read_observations()
                        .iter()
                        .filter_map(|obs| {
                            if obs.is_strong_ref_support()
                                && obs.read_position == ReadPosition::Major
                            {
                                Some(obs.prob_mapping())
                            } else {
                                None
                            }
                        })
                        .collect_vec(),
                )
                .exp();

                let expected_major_rate = LogProb::ln_sum_exp(
                    &pileup
                        .read_observations()
                        .iter()
                        .filter_map(|obs| {
                            if obs.is_strong_ref_support() {
                                Some(obs.prob_mapping() + obs.prob_hit_base)
                            } else {
                                None
                            }
                        })
                        .collect_vec(),
                )
                .exp();

                let major_rate = expected_major / expected_all;
                expected_major > 0.0 && (major_rate - expected_major_rate).abs() < 0.05
            } else {
                false
            }
        })
    }
}
