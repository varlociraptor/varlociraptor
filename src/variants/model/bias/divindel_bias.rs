use std::cmp;

use bio::stats::probs::LogProb;
use ordered_float::NotNan;

use crate::variants::evidence::observation::{Observation, ReadPosition};
use crate::variants::model::bias::Bias;

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord, Hash)]
pub(crate) enum DivIndelBias {
    None,
    Some {
        other_rate: NotNan<f64>,
        min_other_rate: NotNan<f64>,
    },
}

impl DivIndelBias {
    pub(crate) fn values(min_other_rate: f64) -> Vec<Self> {
        vec![
            DivIndelBias::None,
            DivIndelBias::Some {
                other_rate: NotNan::new(0.0).unwrap(),
                min_other_rate: NotNan::new(min_other_rate).unwrap(),
            },
        ]
    }
}

impl Default for DivIndelBias {
    fn default() -> Self {
        DivIndelBias::None
    }
}

impl Bias for DivIndelBias {
    fn prob(&self, observation: &Observation<ReadPosition>) -> LogProb {
        match self {
            DivIndelBias::None => {
                if observation.has_alt_indel_operations {
                    LogProb::ln_zero()
                } else {
                    LogProb::ln_one()
                }
            }
            DivIndelBias::Some { other_rate, .. } => {
                if **other_rate == 0.0 {
                    // METHOD: if there are no other operations there is no artifact.
                    LogProb::ln_zero()
                } else {
                    if observation.has_alt_indel_operations {
                        LogProb(other_rate.ln())
                    } else {
                        LogProb((1.0 - **other_rate).ln())
                    }
                }
            }
        }
    }

    fn prob_any(&self, _observation: &Observation<ReadPosition>) -> LogProb {
        LogProb::ln_one() // TODO check this
    }

    fn is_artifact(&self) -> bool {
        *self != DivIndelBias::None
    }

    fn is_informative(&self, pileups: &[Vec<Observation<ReadPosition>>]) -> bool {
        if !self.is_artifact() {
            return true;
        }
        // METHOD: this bias is only relevant if there is at least one recorded indel operation (indel operations are only recorded for some variants).
        pileups
            .iter()
            .any(|pileup| pileup.iter().any(|obs| obs.has_alt_indel_operations))
    }

    fn is_possible(&self, pileups: &[Vec<Observation<ReadPosition>>]) -> bool {
        if !self.is_artifact() {
            return true;
        }

        pileups.iter().any(|pileup| {
            pileup.iter().any(|observation| match self {
                DivIndelBias::Some { .. } => observation.has_alt_indel_operations,
                DivIndelBias::None => self.prob(observation) != LogProb::ln_zero(),
            })
        })
    }

    fn is_bias_evidence(&self, observation: &Observation<ReadPosition>) -> bool {
        observation.has_alt_indel_operations
    }

    fn min_strong_evidence_ratio(&self) -> f64 {
        if let DivIndelBias::Some { other_rate, .. } = self {
            0.66666 * **other_rate
        } else {
            unreachable!();
        }
    }

    fn learn_parameters(&mut self, pileups: &[Vec<Observation<ReadPosition>>]) {
        // METHOD: by default, there is nothing to learn, however, a bias can use this to
        // infer some parameters over which we would otherwise need to integrate (which would hamper
        // performance too much).
        if let DivIndelBias::Some {
            ref mut other_rate,
            min_other_rate,
        } = self
        {
            let strong_all = pileups
                .iter()
                .map(|pileup| pileup.iter().filter(&Self::is_strong_obs))
                .flatten()
                .count();
            let strong_other = pileups
                .iter()
                .map(|pileup| {
                    pileup
                        .iter()
                        .filter(|obs| Self::is_strong_obs(obs) && obs.has_alt_indel_operations)
                })
                .flatten()
                .count();

            let rate = NotNan::new(if strong_all > 0 {
                strong_other as f64 / strong_all as f64
            } else {
                0.0
            })
            .unwrap();
            dbg!((rate, &min_other_rate));

            *other_rate = cmp::max(rate, *min_other_rate);
        }
    }
}
