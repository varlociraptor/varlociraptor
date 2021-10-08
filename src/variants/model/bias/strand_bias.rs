use bio::stats::Prob;
use bio::stats::probs::LogProb;
use ordered_float::NotNan;

use crate::utils::PROB_05;
use crate::variants::evidence::observation::{Observation, ReadPosition, Strand};
use crate::variants::model::bias::Bias;

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Debug, Ord, EnumIter, Hash)]
pub(crate) enum StrandBias {
    None {
        forward_rate: NotNan<f64>,
    },
    Forward,
    Reverse,
}

impl Default for StrandBias {
    fn default() -> Self {
        StrandBias::None { forward_rate: NotNan::new(0.5).unwrap() }
    }
}

impl Bias for StrandBias {
    fn prob(&self, observation: &Observation<ReadPosition>) -> LogProb {
        match (self, observation.strand) {
            (StrandBias::Forward, Strand::Forward) => LogProb::ln_one(),
            (StrandBias::Reverse, Strand::Forward) => LogProb::ln_zero(),
            (StrandBias::Forward, Strand::Reverse) => LogProb::ln_zero(),
            (StrandBias::Reverse, Strand::Reverse) => LogProb::ln_one(),
            (StrandBias::Forward, Strand::Both) => LogProb::ln_zero(),
            (StrandBias::Reverse, Strand::Both) => LogProb::ln_zero(),
            (StrandBias::None { .. }, Strand::Both) => observation.prob_double_overlap,
            (_, Strand::None) => unreachable!(),
            (StrandBias::None { forward_rate }, observed) => {
                let rate = match observed {
                    Strand::Forward => **forward_rate,
                    Strand::Reverse => 1.0 - **forward_rate,
                    Strand::Both | Strand::None => unreachable!(),
                };
                LogProb::from(Prob(rate)) + observation.prob_single_overlap
            },
        }
    }

    fn prob_any(&self, _observation: &Observation<ReadPosition>) -> LogProb {
        *PROB_05
    }

    fn is_artifact(&self) -> bool {
        if let StrandBias::None { .. } = self {
            false
        } else {
            true
        }
    }

    fn learn_parameters(&mut self, pileups: &[Vec<Observation<ReadPosition>>]) {
        let estimate_forward_rate = || {
            let strong_all = pileups
                .iter()
                .map(|pileup| pileup.iter().filter(|obs| obs.is_strong_ref_support()))
                .flatten()
                .count();
            let strong_forward = pileups
                .iter()
                .map(|pileup| {
                    pileup
                        .iter()
                        .filter(|obs| obs.is_strong_ref_support() && obs.strand == Strand::Forward)
                })
                .flatten()
                .count();
            if strong_all > 10 && strong_forward > 0 && strong_forward != strong_all {
                NotNan::new(strong_forward as f64 / strong_all as f64).unwrap()
            } else {
                // If there are not enough strong ref reads, we just assume the default 0.5.
                NotNan::new(0.5).unwrap()
            }
        };
        if let StrandBias::None { ref mut forward_rate } = self {
            *forward_rate = estimate_forward_rate();
        }
    }
}