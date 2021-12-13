use bio::stats::probs::LogProb;
use bio::stats::Prob;

use ordered_float::NotNan;

use crate::utils::PROB_05;
use crate::variants::evidence::observation::{Observation, ReadPosition, Strand};
use crate::variants::model::bias::Bias;

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Debug, Ord, EnumIter, Hash)]
pub(crate) enum StrandBias {
    None { forward_rate: NotNan<f64> },
    Forward,
    Reverse,
}

impl Default for StrandBias {
    fn default() -> Self {
        StrandBias::None {
            forward_rate: NotNan::new(0.5).unwrap(), // TODO remove
        }
    }
}

impl Bias for StrandBias {
    fn prob_alt(&self, observation: &Observation<ReadPosition>) -> LogProb {
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
            }
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

    fn is_informative(&self, pileups: &[Vec<Observation<ReadPosition>>]) -> bool {
        // METHOD: if all reads come from the forward or reverse strand,
        // we cannot estimate a strand bias and None is the only informative one.
        !self.is_artifact() || Self::estimate_forward_rate(pileups).is_some()
    }

    fn learn_parameters(&mut self, pileups: &[Vec<Observation<ReadPosition>>]) {
        if let StrandBias::None {
            ref mut forward_rate,
        } = self
        {
            // METHOD: either we can estimate the forward rate, or non-artifact biases are discarded by
            // is_informative(). In that case, it is safe to just fall back to a forward rate of 0.5.
            *forward_rate =
                Self::estimate_forward_rate(pileups).unwrap_or(NotNan::new(0.5).unwrap());
        }
    }
}

impl StrandBias {
    fn estimate_forward_rate(pileups: &[Vec<Observation<ReadPosition>>]) -> Option<NotNan<f64>> {
        let strong_all = pileups
            .iter()
            .map(|pileup| pileup.iter().filter(|obs| obs.is_strong_ref_support() && obs.strand != Strand::Both))
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
        
        if strong_all > 2 {
            let forward_fraction = strong_forward as f64 / strong_all as f64;
            if  forward_fraction >= 0.4 && forward_fraction <= 0.6 {
                return Some(NotNan::new(0.5).unwrap())
            }
        } 
        None
    }
}
