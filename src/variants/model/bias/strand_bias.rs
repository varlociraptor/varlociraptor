use bio::stats::probs::LogProb;

use crate::utils::PROB05;
use crate::variants::evidence::observation::{Observation, Strand};
use crate::variants::model::bias::Bias;

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord)]
pub(crate) enum StrandBias {
    None,
    Forward,
    Reverse,
}

impl Default for StrandBias {
    fn default() -> Self {
        StrandBias::None
    }
}

impl Bias for StrandBias {
    fn prob(&self, observation: &Observation) -> LogProb {
        match (self, observation.strand) {
            (StrandBias::Forward, Strand::Forward) => LogProb::ln_one(),
            (StrandBias::Reverse, Strand::Forward) => LogProb::ln_zero(),
            (StrandBias::Forward, Strand::Reverse) => LogProb::ln_zero(),
            (StrandBias::Reverse, Strand::Reverse) => LogProb::ln_one(),
            (StrandBias::Forward, Strand::Both) => LogProb::ln_zero(),
            (StrandBias::Reverse, Strand::Both) => LogProb::ln_zero(),
            (StrandBias::None, Strand::Both) => observation.prob_double_overlap,
            (StrandBias::None, _) => *PROB05 + observation.prob_single_overlap,
            (_, Strand::None) => unreachable!(),
        }
    }

    fn prob_any(&self) -> LogProb {
        *PROB05
    }

    fn is_artifact(&self) -> bool {
        *self != StrandBias::None
    }
}
