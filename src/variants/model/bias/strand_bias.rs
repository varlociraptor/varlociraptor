use bio::stats::probs::LogProb;

use crate::utils::PROB_HALF;
use crate::variants::evidence::observation::{Observation, ReadPosition, Strand};
use crate::variants::model::bias::Bias;

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord, EnumIter)]
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
    fn prob(&self, observation: &Observation<ReadPosition>) -> LogProb {
        match (self, observation.strand) {
            (StrandBias::Forward, Strand::Forward) => LogProb::ln_one(),
            (StrandBias::Reverse, Strand::Forward) => LogProb::ln_zero(),
            (StrandBias::Forward, Strand::Reverse) => LogProb::ln_zero(),
            (StrandBias::Reverse, Strand::Reverse) => LogProb::ln_one(),
            (StrandBias::Forward, Strand::Both) => LogProb::ln_zero(),
            (StrandBias::Reverse, Strand::Both) => LogProb::ln_zero(),
            (StrandBias::None, Strand::Both) => observation.prob_double_overlap,
            (StrandBias::None, _) => *PROB_HALF + observation.prob_single_overlap,
            (_, Strand::None) => unreachable!(),
        }
    }

    fn prob_any(_observation: &Observation<ReadPosition>) -> LogProb {
        *PROB_HALF
    }

    fn is_artifact(&self) -> bool {
        *self != StrandBias::None
    }
}
