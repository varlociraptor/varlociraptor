use bio::stats::probs::LogProb;

use crate::variants::evidence::observation::{IndelOperations, Observation, ReadPosition};
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
    fn prob(&self, observation: &Observation<ReadPosition>) -> LogProb {
        match (self, observation.read_position) {
            (ReadPositionBias::None, _) => observation.prob_hit_base, // normal
            (ReadPositionBias::Some, ReadPosition::Major) => LogProb::ln_one(), // bias
            (ReadPositionBias::Some, ReadPosition::Some) => LogProb::ln_zero(), // no bias
        }
    }

    fn prob_any(&self, observation: &Observation<ReadPosition>) -> LogProb {
        observation.prob_hit_base
    }

    fn is_artifact(&self) -> bool {
        *self != ReadPositionBias::None
    }
}
