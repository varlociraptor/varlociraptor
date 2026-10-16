use bio::stats::probs::LogProb;

use crate::variants::evidence::observations::read_observation::{
    FragmentEnd, ProcessedReadObservation,
};
use crate::variants::model::bias::Bias;

#[derive(Copy, Clone, Default, PartialOrd, PartialEq, Eq, Debug, Ord, EnumIter, Hash)]
pub(crate) enum FragmentEndBias {
    #[default]
    None,
    Some,
}

impl Bias for FragmentEndBias {
    fn prob_alt(&self, observation: &ProcessedReadObservation) -> LogProb {
        match (self, observation.fragment_end) {
            (FragmentEndBias::None, FragmentEnd::Major) => LogProb::ln_one(),
            (FragmentEndBias::None, FragmentEnd::Some) => {
                LogProb::ln_one()
            }
            (FragmentEndBias::Some, FragmentEnd::Major) => LogProb::ln_one(), // bias
            (FragmentEndBias::Some, FragmentEnd::Some) => LogProb::ln_zero(), // no bias
        }
    }

    fn prob_any(&self, _: &ProcessedReadObservation) -> LogProb {
        LogProb::ln_one()
    }

    fn is_artifact(&self) -> bool {
        *self != FragmentEndBias::None
    }
}
