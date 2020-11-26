use bio::stats::probs::LogProb;

use crate::utils::PROB_05;
use crate::variants::evidence::observation::{Observation, ReadOrientation};
use crate::variants::model::bias::Bias;

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord, EnumIter)]
pub(crate) enum ReadOrientationBias {
    None,
    F1R2,
    F2R1,
}

impl Default for ReadOrientationBias {
    fn default() -> Self {
        ReadOrientationBias::None
    }
}

impl Bias for ReadOrientationBias {
    fn prob(&self, observation: &Observation) -> LogProb {
        match (self, observation.read_orientation) {
            (ReadOrientationBias::None, ReadOrientation::F1R2) => *PROB_05, // normal
            (ReadOrientationBias::None, ReadOrientation::F2R1) => *PROB_05, // normal
            (ReadOrientationBias::F1R2, ReadOrientation::F1R2) => LogProb::ln_one(), // bias
            (ReadOrientationBias::F2R1, ReadOrientation::F2R1) => LogProb::ln_one(), // bias
            (ReadOrientationBias::F1R2, ReadOrientation::F2R1) => LogProb::ln_zero(), // no bias
            (ReadOrientationBias::F2R1, ReadOrientation::F1R2) => LogProb::ln_zero(), // no bias
            (_, ReadOrientation::None) => LogProb::ln_one(), // If we have no information, everything is equally possible.
            (ReadOrientationBias::None, _) => LogProb::ln_one(), // no bias and "other" observations
            _ => LogProb::ln_zero(), // other observations make any bias impossible
        }
    }

    fn prob_any(&self) -> LogProb {
        *PROB_05
    }

    fn is_artifact(&self) -> bool {
        *self != ReadOrientationBias::None
    }
}
