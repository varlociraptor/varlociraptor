use bio::stats::probs::LogProb;

use crate::utils::{PROB_025, PROB_033};
use crate::variants::evidence::observation::{IndelOperations, Observation, ReadPosition};
use crate::variants::model::bias::Bias;

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord, EnumIter, Hash)]
pub(crate) enum DivIndelBias {
    None,
    Some,
}

impl Default for DivIndelBias {
    fn default() -> Self {
        DivIndelBias::None
    }
}

impl Bias for DivIndelBias {
    fn prob(&self, observation: &Observation<ReadPosition, IndelOperations>) -> LogProb {
        match (self, observation.indel_operations) {
            (DivIndelBias::None, IndelOperations::Primary) => *PROB_033,
            (DivIndelBias::None, IndelOperations::Secondary) => *PROB_033,
            (DivIndelBias::None, IndelOperations::Other) => LogProb::ln_zero(),
            (DivIndelBias::None, IndelOperations::None) => *PROB_033,
            (DivIndelBias::Some, IndelOperations::Primary) => *PROB_025,
            (DivIndelBias::Some, IndelOperations::Secondary) => *PROB_025,
            (DivIndelBias::Some, IndelOperations::Other) => *PROB_025,
            (DivIndelBias::Some, IndelOperations::None) => *PROB_025,
        }
    }

    fn prob_any(&self, _observation: &Observation<ReadPosition, IndelOperations>) -> LogProb {
        *PROB_025
    }

    fn is_artifact(&self) -> bool {
        *self != DivIndelBias::None
    }

    fn is_informative(&self, pileups: &[Vec<Observation<ReadPosition, IndelOperations>>]) -> bool {
        if !self.is_artifact() {
            return true;
        }
        // METHOD: this bias is only relevant if there is at least one recorded indel operation (indel operations are only recorded for some variants).
        pileups.iter().any(|pileup| {
            pileup
                .iter()
                .any(|obs| obs.indel_operations != IndelOperations::None)
        })
    }
}
