use bio::stats::probs::LogProb;

use crate::variants::evidence::observations::pileup::Pileup;
use crate::variants::evidence::observations::read_observation::ProcessedReadObservation;
use crate::variants::model::bias::Bias;

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord, EnumIter, Hash, Default)]
pub(crate) enum SoftclipBias {
    #[default]
    None,
    Some,
}

impl Bias for SoftclipBias {
    fn prob_alt(&self, observation: &ProcessedReadObservation) -> LogProb {
        match (self, observation.softclipped) {
            (SoftclipBias::Some, true) => LogProb::ln_one(),
            (SoftclipBias::Some, false) => LogProb::ln_zero(),
            (SoftclipBias::None, _) => LogProb::ln_one(),
        }
    }

    fn prob_any(&self, _observation: &ProcessedReadObservation) -> LogProb {
        LogProb::ln_one()
    }

    fn is_artifact(&self) -> bool {
        *self != SoftclipBias::None
    }

    fn is_informative(&self, pileups: &[Pileup]) -> bool {
        if !self.is_artifact() {
            return true;
        }
        // METHOD: this bias is only relevant if there is at least one softclip.
        pileups
            .iter()
            .any(|pileup| pileup.read_observations().iter().any(|obs| obs.softclipped))
    }
}
