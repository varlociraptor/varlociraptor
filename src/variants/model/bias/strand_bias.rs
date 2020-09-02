use bio::stats::probs::LogProb;

use crate::utils::PROB05;
use crate::variants::evidence::observation::Observation;
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
        let obs_strand = (observation.forward_strand, observation.reverse_strand);

        match (self, obs_strand) {
            (StrandBias::Forward, (true, false)) => LogProb::ln_one(),
            (StrandBias::Reverse, (true, false)) => LogProb::ln_zero(),
            (StrandBias::Forward, (false, true)) => LogProb::ln_zero(),
            (StrandBias::Reverse, (false, true)) => LogProb::ln_one(),
            (StrandBias::Forward, (true, true)) => LogProb::ln_zero(),
            (StrandBias::Reverse, (true, true)) => LogProb::ln_zero(),
            (StrandBias::None, _) => {
                if observation.forward_strand != observation.reverse_strand {
                    *PROB05 + observation.prob_single_overlap
                } else {
                    observation.prob_double_overlap
                }
            }
            (_, (false, false)) => unreachable!(),
        }
    }

    fn prob_any(&self) -> LogProb {
        *PROB05
    }

    fn is_artifact(&self) -> bool {
        *self != StrandBias::None
    }
}
