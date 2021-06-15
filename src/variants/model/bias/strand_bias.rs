use bio::stats::probs::LogProb;

use crate::utils::PROB_05;
use crate::variants::evidence::observation::{Observation, ReadPosition, Strand};
use crate::variants::model::bias::Bias;

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord, EnumIter, Hash)]
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
            (StrandBias::None, _) => *PROB_05 + observation.prob_single_overlap,
            (_, Strand::None) => unreachable!(),
        }
    }

    fn prob_any(&self, _observation: &Observation<ReadPosition>) -> LogProb {
        *PROB_05
    }

    fn is_artifact(&self) -> bool {
        *self != StrandBias::None
    }

    fn is_informative(&self, pileups: &[Vec<Observation<ReadPosition>>]) -> bool {
        if let StrandBias::None = *self {
            return true;
        }
        // METHOD: strand bias is only informative if we observe both strands in the pileup.
        // Otherwise, we likely have a situation caused by targeted sequencing where we should
        // not infer the strand bias (e.g. only covered by a single amplicon or at the beginning
        // of a target region).
        pileups
            .iter()
            .flatten()
            .any(|observation| observation.strand == Strand::Forward)
            && pileups
                .iter()
                .flatten()
                .any(|observation| observation.strand == Strand::Reverse)
    }
}
