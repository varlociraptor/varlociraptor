use bio::stats::probs::LogProb;

use crate::utils::PROB_05;
use crate::variants::evidence::observation::{AltLocus, ProcessedObservation};
use crate::variants::model::bias::Bias;

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord, EnumIter, Hash)]
pub(crate) enum AltLocusBias {
    None,
    Some,
}

impl Default for AltLocusBias {
    fn default() -> Self {
        AltLocusBias::None
    }
}

impl Bias for AltLocusBias {
    fn prob_alt(&self, observation: &ProcessedObservation) -> LogProb {
        // METHOD: a read pointing to the major alt locus
        // are indicative for the variant to come from a different (distant) allele.
        // If all alt reads agree on this, we consider the bias to be present.
        match (self, observation.is_max_mapq, observation.alt_locus) {
            (AltLocusBias::None, _, _) => *PROB_05, // normal
            (AltLocusBias::Some, true, AltLocus::Some | AltLocus::None) => LogProb::ln_zero(), // no bias
            (AltLocusBias::Some, true, AltLocus::Major) => LogProb::ln_one(), // bias
            (AltLocusBias::Some, false, AltLocus::Major) => LogProb::ln_one(),              // bias
            (AltLocusBias::Some, false, AltLocus::Some | AltLocus::None) => LogProb::ln_zero(),              // no bias
        }
    }

    fn prob_any(&self, observation: &ProcessedObservation) -> LogProb {
        *PROB_05
    }

    fn is_artifact(&self) -> bool {
        *self != AltLocusBias::None
    }

    fn is_informative(&self, pileups: &[Vec<ProcessedObservation>]) -> bool {
        // METHOD: we consider this bias if there is at least one non maximum MAPQ read.
        !self.is_artifact() || pileups.iter().any(|pileup| {
            pileup.iter().any(|obs| !obs.is_max_mapq)
        })
    }
}
