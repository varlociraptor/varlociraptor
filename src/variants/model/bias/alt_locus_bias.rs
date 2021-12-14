use bio::stats::probs::LogProb;

use crate::utils::PROB_05;
use crate::variants::evidence::observation::{ProcessedObservation, AltLocus};
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
        // METHOD: either a reduced MAPQ or a read pointing to the major alt locus (or both)
        // are indicative for the variant to come from a different (distant) allele.
        // If all alt reads agree on this, we consider the bias to be present.
        match (self, observation.is_max_mapq, observation.alt_locus) {
            (AltLocusBias::None, _, _) => *PROB_05, // normal
            (AltLocusBias::Some, true, AltLocus::Some | AltLocus::None) => LogProb::ln_zero(), // no bias
            (AltLocusBias::Some, true, AltLocus::Major) => LogProb::ln_one(), // bias
            (AltLocusBias::Some, false, _) => LogProb::ln_one(), // bias
        }
    }

    fn prob_any(&self, observation: &ProcessedObservation) -> LogProb {
        *PROB_05
    }

    fn is_artifact(&self) -> bool {
        *self != AltLocusBias::None
    }
}
