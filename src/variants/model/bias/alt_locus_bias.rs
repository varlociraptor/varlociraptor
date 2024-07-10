use bio::stats::probs::LogProb;

use crate::utils::PROB_05;
use crate::variants::evidence::observations::pileup::Pileup;
use crate::variants::evidence::observations::read_observation::{
    AltLocus, ProcessedReadObservation,
};
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
    fn prob_alt(&self, observation: &ProcessedReadObservation) -> LogProb {
        // METHOD: a read pointing to the major alt locus
        // are indicative for the variant to come from a different (distant) allele.
        // If all alt reads agree on this, we consider the bias to be present.
        match (self, observation.is_max_mapq, observation.alt_locus) {
            (AltLocusBias::None, _, _) => *PROB_05, // normal
            (AltLocusBias::Some, true, AltLocus::Some | AltLocus::None) => LogProb::ln_zero(), // no bias
            (AltLocusBias::Some, true, AltLocus::Major) => LogProb::ln_one(), // bias
            (AltLocusBias::Some, false, AltLocus::Major) => LogProb::ln_one(), // bias
            (AltLocusBias::Some, false, AltLocus::Some | AltLocus::None) => LogProb::ln_zero(), // no bias
        }
    }

    fn prob_ref(&self, observation: &ProcessedReadObservation) -> LogProb {
        // METHOD: ref reads should not point to the alt locus. The reason is that in that case,
        // the homology does not appear to be variant specific, and hence the normal MAPQs
        // should be able to capture it.
        match (self, observation.is_max_mapq, observation.alt_locus) {
            (AltLocusBias::None, _, _) => *PROB_05, // normal
            (AltLocusBias::Some, true, AltLocus::Some | AltLocus::None) => LogProb::ln_one(), // no bias
            (AltLocusBias::Some, true, AltLocus::Major) => LogProb::ln_zero(), // bias
            (AltLocusBias::Some, false, AltLocus::Major) => LogProb::ln_zero(), // bias
            (AltLocusBias::Some, false, AltLocus::Some | AltLocus::None) => LogProb::ln_one(), // no bias
        }
    }

    fn prob_any(&self, _observation: &ProcessedReadObservation) -> LogProb {
        *PROB_05
    }

    fn is_artifact(&self) -> bool {
        *self != AltLocusBias::None
    }

    fn is_informative(&self, pileups: &[Pileup]) -> bool {
        // METHOD: we consider this bias if more than 10% of the reads does not have the maximum MAPQ.
        if !self.is_artifact() {
            return true;
        }

        let n: usize = pileups
            .iter()
            .map(|pileup| {
                pileup
                    .read_observations()
                    .iter()
                    .filter(|obs| obs.is_strong_alt_support())
                    .count()
            })
            .sum();
        let non_max_mapq: usize = pileups
            .iter()
            .map(|pileup| {
                pileup
                    .read_observations()
                    .iter()
                    .filter(|obs| !obs.is_max_mapq && obs.is_strong_alt_support())
                    .count()
            })
            .sum();
        n > 0 && non_max_mapq as f64 > (n as f64 * 0.01) && (n - non_max_mapq) < 10
    }
}
