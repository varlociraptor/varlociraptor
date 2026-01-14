use bio::stats::probs::LogProb;

use crate::utils::PROB_05;
use crate::variants::evidence::observations::pileup::Pileup;
use crate::variants::evidence::observations::read_observation::{
    AltLocus, ProcessedReadObservation, ReadObservation, ReadPosition,
};
use crate::variants::model::bias::Bias;

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord, EnumIter, Hash, Default)]
pub(crate) enum AltLocusBias {
    #[default]
    None,
    Some {
        has_alt_loci: bool,
    },
}

impl AltLocusBias {
    fn get_counts<F>(pileups: &[Pileup], filterfunc: F) -> (usize, usize)
    where
        F: Fn(&ReadObservation<ReadPosition, AltLocus>) -> bool,
    {
        let n: usize = pileups
            .iter()
            .map(|pileup| {
                pileup
                    .read_observations()
                    .iter()
                    .filter(|obs| filterfunc(*obs))
                    .count()
            })
            .sum();
        let non_max_mapq: usize = pileups
            .iter()
            .map(|pileup| {
                pileup
                    .read_observations()
                    .iter()
                    .filter(|obs| !obs.is_max_mapq && filterfunc(*obs))
                    .count()
            })
            .sum();
        (n, non_max_mapq)
    }

    fn has_alt_loci(pileups: &[Pileup]) -> bool {
        pileups
            .iter()
            .map(|pileup| {
                pileup
                    .read_observations()
                    .iter()
                    .filter(|obs| obs.alt_locus != AltLocus::None)
                    .count()
            })
            .sum::<usize>()
            > 0
    }
}

impl Bias for AltLocusBias {
    fn prob_alt(&self, observation: &ProcessedReadObservation) -> LogProb {
        // METHOD: a read pointing to the major alt locus
        // are indicative for the variant to come from a different (distant) allele.
        // If all alt reads agree on this, we consider the bias to be present.
        if let AltLocusBias::Some { has_alt_loci } = self {
            if *has_alt_loci {
                match observation.alt_locus {
                    AltLocus::Some | AltLocus::None => LogProb::ln_zero(),
                    AltLocus::Major => LogProb::ln_one(),
                }
            } else {
                if observation.is_max_mapq {
                    LogProb::ln_zero()
                } else {
                    LogProb::ln_one()
                }
            }
        } else {
            // AltLocus::None
            *PROB_05
        }
    }

    fn prob_ref(&self, observation: &ProcessedReadObservation) -> LogProb {
        // METHOD: ref reads should not point to the alt locus. The reason is that in that case,
        // the homology does not appear to be variant specific, and hence the normal MAPQs
        // should be able to capture it.
        // We thus invert the logic of prob_alt here.
        if let AltLocusBias::Some { has_alt_loci } = self {
            if *has_alt_loci {
                match observation.alt_locus {
                    AltLocus::Some | AltLocus::None => LogProb::ln_one(), // no bias
                    AltLocus::Major => LogProb::ln_zero(),
                }
            } else {
                // METHOD: above logic only applies for a concrete ALT locus.
                // In case we can only rely on MAPQ we should be a bit more conservative,
                // allowing basically anything to occur.
                *PROB_05
            }
        } else {
            // AltLocus::None
            *PROB_05
        }
    }

    fn prob_any(&self, _observation: &ProcessedReadObservation) -> LogProb {
        *PROB_05
    }

    fn is_artifact(&self) -> bool {
        *self != AltLocusBias::None
    }

    fn learn_parameters(&mut self, pileups: &[Pileup]) {
        if let AltLocusBias::Some {
            ref mut has_alt_loci,
        } = self
        {
            *has_alt_loci = Self::has_alt_loci(pileups);
        }
    }

    fn is_informative(&self, pileups: &[Pileup]) -> bool {
        // METHOD: we consider this bias if more than 10% of the ALT reads does not have the maximum MAPQ
        // and either there is any AltLocus::Some/Major observation or if the REF does have at least 10% max MAPQ observations.
        if !self.is_artifact() {
            return true;
        }

        let (n_alt, non_max_mapq_alt) =
            Self::get_counts(pileups, |obs| obs.is_strong_alt_support());
        let (n_ref, non_max_mapq_ref) =
            Self::get_counts(pileups, |obs| obs.is_strong_ref_support());

        let enough_alt_locus_obs = n_alt > 0
            && non_max_mapq_alt as f64 > (n_alt as f64 * 0.01)
            && (n_alt - non_max_mapq_alt) < 10;
        let enough_max_mapq_in_ref =
            n_ref > 0 && ((non_max_mapq_ref as f64) < (n_ref as f64 * 0.9));
        let has_alt_loci = Self::has_alt_loci(pileups);

        enough_alt_locus_obs && (has_alt_loci || enough_max_mapq_in_ref)
    }
}
