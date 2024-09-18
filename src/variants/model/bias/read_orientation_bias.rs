use bio::stats::probs::LogProb;
use bio_types::sequence::SequenceReadPairOrientation;

use crate::utils::PROB_05;
use crate::variants::evidence::observations::pileup::Pileup;
use crate::variants::evidence::observations::read_observation::ProcessedReadObservation;
use crate::variants::model::bias::Bias;

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord, EnumIter, Hash, Default)]
pub(crate) enum ReadOrientationBias {
    #[default]
    None,
    F1R2,
    F2R1,
}

impl Bias for ReadOrientationBias {
    fn prob_alt(&self, observation: &ProcessedReadObservation) -> LogProb {
        match (self, observation.read_orientation) {
            (ReadOrientationBias::None, SequenceReadPairOrientation::F1R2) => *PROB_05, // normal
            (ReadOrientationBias::None, SequenceReadPairOrientation::F2R1) => *PROB_05, // normal
            (ReadOrientationBias::F1R2, SequenceReadPairOrientation::F1R2) => LogProb::ln_one(), // bias
            (ReadOrientationBias::F2R1, SequenceReadPairOrientation::F2R1) => LogProb::ln_one(), // bias
            (ReadOrientationBias::F1R2, SequenceReadPairOrientation::F2R1) => LogProb::ln_zero(), // no bias
            (ReadOrientationBias::F2R1, SequenceReadPairOrientation::F1R2) => LogProb::ln_zero(), // no bias
            _ => *PROB_05, // For None and nonstandard orientations, the true one can be either F1R2 or F2R1, hence 0.5.
        }
    }

    fn prob_any(&self, _observation: &ProcessedReadObservation) -> LogProb {
        *PROB_05
    }

    fn is_artifact(&self) -> bool {
        *self != ReadOrientationBias::None
    }

    fn is_informative(&self, pileups: &[Pileup]) -> bool {
        if let ReadOrientationBias::None = *self {
            return true;
        }
        // METHOD: read orientation bias is only informative if the majority of the reads
        // provide orientation information. This way, we omit bias estimation in totally
        // uncertain cases, where the bias is actually not observed.
        // This behavior is important, because otherwise pathological cases like
        // 98Ns**^68Ns***18Vs**^4Vs***2Ns-*^1Np+<*1Np+>*1Ns-**, would weakly vote for a bias,
        // although almost all reads do not provide orientation info and the only that do are
        // reference reads.
        let n_uncertain: usize = pileups
            .iter()
            .flat_map(|pileup| pileup.read_observations().iter())
            .map(|observation| {
                if !(observation.read_orientation == SequenceReadPairOrientation::F1R2
                    || observation.read_orientation == SequenceReadPairOrientation::F2R1)
                {
                    1
                } else {
                    0
                }
            })
            .sum();
        let n: usize = pileups
            .iter()
            .map(|pileup| pileup.read_observations().len())
            .sum();
        let enough_information = (n_uncertain as f64) < (n as f64 / 2.0);

        // METHOD: read orientation bias needs a uniform distribution of F1R2 and F2R1 among the
        // reference reads. Otherwise, we cannot reliably detect whether there is something odd
        // in the alt reads.
        let strong_ref_total_count = pileups
            .iter()
            .flat_map(|pileup| pileup.read_observations().iter())
            .filter(|observation| {
                observation.is_strong_ref_support()
                    && (observation.read_orientation == SequenceReadPairOrientation::F1R2
                        || observation.read_orientation == SequenceReadPairOrientation::F2R1)
            })
            .count();
        let strong_ref_f1r2 = pileups
            .iter()
            .flat_map(|pileup| pileup.read_observations().iter())
            .filter(|observation| {
                observation.is_strong_ref_support()
                    && (observation.read_orientation == SequenceReadPairOrientation::F1R2)
            })
            .count();
        let uniform_distribution = if strong_ref_total_count > 2 {
            let fraction = strong_ref_f1r2 as f64 / strong_ref_total_count as f64;
            // TODO use strong_ref_total_count and binomial to calculate a confidence interval
            (0.3..=0.7).contains(&fraction)
        } else {
            false
        };

        dbg!((
            enough_information,
            uniform_distribution,
            strong_ref_total_count,
            strong_ref_f1r2
        ));

        enough_information && uniform_distribution
    }
}
