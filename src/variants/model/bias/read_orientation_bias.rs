use bio::stats::probs::LogProb;
use bio_types::sequence::SequenceReadPairOrientation;

use crate::utils::PROB_HALF;
use crate::variants::evidence::observation::{Observation, ReadPosition};
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
    fn prob(&self, observation: &Observation<ReadPosition>) -> LogProb {
        match (self, observation.read_orientation) {
            (ReadOrientationBias::None, SequenceReadPairOrientation::F1R2) => *PROB_HALF, // normal
            (ReadOrientationBias::None, SequenceReadPairOrientation::F2R1) => *PROB_HALF, // normal
            (ReadOrientationBias::F1R2, SequenceReadPairOrientation::F1R2) => LogProb::ln_one(), // bias
            (ReadOrientationBias::F2R1, SequenceReadPairOrientation::F2R1) => LogProb::ln_one(), // bias
            (ReadOrientationBias::F1R2, SequenceReadPairOrientation::F2R1) => LogProb::ln_zero(), // no bias
            (ReadOrientationBias::F2R1, SequenceReadPairOrientation::F1R2) => LogProb::ln_zero(), // no bias
            _ => *PROB_HALF, // For None and nonstandard orientations, the true one can be either F1R2 or F2R1, hence 0.5.
        }
    }

    fn prob_any(&self, _observation: &Observation<ReadPosition>) -> LogProb {
        *PROB_HALF
    }

    fn is_artifact(&self) -> bool {
        *self != ReadOrientationBias::None
    }

    fn is_possible(&self, pileups: &[Vec<Observation<ReadPosition>>]) -> bool {
        pileups.iter().any(|pileup| {
            pileup.iter().any(|observation| {
                if let ReadOrientationBias::None = *self {
                    true
                } else {
                    // METHOD: we only statistically consider a bias if there is at least one ALT observation
                    // that votes for it without uncertainty. This way, we omit bias estimation in totally
                    // uncertain cases, where the bias is actually not observed.
                    // This behavior is important, because otherwise pathological cases like
                    // 98Ns**^68Ns***18Vs**^4Vs***2Ns-*^1Np+<*1Np+>*1Ns-**, would weakly vote for a bias,
                    // although 90% of the reads do not provide orientation info and the only that do are
                    // reference reads.
                    self.prob(observation) == LogProb::ln_one()
                        && *observation.bayes_factor_alt() > 1.0
                }
            })
        })
    }
}
