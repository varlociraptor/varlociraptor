use bio::stats::probs::LogProb;
use bio_types::sequence::SequenceReadPairOrientation;

use crate::grammar::{VAFRange, VAFSpectrum};
use crate::utils::PROB_05;
use crate::variants::evidence::observation::{Observation, ReadPosition};
use crate::variants::model::bias::Bias;
use crate::variants::model::AlleleFreq;

use super::{BiasPriorHyperLikelihood, EqualDistributionAssumption};

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord, EnumIter, Hash)]
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
    type Assumption = EqualDistributionAssumption;

    fn prob(&self, observation: &Observation<ReadPosition>) -> LogProb {
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

    fn prob_any(&self, _observation: &Observation<ReadPosition>) -> LogProb {
        *PROB_05
    }

    fn is_artifact(&self) -> bool {
        *self != ReadOrientationBias::None
    }

    fn is_informative(&self, pileups: &[Vec<Observation<ReadPosition>>]) -> bool {
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
            .flatten()
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
        let n: usize = pileups.iter().map(|pileup| pileup.len()).sum();
        (n_uncertain as f64) < (n as f64 / 2.0)
    }
}

impl BiasPriorHyperLikelihood for ReadOrientationBias {
    fn hyper_likelihood(
        f1r2_allele_freq: AlleleFreq,
        observation: &Observation<ReadPosition>,
    ) -> LogProb {
        match observation.read_orientation {
            SequenceReadPairOrientation::F1R2 => LogProb(f1r2_allele_freq.ln()),
            SequenceReadPairOrientation::F2R1 => LogProb(f1r2_allele_freq.ln()).ln_one_minus_exp(),
            _ => LogProb::ln_one(),
        }
    }
}
