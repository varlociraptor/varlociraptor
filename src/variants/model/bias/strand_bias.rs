use bio::stats::hmm::forward;
use bio::stats::probs::LogProb;

use crate::grammar::{VAFRange, VAFSpectrum};
use crate::utils::PROB_05;
use crate::variants::evidence::observation::{Observation, ReadPosition, Strand};
use crate::variants::model::bias::Bias;
use crate::variants::model::AlleleFreq;

use super::{BiasPriorHyperLikelihood, EqualDistributionAssumption};

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
    type Assumption = EqualDistributionAssumption;

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
}

impl BiasPriorHyperLikelihood for StrandBias {
    fn hyper_likelihood(
        forward_allele_freq: AlleleFreq,
        observation: &Observation<ReadPosition>,
    ) -> LogProb {
        match observation.strand {
            Strand::Forward => LogProb(forward_allele_freq.ln()),
            Strand::Reverse => LogProb(forward_allele_freq.ln()).ln_one_minus_exp(),
            Strand::Both => LogProb::ln_one(),
            Strand::None => LogProb::ln_one(),
        }
    }
}
