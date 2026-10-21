// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use bio::stats::{LogProb};

use crate::utils::NUMERICAL_EPSILON;
use crate::variants::evidence::observations::read_observation::ProcessedReadObservation;
use crate::variants::model::bias::Artifacts;
use crate::variants::model::AlleleFreq;

pub(crate) mod contaminated_sample;
pub(crate) mod single_sample;

#[derive(PartialEq, Eq, Debug, Clone, Hash)]
pub(crate) struct Event {
    pub(crate) allele_freq: AlleleFreq,
    pub(crate) artifacts: Artifacts,
    pub(crate) is_discrete: bool,
}

impl Event {
    pub(crate) fn absent() -> Self {
        Event {
            allele_freq: AlleleFreq(0.0),
            artifacts: Artifacts::none(),
            is_discrete: true,
        }
    }

    pub(crate) fn is_artifact(&self) -> bool {
        self.artifacts.is_artifact()
    }

    pub(crate) fn is_absent(&self) -> bool {
        *self.allele_freq == 0.0
    }
}

fn prob_sample_alt(observation: &ProcessedReadObservation, allele_freq: LogProb) -> LogProb {
    if allele_freq != LogProb::ln_one() {
        // The effective sample probability for the alt allele is the allele frequency times
        // the probability to obtain a feasible fragment (prob_sample_alt).
        (allele_freq + observation.prob_sample_alt).cap_numerical_overshoot(NUMERICAL_EPSILON)
    } else {
        // If allele frequency is 1.0, sampling bias does have no effect because all reads
        // should come from the alt allele.
        allele_freq
    }
}

/// Calculate likelihood of allele freq given observation in a single sample assuming that the
/// underlying fragment/read is mapped correctly.
fn likelihood_mapping(
    allele_freq: LogProb,
    biases: &Artifacts,
    observation: &ProcessedReadObservation,
) -> LogProb {
    // Step 1: calculate probability to sample from alt allele
    let prob_sample_alt = prob_sample_alt(observation, allele_freq);
    let prob_sample_ref = prob_sample_alt.ln_one_minus_exp();

    let prob_bias_alt = biases.prob_alt(observation);
    let prob_bias_ref = biases.prob_ref(observation);

    // Step 2: read comes from case sample and is correctly mapped
    let prob = LogProb::ln_sum_exp(&[
        // alt allele
        prob_sample_alt + prob_bias_alt + observation.prob_alt(),
        // ref allele
        prob_sample_ref + observation.prob_ref() + prob_bias_ref,
    ][..]);
    assert!(!prob.is_nan());

    prob
}