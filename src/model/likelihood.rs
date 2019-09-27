// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::BTreeMap;

use bio::stats::{bayesian::model::Likelihood, LogProb, Prob};

use crate::model::evidence::Observation;
use crate::model::sample::Pileup;
use crate::model::{AlleleFreq, StrandBias};
use crate::utils::NUMERICAL_EPSILON;

pub type ContaminatedSampleCache = BTreeMap<ContaminatedSampleEvent, LogProb>;
pub type SingleSampleCache = BTreeMap<Event, LogProb>;

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone)]
pub struct Event {
    pub allele_freq: AlleleFreq,
    pub strand_bias: StrandBias,
}

fn prob_sample_alt(observation: &Observation, allele_freq: LogProb) -> LogProb {
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

pub trait ContaminatedSamplePairView<T> {
    fn primary(&self) -> &T;
    fn secondary(&self) -> &T;
}

impl<T> ContaminatedSamplePairView<T> for Vec<T> {
    fn primary(&self) -> &T {
        &self[0]
    }

    fn secondary(&self) -> &T {
        &self[1]
    }
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone)]
pub struct ContaminatedSampleEvent {
    pub primary: Event,
    pub secondary: Event,
}

/// Variant calling model, taking purity and allele frequencies into account.
#[derive(Clone, Copy, Debug)]
pub struct ContaminatedSampleLikelihoodModel {
    /// Purity of the case sample.
    purity: LogProb,
    impurity: LogProb,
}

impl Default for ContaminatedSampleLikelihoodModel {
    fn default() -> Self {
        ContaminatedSampleLikelihoodModel::new(1.0)
    }
}

impl ContaminatedSampleLikelihoodModel {
    /// Create new model.
    pub fn new(purity: f64) -> Self {
        assert!(purity > 0.0 && purity <= 1.0);
        let purity = LogProb(purity.ln());
        ContaminatedSampleLikelihoodModel {
            purity: purity,
            impurity: purity.ln_one_minus_exp(),
        }
    }

    fn likelihood_observation(
        &self,
        allele_freq_primary: LogProb,
        allele_freq_secondary: LogProb,
        strand_bias_primary: StrandBias,
        strand_bias_secondary: StrandBias,
        observation: &Observation,
    ) -> LogProb {
        // Step 1: likelihoods for the mapping case.
        // Case 1: read comes from primary sample and is correctly mapped
        let prob_primary =
            self.purity + likelihood_mapping(allele_freq_primary, strand_bias_primary, observation);
        // Case 2: read comes from secondary sample and is correctly mapped
        let prob_secondary = self.impurity
            + likelihood_mapping(allele_freq_secondary, strand_bias_secondary, observation);

        // Step 4: total probability
        // Important note: we need to multiply a probability for a hypothetical missed allele
        // in the mismapping case. Otherwise, it can happen that mismapping dominates subtle
        // differences in the likelihood for alt and ref allele with low probabilities and very
        // low allele frequencies, such that we loose sensitivity for those.
        let total = (observation.prob_mapping + prob_secondary.ln_add_exp(prob_primary))
            .ln_add_exp(observation.prob_mismapping + observation.prob_missed_allele);
        assert!(!total.is_nan());
        total
    }
}

impl Likelihood<ContaminatedSampleCache> for ContaminatedSampleLikelihoodModel {
    type Event = ContaminatedSampleEvent;
    type Data = Pileup;

    fn compute(
        &self,
        events: &Self::Event,
        pileup: &Self::Data,
        cache: &mut ContaminatedSampleCache,
    ) -> LogProb {
        if cache.contains_key(events) {
            *cache.get(events).unwrap()
        } else {
            let ln_af_primary = LogProb(events.primary.allele_freq.ln());
            let ln_af_secondary = LogProb(events.secondary.allele_freq.ln());

            // calculate product of per-oservation likelihoods in log space
            let likelihood = pileup.iter().fold(LogProb::ln_one(), |prob, obs| {
                let lh = self.likelihood_observation(
                    ln_af_primary,
                    ln_af_secondary,
                    events.primary.strand_bias,
                    events.secondary.strand_bias,
                    obs,
                );
                prob + lh
            });

            assert!(!likelihood.is_nan());
            cache.insert(events.clone(), likelihood);

            likelihood
        }
    }
}

lazy_static! {
    static ref PROB_STRAND: LogProb = LogProb(0.5_f64.ln());
}

/// Likelihood model for single sample.
#[derive(Clone, Copy, Debug, Default)]
pub struct SampleLikelihoodModel {}

impl SampleLikelihoodModel {
    /// Create new model.
    pub fn new() -> Self {
        SampleLikelihoodModel {}
    }

    /// Likelihood to observe a read given allele frequency for a single sample.
    fn likelihood_observation(
        &self,
        allele_freq: LogProb,
        strand_bias: StrandBias,
        observation: &Observation,
    ) -> LogProb {
        // Step 1: likelihood for the mapping case.
        let prob = likelihood_mapping(allele_freq, strand_bias, observation);

        // Step 2: total probability
        // Important note: we need to multiply a probability for a hypothetical missed allele
        // in the mismapping case. Otherwise, it can happen that mismapping dominates subtle
        // differences in the likelihood for alt and ref allele with low probabilities and very
        // low allele frequencies, such that we loose sensitivity for those.
        let total = (observation.prob_mapping + prob)
            .ln_add_exp(observation.prob_mismapping + observation.prob_missed_allele);
        assert!(!total.is_nan());
        total
    }
}

/// Calculate likelihood of allele freq given observation in a single sample assuming that the
/// underlying fragment/read is mapped correctly.
fn likelihood_mapping(
    allele_freq: LogProb,
    strand_bias: StrandBias,
    observation: &Observation,
) -> LogProb {
    // Step 1: calculate probability to sample from alt allele
    let prob_sample_alt = prob_sample_alt(observation, allele_freq);
    let prob_sample_ref = prob_sample_alt.ln_one_minus_exp();

    let obs_strand = (observation.forward_strand, observation.reverse_strand);

    let prob_strand = match (strand_bias, obs_strand) {
        (StrandBias::Forward, (true, false)) => LogProb::ln_one(),
        (StrandBias::Reverse, (true, false)) => LogProb::ln_zero(),
        (StrandBias::Forward, (false, true)) => LogProb::ln_zero(),
        (StrandBias::Reverse, (false, true)) => LogProb::ln_one(),
        (StrandBias::Forward, (true, true)) => LogProb::ln_zero(),
        (StrandBias::Reverse, (true, true)) => LogProb::ln_zero(),
        (StrandBias::None, _) => {
            if observation.forward_strand != observation.reverse_strand {
                LogProb::from(Prob(0.5)) + observation.prob_single_overlap
            } else {
                observation.prob_double_overlap
            }
        }
        (_, (false, false)) => unreachable!(),
    };

    // Step 2: read comes from case sample and is correctly mapped
    let prob = LogProb::ln_sum_exp(&[
        // alt allele
        prob_sample_alt + prob_strand + observation.prob_alt,
        // ref allele (we don't care about the strand)
        prob_sample_ref + observation.prob_ref,
    ]);
    assert!(!prob.is_nan());

    prob
}

impl Likelihood<SingleSampleCache> for SampleLikelihoodModel {
    type Event = Event;
    type Data = Pileup;

    /// Likelihood to observe a pileup given allele frequencies for case and control.
    fn compute(&self, event: &Event, pileup: &Pileup, cache: &mut SingleSampleCache) -> LogProb {
        if cache.contains_key(event) {
            *cache.get(event).unwrap()
        } else {
            let ln_af = LogProb(event.allele_freq.ln());

            // calculate product of per-read likelihoods in log space
            let likelihood = pileup.iter().fold(LogProb::ln_one(), |prob, obs| {
                let lh = self.likelihood_observation(ln_af, event.strand_bias, obs);
                prob + lh
            });

            assert!(!likelihood.is_nan());

            cache.insert(event.clone(), likelihood);

            likelihood
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::likelihood;
    use crate::model::tests::observation;
    use crate::model::StrandBias;
    use bio::stats::LogProb;
    use itertools_num::linspace;

    fn event(allele_freq: f64) -> Event {
        Event {
            allele_freq: AlleleFreq(allele_freq),
            strand_bias: StrandBias::None,
        }
    }

    #[test]
    fn test_likelihood_observation_absent_single() {
        let observation = observation(LogProb::ln_one(), LogProb::ln_zero(), LogProb::ln_one());

        let model = SampleLikelihoodModel::new();

        let lh = model.likelihood_observation(
            LogProb(AlleleFreq(0.0).ln()),
            StrandBias::None,
            &observation,
        );
        assert_relative_eq!(*lh, *LogProb::ln_one());
    }

    #[test]
    fn test_likelihood_observation_absent() {
        let model = ContaminatedSampleLikelihoodModel::new(1.0);
        let observation = observation(LogProb::ln_one(), LogProb::ln_zero(), LogProb::ln_one());

        let lh = model.likelihood_observation(
            LogProb(AlleleFreq(0.0).ln()),
            LogProb(AlleleFreq(0.0).ln()),
            StrandBias::None,
            StrandBias::None,
            &observation,
        );
        assert_relative_eq!(*lh, *LogProb::ln_one());
    }

    #[test]
    fn test_likelihood_pileup_absent() {
        let model = ContaminatedSampleLikelihoodModel::new(1.0);
        let mut observations = Vec::new();
        for _ in 0..10 {
            observations.push(observation(
                LogProb::ln_one(),
                LogProb::ln_zero(),
                LogProb::ln_one(),
            ));
        }
        let mut cache = likelihood::ContaminatedSampleCache::default();

        let lh = model.compute(
            &ContaminatedSampleEvent {
                primary: event(0.0),
                secondary: event(0.0),
            },
            &observations,
            &mut cache,
        );
        assert_relative_eq!(*lh, *LogProb::ln_one());
    }

    #[test]
    fn test_likelihood_pileup_absent_single() {
        let model = SampleLikelihoodModel::new();
        let mut observations = Vec::new();
        for _ in 0..10 {
            observations.push(observation(
                LogProb::ln_one(),
                LogProb::ln_zero(),
                LogProb::ln_one(),
            ));
        }
        let mut cache = likelihood::SingleSampleCache::default();
        let evt = event(0.0);
        let lh = model.compute(&evt, &observations, &mut cache);
        assert_relative_eq!(*lh, *LogProb::ln_one());
        assert!(cache.contains_key(&evt))
    }

    #[test]
    fn test_likelihood_observation_case_control() {
        let model = ContaminatedSampleLikelihoodModel::new(1.0);
        let observation = observation(LogProb::ln_one(), LogProb::ln_one(), LogProb::ln_zero());

        let lh = model.likelihood_observation(
            LogProb(AlleleFreq(1.0).ln()),
            LogProb(AlleleFreq(0.0).ln()),
            StrandBias::None,
            StrandBias::None,
            &observation,
        );
        assert_relative_eq!(*lh, *LogProb::ln_one());

        let lh = model.likelihood_observation(
            LogProb(AlleleFreq(0.0).ln()),
            LogProb(AlleleFreq(0.0).ln()),
            StrandBias::None,
            StrandBias::None,
            &observation,
        );
        assert_relative_eq!(*lh, *LogProb::ln_zero());

        let lh = model.likelihood_observation(
            LogProb(AlleleFreq(0.5).ln()),
            LogProb(AlleleFreq(0.0).ln()),
            StrandBias::None,
            StrandBias::None,
            &observation,
        );
        assert_relative_eq!(*lh, 0.5f64.ln());

        let lh = model.likelihood_observation(
            LogProb(AlleleFreq(0.5).ln()),
            LogProb(AlleleFreq(0.5).ln()),
            StrandBias::None,
            StrandBias::None,
            &observation,
        );
        assert_relative_eq!(*lh, 0.5f64.ln());

        let lh = model.likelihood_observation(
            LogProb(AlleleFreq(0.1).ln()),
            LogProb(AlleleFreq(0.0).ln()),
            StrandBias::None,
            StrandBias::None,
            &observation,
        );
        assert_relative_eq!(*lh, 0.1f64.ln());

        // test with 50% purity
        let model = ContaminatedSampleLikelihoodModel::new(0.5);

        let lh = model.likelihood_observation(
            LogProb(AlleleFreq(0.0).ln()),
            LogProb(AlleleFreq(1.0).ln()),
            StrandBias::None,
            StrandBias::None,
            &observation,
        );
        assert_relative_eq!(*lh, 0.5f64.ln(), epsilon = 0.0000000001);
    }

    #[test]
    fn test_likelihood_observation_single_sample() {
        let model = SampleLikelihoodModel::new();

        let observation = observation(
            // prob_mapping
            LogProb::ln_one(),
            // prob_alt
            LogProb::ln_one(),
            // prob_ref
            LogProb::ln_zero(),
        );

        let lh = model.likelihood_observation(
            LogProb(AlleleFreq(1.0).ln()),
            StrandBias::None,
            &observation,
        );
        assert_relative_eq!(*lh, *LogProb::ln_one());

        let lh = model.likelihood_observation(
            LogProb(AlleleFreq(0.0).ln()),
            StrandBias::None,
            &observation,
        );
        assert_relative_eq!(*lh, *LogProb::ln_zero());

        let lh = model.likelihood_observation(
            LogProb(AlleleFreq(0.5).ln()),
            StrandBias::None,
            &observation,
        );
        assert_relative_eq!(*lh, 0.5f64.ln());

        let lh = model.likelihood_observation(
            LogProb(AlleleFreq(0.1).ln()),
            StrandBias::None,
            &observation,
        );
        assert_relative_eq!(*lh, 0.1f64.ln());
    }

    #[test]
    fn test_likelihood_pileup() {
        let model = ContaminatedSampleLikelihoodModel::new(1.0);
        let mut observations = Vec::new();
        for _ in 0..5 {
            observations.push(observation(
                LogProb::ln_one(),
                LogProb::ln_one(),
                LogProb::ln_zero(),
            ));
        }
        for _ in 0..5 {
            observations.push(observation(
                LogProb::ln_one(),
                LogProb::ln_zero(),
                LogProb::ln_one(),
            ));
        }
        let mut cache = likelihood::ContaminatedSampleCache::default();
        let lh = model.compute(
            &ContaminatedSampleEvent {
                primary: event(0.5),
                secondary: event(0.0),
            },
            &observations,
            &mut cache,
        );
        for af in linspace(0.0, 1.0, 10) {
            if af != 0.5 {
                let evt = ContaminatedSampleEvent {
                    primary: event(af),
                    secondary: event(0.0),
                };
                let l = model.compute(&evt, &observations, &mut cache);
                assert!(lh > l);
                assert!(cache.contains_key(&evt));
            }
        }
    }
}
