// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use bio::stats::{bayesian::model::Likelihood, LogProb};
use lru::LruCache;

use crate::utils::NUMERICAL_EPSILON;
use crate::variants::evidence::observation::{Observation, ReadPosition};
use crate::variants::model::bias::Artifacts;
use crate::variants::model::AlleleFreq;
use crate::variants::sample::Pileup;

pub(crate) type ContaminatedSampleCache = LruCache<ContaminatedSampleEvent, LogProb>;
pub(crate) type SingleSampleCache = LruCache<Event, LogProb>;

#[derive(PartialEq, Eq, Debug, Clone, Hash)]
pub(crate) struct Event {
    pub(crate) allele_freq: AlleleFreq,
    pub(crate) artifacts: Artifacts,
    pub(crate) is_discrete: bool,
}

fn prob_sample_alt(observation: &Observation<ReadPosition>, allele_freq: LogProb) -> LogProb {
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

pub(crate) trait ContaminatedSamplePairView<T> {
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

#[derive(PartialEq, Eq, Debug, Clone, Hash)]
pub(crate) struct ContaminatedSampleEvent {
    pub(crate) primary: Event,
    pub(crate) secondary: Event,
}

/// Variant calling model, taking purity and allele frequencies into account.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ContaminatedSampleLikelihoodModel {
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
    pub(crate) fn new(purity: f64) -> Self {
        assert!(purity > 0.0 && purity <= 1.0);
        let purity = LogProb(purity.ln());
        ContaminatedSampleLikelihoodModel {
            purity,
            impurity: purity.ln_one_minus_exp(),
        }
    }

    fn likelihood_observation(
        &self,
        allele_freq_primary: LogProb,
        allele_freq_secondary: LogProb,
        biases_primary: &Artifacts,
        biases_secondary: &Artifacts,
        observation: &Observation<ReadPosition>,
    ) -> LogProb {
        // Step 1: likelihoods for the mapping case.
        // Case 1: read comes from primary sample and is correctly mapped
        let prob_primary =
            self.purity + likelihood_mapping(allele_freq_primary, biases_primary, observation);
        // Case 2: read comes from secondary sample and is correctly mapped
        let prob_secondary = self.impurity
            + likelihood_mapping(allele_freq_secondary, biases_secondary, observation);

        // Step 4: total probability
        // Important note: we need to multiply a probability for a hypothetical missed allele
        // in the mismapping case. Otherwise, it can happen that mismapping dominates subtle
        // differences in the likelihood for alt and ref allele with low probabilities and very
        // low allele frequencies, such that we loose sensitivity for those.
        let total = (observation.prob_mapping() + prob_secondary.ln_add_exp(prob_primary))
            .ln_add_exp(
                observation.prob_mismapping()
                    + observation.prob_missed_allele
                    + biases_primary.prob_any(observation),
            );
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
        if let Some(prob) = cache.get(events) {
            *prob
        } else {
            let ln_af_primary = LogProb(events.primary.allele_freq.ln());
            let ln_af_secondary = LogProb(events.secondary.allele_freq.ln());

            // calculate product of per-observation likelihoods in log space
            let likelihood = pileup.iter().fold(LogProb::ln_one(), |prob, obs| {
                let lh = self.likelihood_observation(
                    ln_af_primary,
                    ln_af_secondary,
                    &events.primary.artifacts,
                    &events.secondary.artifacts,
                    obs,
                );
                prob + lh
            });

            assert!(!likelihood.is_nan());

            // METHOD: No caching for events with continuous VAFs as they are unlikely to reoccur.
            cache.put(events.clone(), likelihood);

            likelihood
        }
    }
}

/// Likelihood model for single sample.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct SampleLikelihoodModel {}

impl SampleLikelihoodModel {
    /// Create new model.
    pub(crate) fn new() -> Self {
        SampleLikelihoodModel {}
    }

    /// Likelihood to observe a read given allele frequency for a single sample.
    fn likelihood_observation(
        &self,
        allele_freq: LogProb,
        biases: &Artifacts,
        observation: &Observation<ReadPosition>,
    ) -> LogProb {
        // Step 1: likelihood for the mapping case.
        let prob = likelihood_mapping(allele_freq, biases, observation);

        // Step 2: total probability
        // Important note: we need to multiply a probability for a hypothetical missed allele
        // in the mismapping case. Otherwise, it can happen that mismapping dominates subtle
        // differences in the likelihood for alt and ref allele with low probabilities and very
        // low allele frequencies, such that we loose sensitivity for those.
        let total = (observation.prob_mapping() + prob).ln_add_exp(
            observation.prob_mismapping()
                + observation.prob_missed_allele
                + biases.prob_any(observation),
        );
        assert!(!total.is_nan());
        total
    }
}

/// Calculate likelihood of allele freq given observation in a single sample assuming that the
/// underlying fragment/read is mapped correctly.
fn likelihood_mapping(
    allele_freq: LogProb,
    biases: &Artifacts,
    observation: &Observation<ReadPosition>,
) -> LogProb {
    // Step 1: calculate probability to sample from alt allele
    let prob_sample_alt = prob_sample_alt(observation, allele_freq);
    let prob_sample_ref = prob_sample_alt.ln_one_minus_exp();

    let prob_bias_alt = biases.prob_alt(observation);
    let prob_bias_ref = biases.prob_ref(observation);

    // Step 2: read comes from case sample and is correctly mapped
    let prob = LogProb::ln_sum_exp(&[
        // alt allele
        prob_sample_alt + prob_bias_alt + observation.prob_alt,
        // ref allele
        prob_sample_ref + observation.prob_ref + prob_bias_ref,
    ]);
    assert!(!prob.is_nan());

    prob
}

impl Likelihood<SingleSampleCache> for SampleLikelihoodModel {
    type Event = Event;
    type Data = Pileup;

    /// Likelihood to observe a pileup given allele frequencies for case and control.
    fn compute(&self, event: &Event, pileup: &Pileup, cache: &mut SingleSampleCache) -> LogProb {
        if let Some(prob) = cache.get(event) {
            *prob
        } else {
            let ln_af = LogProb(event.allele_freq.ln());

            // calculate product of per-read likelihoods in log space
            let likelihood = pileup.iter().fold(LogProb::ln_one(), |prob, obs| {
                let lh = self.likelihood_observation(ln_af, &event.artifacts, obs);
                prob + lh
            });

            assert!(!likelihood.is_nan());

            cache.put(event.clone(), likelihood);

            likelihood
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variants::model::bias::Artifacts;
    use crate::variants::model::likelihood;
    use crate::variants::model::tests::observation;
    use bio::stats::LogProb;
    use itertools_num::linspace;

    fn biases() -> Artifacts {
        Artifacts::none()
    }

    fn event(allele_freq: f64) -> Event {
        Event {
            allele_freq: AlleleFreq(allele_freq),
            artifacts: biases(),
            is_discrete: true,
        }
    }

    #[test]
    fn test_likelihood_observation_absent_single() {
        let observation = observation(LogProb::ln_one(), LogProb::ln_zero(), LogProb::ln_one());

        let model = SampleLikelihoodModel::new();

        let lh =
            model.likelihood_observation(LogProb(AlleleFreq(0.0).ln()), &biases(), &observation);
        assert_relative_eq!(*lh, *biases().prob_ref(&observation));
    }

    #[test]
    fn test_likelihood_observation_absent() {
        let model = ContaminatedSampleLikelihoodModel::new(1.0);
        let observation = observation(LogProb::ln_one(), LogProb::ln_zero(), LogProb::ln_one());

        let lh = model.likelihood_observation(
            LogProb(AlleleFreq(0.0).ln()),
            LogProb(AlleleFreq(0.0).ln()),
            &biases(),
            &biases(),
            &observation,
        );
        assert_relative_eq!(*lh, *biases().prob_ref(&observation));
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
        let mut cache = likelihood::ContaminatedSampleCache::new(100);

        let lh = model.compute(
            &ContaminatedSampleEvent {
                primary: event(0.0),
                secondary: event(0.0),
            },
            &observations,
            &mut cache,
        );
        assert_relative_eq!(
            *lh,
            *observations
                .iter()
                .map(|observation| biases().prob_ref(&observation))
                .sum::<LogProb>()
        );
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
        let mut cache = likelihood::SingleSampleCache::new(100);
        let evt = event(0.0);
        let lh = model.compute(&evt, &observations, &mut cache);
        assert_relative_eq!(
            *lh,
            *observations
                .iter()
                .map(|observation| biases().prob_ref(&observation))
                .sum::<LogProb>()
        );
        assert!(cache.get(&evt).is_some())
    }

    #[test]
    #[allow(clippy::float_cmp)]
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
        let mut cache = likelihood::ContaminatedSampleCache::new(100);
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
                assert!(cache.get(&evt).is_some());
            }
        }
    }
}
