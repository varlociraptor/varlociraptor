// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::HashMap;

use bio::stats::{bayesian::model::Likelihood, LogProb};

use crate::grammar::VAFRange;
use crate::utils::NUMERICAL_EPSILON;
use crate::variants::evidence::observation::{Observation, ReadPosition};
use crate::variants::model::bias::Biases;
use crate::variants::model::AlleleFreq;
use crate::variants::sample::Pileup;

pub(crate) type ContaminatedSampleCache = HashMap<ContaminatedSampleEvent, LogProb>;
pub(crate) type SingleSampleCache = HashMap<Event, LogProb>;

#[derive(PartialEq, Eq, Debug, Clone, Hash)]
pub(crate) enum Event {
    AlleleFreq(AlleleFreq),
    Biases(Biases),
}

impl Event {
    pub(crate) fn allele_freq(&self) -> Option<AlleleFreq> {
        if let Event::AlleleFreq(vaf) = self {
            Some(*vaf)
        } else {
            None
        }
    }

    pub(crate) fn is_artifact(&self) -> bool {
        if let Event::Biases(_) = self {
            true
        } else {
            false
        }
    }
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
    pub(crate) secondary: Option<Event>,
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
        allele_freq_secondary: Option<LogProb>,
        biases: &Biases,
        observation: &Observation<ReadPosition>,
    ) -> LogProb {
        // Obtain purity and impurity.
        // METHOD: Ignore purity/impurity in case of artifact event. The reason is
        // that artifacts happen during the sequencing and mapping process, and are not
        // propagated between samples via contamination.
        let (purity, impurity, allele_freq_secondary) =
            if let Some(allele_freq_secondary) = allele_freq_secondary {
                (self.purity, self.impurity, allele_freq_secondary)
            } else {
                (LogProb::ln_one(), LogProb::ln_zero(), LogProb::ln_zero())
            };

        // Step 1: likelihoods for the mapping case.
        // Case 1: read comes from primary sample and is correctly mapped
        let prob_primary = purity + likelihood_mapping(allele_freq_primary, biases, observation);
        // Case 2: read comes from secondary sample and is correctly mapped
        let prob_secondary =
            impurity + likelihood_mapping(allele_freq_secondary, biases, observation);

        // Step 4: total probability
        // Important note: we need to multiply a probability for a hypothetical missed allele
        // in the mismapping case. Otherwise, it can happen that mismapping dominates subtle
        // differences in the likelihood for alt and ref allele with low probabilities and very
        // low allele frequencies, such that we loose sensitivity for those.
        let total = (observation.prob_mapping() + prob_secondary.ln_add_exp(prob_primary))
            .ln_add_exp(
                observation.prob_mismapping()
                    + observation.prob_missed_allele
                    + biases.prob_any(observation),
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
        if cache.contains_key(events) {
            *cache.get(events).unwrap()
        } else {
            let likelihood = match events {
                ContaminatedSampleEvent { primary: Event::AlleleFreq(primary_vaf), secondary: Some(Event::AlleleFreq(secondary_vaf)) } => {
                    let ln_af_primary = LogProb(primary_vaf.ln());
                    let ln_af_secondary = LogProb(secondary_vaf.ln());
                    let biases = Biases::none();
                    pileup.iter().map(|obs| {
                        self.likelihood_observation(
                            ln_af_primary,
                            Some(ln_af_secondary),
                            &biases,
                            obs,
                        )
                    }).sum()
                }
                ContaminatedSampleEvent { primary: Event::Biases(biases), .. } => {
                    // METHOD: for biases, we integrate over [0,1], since any VAF is possible.
                    // Thereby, contamination is not considered, because the artifact indicated by the bias is generated
                    // during sequencing or mapping and thereby not propagated from sample to sample via contamination.
                    let (min_vaf, max_vaf) = VAFRange::present_observable_bounds(pileup.len());
                    let density = |i, vaf: f64| {
                        let vaf = LogProb(vaf.ln());
                        pileup.iter().map(|obs| {
                            self.likelihood_observation(
                                vaf,
                                None,
                                biases,
                                obs,
                            )
                        }).sum()
                    };

                    LogProb::ln_simpsons_integrate_exp(density, *min_vaf, *max_vaf, 11)
                },
                _ => unreachable!("bug: contaminated sample event where primary is not a bias but no vaf for secondary given"),
            };

            assert!(!likelihood.is_nan());
            cache.insert(events.clone(), likelihood);
            dbg!((events, likelihood));

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
        biases: &Biases,
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
    biases: &Biases,
    observation: &Observation<ReadPosition>,
) -> LogProb {
    // Step 1: calculate probability to sample from alt allele
    let prob_sample_alt = prob_sample_alt(observation, allele_freq);
    let prob_sample_ref = prob_sample_alt.ln_one_minus_exp();

    let prob_bias = biases.prob(observation);
    let prob_any_bias = biases.prob_any(observation);

    // Step 2: read comes from case sample and is correctly mapped
    let prob = LogProb::ln_sum_exp(&[
        // alt allele
        prob_sample_alt + prob_bias + observation.prob_alt,
        // ref allele (we don't care about the strand)
        prob_sample_ref + observation.prob_ref + prob_any_bias,
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
            let likelihood = match event {
                Event::AlleleFreq(allele_freq) => {
                    let ln_af = LogProb(allele_freq.ln());
                    let biases = Biases::none();
                    pileup
                        .iter()
                        .map(|obs| self.likelihood_observation(ln_af, &biases, obs))
                        .sum()
                }
                Event::Biases(biases) => {
                    // METHOD: for biases, we integrate over [0,1], since any VAF is possible.
                    let (min_vaf, max_vaf) = VAFRange::present_observable_bounds(pileup.len());
                    let density = |i, vaf: f64| {
                        let vaf = LogProb(vaf.ln());
                        pileup
                            .iter()
                            .map(|obs| self.likelihood_observation(vaf, biases, obs))
                            .sum()
                    };

                    LogProb::ln_simpsons_integrate_exp(density, *min_vaf, *max_vaf, 11)
                }
            };

            assert!(!likelihood.is_nan());

            cache.insert(event.clone(), likelihood);

            likelihood
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variants::model::bias::Biases;
    use crate::variants::model::likelihood;
    use crate::variants::model::tests::observation;
    use bio::stats::LogProb;
    use itertools_num::linspace;

    fn biases() -> Biases {
        Biases::none()
    }

    fn event(allele_freq: f64) -> Event {
        Event::AlleleFreq(AlleleFreq(allele_freq))
    }

    #[test]
    fn test_likelihood_observation_absent_single() {
        let observation = observation(LogProb::ln_one(), LogProb::ln_zero(), LogProb::ln_one());

        let model = SampleLikelihoodModel::new();

        let lh =
            model.likelihood_observation(LogProb(AlleleFreq(0.0).ln()), &biases(), &observation);
        assert_relative_eq!(*lh, *biases().prob_any(&observation));
    }

    #[test]
    fn test_likelihood_observation_absent() {
        let model = ContaminatedSampleLikelihoodModel::new(1.0);
        let observation = observation(LogProb::ln_one(), LogProb::ln_zero(), LogProb::ln_one());

        let lh = model.likelihood_observation(
            LogProb(AlleleFreq(0.0).ln()),
            Some(LogProb(AlleleFreq(0.0).ln())),
            &biases(),
            &observation,
        );
        assert_relative_eq!(*lh, *biases().prob_any(&observation));
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
                secondary: Some(event(0.0)),
            },
            &observations,
            &mut cache,
        );
        assert_relative_eq!(
            *lh,
            *observations
                .iter()
                .map(|observation| biases().prob_any(&observation))
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
        let mut cache = likelihood::SingleSampleCache::default();
        let evt = event(0.0);
        let lh = model.compute(&evt, &observations, &mut cache);
        assert_relative_eq!(
            *lh,
            *observations
                .iter()
                .map(|observation| biases().prob_any(&observation))
                .sum::<LogProb>()
        );
        assert!(cache.contains_key(&evt))
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
        let mut cache = likelihood::ContaminatedSampleCache::default();
        let lh = model.compute(
            &ContaminatedSampleEvent {
                primary: event(0.5),
                secondary: Some(event(0.0)),
            },
            &observations,
            &mut cache,
        );
        for af in linspace(0.0, 1.0, 10) {
            if af != 0.5 {
                let evt = ContaminatedSampleEvent {
                    primary: event(af),
                    secondary: Some(event(0.0)),
                };
                let l = model.compute(&evt, &observations, &mut cache);
                assert!(lh > l);
                assert!(cache.contains_key(&evt));
            }
        }
    }
}
