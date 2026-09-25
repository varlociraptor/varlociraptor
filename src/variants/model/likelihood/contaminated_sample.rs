use bio::stats::{bayesian::model::Likelihood, LogProb};
use lru::LruCache;
use itertools::Itertools;

use crate::variants::evidence::observations::pileup::Pileup;
use crate::variants::evidence::observations::read_observation::ProcessedReadObservation;
use crate::variants::model::bias::Artifacts;
use crate::variants::model::AlleleFreq;
use crate::variants::model::likelihood::{Event, likelihood_mapping};

pub(crate) type ContaminatedSampleCache = LruCache<ContaminatedSampleEvent, LogProb>;

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
        observation: &ProcessedReadObservation,
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

    /// Estimate the allele frequency as expected value of the posterior distribution with a uniform prior.
    pub(crate) fn estimate_allele_freq(&self, pileup: &Pileup, contaminant_allele_freq: AlleleFreq) -> AlleleFreq {
        // METHOD: Calculate expected depth and expected alt depth, and correct it with the purity.
        // This is a good estimate of the mode of the posterior distribution.
        let depth = LogProb::ln_sum_exp(&pileup.read_observations().iter().map(|obs| obs.prob_mapping()).collect_vec());
        let alt_depth = LogProb::ln_sum_exp(&pileup.read_observations().iter().map(|obs| obs.prob_mapping() + obs.prob_missed_allele + obs.prob_alt()).collect_vec());
        let expected_observed_allele_freq = (alt_depth - depth).exp();
        // METHOD: it holds purity * sample_vaf + impurity * contaminant_vaf = observed_vaf
        // We resolve this to sample_vaf below
        AlleleFreq((expected_observed_allele_freq - *contaminant_allele_freq * self.impurity.exp()) / self.purity.exp())
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
            let likelihood =
                pileup
                    .read_observations()
                    .iter()
                    .fold(LogProb::ln_one(), |prob, obs| {
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


#[cfg(test)]
mod tests {
    use super::*;
    use crate::variants::model::bias::Artifacts;
    use crate::variants::model::likelihood;
    use crate::variants::model::likelihood::contaminated_sample::ContaminatedSampleLikelihoodModel;
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
        let mut observations = Pileup::default();
        for _ in 0..10 {
            observations.read_observations_mut().push(observation(
                LogProb::ln_one(),
                LogProb::ln_zero(),
                LogProb::ln_one(),
            ));
        }
        let mut cache = ContaminatedSampleCache::new(100);

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
                .read_observations()
                .iter()
                .map(|observation| biases().prob_ref(observation))
                .sum::<LogProb>()
        );
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn test_likelihood_pileup() {
        let model = ContaminatedSampleLikelihoodModel::new(1.0);
        let mut observations = Pileup::default();
        for _ in 0..5 {
            observations.read_observations_mut().push(observation(
                LogProb::ln_one(),
                LogProb::ln_one(),
                LogProb::ln_zero(),
            ));
        }
        for _ in 0..5 {
            observations.read_observations_mut().push(observation(
                LogProb::ln_one(),
                LogProb::ln_zero(),
                LogProb::ln_one(),
            ));
        }
        let mut cache = ContaminatedSampleCache::new(100);
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
