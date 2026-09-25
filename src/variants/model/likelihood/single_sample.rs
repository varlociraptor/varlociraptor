use bio::stats::{bayesian::model::Likelihood, LogProb};
use lru::LruCache;
use itertools::Itertools;

use crate::variants::evidence::observations::pileup::Pileup;
use crate::variants::evidence::observations::read_observation::ProcessedReadObservation;
use crate::variants::model::bias::Artifacts;
use crate::variants::model::AlleleFreq;
use crate::variants::model::likelihood::{Event, likelihood_mapping};

/// Likelihood model for single sample.

pub(crate) type SingleSampleCache = LruCache<Event, LogProb>;

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
        observation: &ProcessedReadObservation,
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

    /// Estimate the allele frequency as expected value of the posterior distribution with a uniform prior.
    pub(crate) fn estimate_allele_freq(&self, pileup: &Pileup) -> AlleleFreq {
        // Calculate expected depth
        let depth = LogProb::ln_sum_exp(&pileup.read_observations().iter().map(|obs| obs.prob_mapping()).collect_vec());
        let alt_depth = LogProb::ln_sum_exp(&pileup.read_observations().iter().map(|obs| obs.prob_mapping() + obs.prob_missed_allele + obs.prob_alt()).collect_vec());
        AlleleFreq((alt_depth - depth).exp())
    }
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
            let likelihood =
                pileup
                    .read_observations()
                    .iter()
                    .fold(LogProb::ln_one(), |prob, obs| {
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
    use crate::variants::model::likelihood::single_sample::SampleLikelihoodModel;
    use crate::variants::model::tests::observation;
    use bio::stats::LogProb;

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
    fn test_likelihood_pileup_absent_single() {
        let model = SampleLikelihoodModel::new();
        let mut observations = Pileup::default();
        for _ in 0..10 {
            observations.read_observations_mut().push(observation(
                LogProb::ln_one(),
                LogProb::ln_zero(),
                LogProb::ln_one(),
            ));
        }
        let mut cache = SingleSampleCache::new(100);
        let evt = event(0.0);
        let lh = model.compute(&evt, &observations, &mut cache);
        assert_relative_eq!(
            *lh,
            *observations
                .read_observations()
                .iter()
                .map(|observation| biases().prob_ref(observation))
                .sum::<LogProb>()
        );
        assert!(cache.get(&evt).is_some())
    }
}
