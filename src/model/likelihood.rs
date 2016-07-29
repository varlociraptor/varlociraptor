use bio::stats::logprobs;
use bio::stats::LogProb;

use model::sample::Observation;


pub type AlleleFreq = f64;


/// Variant calling model, taking purity and allele frequencies into account.
pub struct LatentVariableModel {
    /// Purity of the case sample.
    purity: LogProb
}


impl LatentVariableModel {

    /// Create new model.
    pub fn new(purity: LogProb) -> Self {
        LatentVariableModel { purity: purity }
    }

    /// Impurity of the case sample (fraction of control cells in the case sample).
    fn impurity(&self) -> LogProb {
        logprobs::ln_1m_exp(self.purity)
    }

    /// Likelihood to observe a read given allele frequencies for case and control.
    fn likelihood_observation(&self,
                       observation: &Observation,
                       allele_freq_case: f64,
                       allele_freq_control: f64) -> LogProb {
        // read comes from control sample and is correctly mapped
        let prob_control = self.impurity() +
                           logprobs::add(allele_freq_control.ln() + observation.prob_alt,
                                         (1.0 - allele_freq_control).ln() + observation.prob_ref);
        // read comes from case sample and is correctly mapped
        let prob_case = self.purity +
                        logprobs::add(allele_freq_case.ln() + observation.prob_alt,
                                      (1.0 - allele_freq_case).ln() + observation.prob_ref);
        // total probability
        let total = logprobs::add(observation.prob_mapping + logprobs::add(prob_control, prob_case),
                                  logprobs::ln_1m_exp(observation.prob_mapping) + observation.prob_mismapped);
        total
    }

    /// Likelihood to observe a pileup given allele frequencies for case and control.
    pub fn likelihood_pileup(&self,
                             pileup: &[Observation],
                             allele_freq_case: f64,
                             allele_freq_control: f64) -> LogProb {
        // calculate product of per-read likelihoods in log space
        let likelihood = pileup.iter().fold(0.0,
            |prob, obs| prob + self.likelihood_observation(obs, allele_freq_case, allele_freq_control));
        likelihood
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::sample::Observation;
    use itertools::linspace;

    #[test]
    fn test_likelihood_observation() {
        let model = LatentVariableModel::new(1.0f64.ln());
        let observation = Observation{
            prob_mapping: 1.0f64.ln(),
            prob_alt: 1.0f64.ln(),
            prob_ref: 0.0f64.ln(),
            prob_mismapped: 1.0f64.ln()
        };

        let lh = model.likelihood_observation(&observation, 1.0, 0.0);
        assert_relative_eq!(lh, 1.0f64.ln());

        let lh = model.likelihood_observation(&observation, 0.0, 0.0);
        assert_relative_eq!(lh, 0.0f64.ln());

        let lh = model.likelihood_observation(&observation, 0.5, 0.0);
        assert_relative_eq!(lh, 0.5f64.ln());

        let lh = model.likelihood_observation(&observation, 0.5, 0.5);
        assert_relative_eq!(lh, 0.5f64.ln());

        let lh = model.likelihood_observation(&observation, 0.1, 0.0);
        assert_relative_eq!(lh, 0.1f64.ln());

        // test with 50% purity
        let model = LatentVariableModel::new(0.5f64.ln());

        let lh = model.likelihood_observation(&observation, 0.0, 1.0);
        assert_relative_eq!(lh, 0.5f64.ln());
    }

    #[test]
    fn test_likelihood_pileup() {
        let model = LatentVariableModel::new(1.0f64.ln());
        let mut observations = Vec::new();
        for _ in 0..5 {
            observations.push(Observation{
                prob_mapping: 1.0f64.ln(),
                prob_alt: 1.0f64.ln(),
                prob_ref: 0.0f64.ln(),
                prob_mismapped: 1.0f64.ln()
            });
        }
        for _ in 0..5 {
            observations.push(Observation{
                prob_mapping: 1.0f64.ln(),
                prob_alt: 0.0f64.ln(),
                prob_ref: 1.0f64.ln(),
                prob_mismapped: 1.0f64.ln()
            });
        }
        let lh = model.likelihood_pileup(&observations, 0.5, 0.0);
        for af in linspace(0.0, 1.0, 10) {
            if af != 0.5 {
                let l = model.likelihood_pileup(&observations, af, 0.0);
                assert!(lh > l);
            }
        }
    }
}
