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
    pub fn new(purity: f64) -> Self {
        assert!(purity > 0.0 && purity <= 1.0);
        LatentVariableModel { purity: LogProb(purity.ln()) }
    }

    /// Impurity of the case sample (fraction of control cells in the case sample).
    fn impurity(&self) -> LogProb {
        self.purity.ln_one_minus_exp()
    }

    /// Likelihood to observe a read given allele frequencies for case and control.
    fn likelihood_observation(&self,
                       observation: &Observation,
                       allele_freq_case: LogProb,
                       allele_freq_control: LogProb) -> LogProb {
        // read comes from control sample and is correctly mapped
        let prob_control = self.impurity() +
                           (allele_freq_control + observation.prob_alt).ln_add_exp(
                                 allele_freq_control.ln_one_minus_exp() + observation.prob_ref
                           );
        // read comes from case sample and is correctly mapped
        let prob_case = self.purity +
                        (allele_freq_case + observation.prob_alt).ln_add_exp(
                                      allele_freq_case.ln_one_minus_exp() + observation.prob_ref
                        );
        // total probability
        let total = (observation.prob_mapping + prob_control.ln_add_exp(prob_case)).ln_add_exp(
                          observation.prob_mapping.ln_one_minus_exp() + observation.prob_mismapped
                    );
        total
    }

    /// Likelihood to observe a pileup given allele frequencies for case and control.
    #[cfg_attr(feature="flame_it", flame)]
    pub fn likelihood_pileup(&self,
                             pileup: &[Observation],
                             allele_freq_case: f64,
                             allele_freq_control: f64) -> LogProb {
        let allele_freq_case = LogProb(allele_freq_case.ln());
        let allele_freq_control = LogProb(allele_freq_control.ln());
        // calculate product of per-read likelihoods in log space
        let likelihood = pileup.iter().fold(LogProb::ln_one(),
            |prob, obs| prob + self.likelihood_observation(obs, allele_freq_case, allele_freq_control));
        likelihood
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::sample::Observation;
    use itertools::linspace;
    use bio::stats::LogProb;

    #[test]
    fn test_likelihood_observation() {
        let model = LatentVariableModel::new(1.0);
        let observation = Observation{
            prob_mapping: LogProb::ln_one(),
            prob_alt: LogProb::ln_one(),
            prob_ref: LogProb::ln_zero(),
            prob_mismapped: LogProb::ln_one()
        };

        let lh = model.likelihood_observation(&observation, LogProb::ln_one(), LogProb::ln_zero());
        assert_relative_eq!(*lh, *LogProb::ln_one());

        let lh = model.likelihood_observation(&observation, LogProb::ln_zero(), LogProb::ln_zero());
        assert_relative_eq!(*lh, *LogProb::ln_zero());

        let lh = model.likelihood_observation(&observation, LogProb(0.5f64.ln()), LogProb::ln_zero());
        assert_relative_eq!(*lh, 0.5f64.ln());

        let lh = model.likelihood_observation(&observation, LogProb(0.5f64.ln()), LogProb(0.5f64.ln()));
        assert_relative_eq!(*lh, 0.5f64.ln());

        let lh = model.likelihood_observation(&observation, LogProb(0.1f64.ln()), LogProb::ln_zero());
        assert_relative_eq!(*lh, 0.1f64.ln());

        // test with 50% purity
        let model = LatentVariableModel::new(0.5);

        let lh = model.likelihood_observation(&observation, LogProb::ln_zero(), LogProb::ln_one());
        assert_relative_eq!(*lh, 0.5f64.ln());
    }

    #[test]
    fn test_likelihood_pileup() {
        let model = LatentVariableModel::new(1.0);
        let mut observations = Vec::new();
        for _ in 0..5 {
            observations.push(Observation{
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_one(),
                prob_ref: LogProb::ln_zero(),
                prob_mismapped: LogProb::ln_one()
            });
        }
        for _ in 0..5 {
            observations.push(Observation{
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_zero(),
                prob_ref: LogProb::ln_one(),
                prob_mismapped: LogProb::ln_one()
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
