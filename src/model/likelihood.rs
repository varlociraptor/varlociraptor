use bio::stats::LogProb;

use model::evidence::Observation;
use model::AlleleFreq;


/// Variant calling model, taking purity and allele frequencies into account.
#[derive(Clone, Copy, Debug)]
pub struct LatentVariableModel {
    /// Purity of the case sample.
    purity: Option<LogProb>
}


impl LatentVariableModel {

    /// Create new model.
    pub fn new(purity: f64) -> Self {
        assert!(purity > 0.0 && purity <= 1.0);
        LatentVariableModel { purity: Some(LogProb(purity.ln())) }
    }

    pub fn with_single_sample() -> Self {
        LatentVariableModel { purity: None }
    }

    /// Impurity of the case sample (fraction of control cells in the case sample).
    fn impurity(&self) -> LogProb {
        self.purity.unwrap().ln_one_minus_exp()
    }

    /// Likelihood to observe a read given allele frequencies for case and control.
    fn likelihood_observation(&self,
                       observation: &Observation,
                       allele_freq_case: AlleleFreq,
                       allele_freq_control: Option<AlleleFreq>,
                       prob_sample_alt: LogProb
    ) -> LogProb {
        // probability to sample observation: AF * placement induced probability
        let prob_sample_alt_case = LogProb(allele_freq_case.ln()) +
                                   prob_sample_alt;

        match (allele_freq_control, self.purity) {
            (Some(allele_freq_control), Some(purity)) => {
                // read comes from control sample and is correctly mapped
                let prob_sample_alt_control = LogProb(allele_freq_control.ln()) +
                                              prob_sample_alt; // TODO correct?
                let prob_control = self.impurity() +
                                   (prob_sample_alt_control + observation.prob_alt).ln_add_exp(
                    prob_sample_alt_control.ln_one_minus_exp() +
                    observation.prob_ref
                );
                assert!(!prob_control.is_nan());
                // read comes from case sample and is correctly mapped
                let prob_case = purity +
                                (prob_sample_alt_case + observation.prob_alt).ln_add_exp(
                    prob_sample_alt_case.ln_one_minus_exp() +
                    observation.prob_ref
                );
                assert!(!prob_case.is_nan());
                // total probability
                let total = (observation.prob_mapping + prob_control.ln_add_exp(prob_case)).ln_add_exp(
                    observation.prob_mapping.ln_one_minus_exp() +
                    observation.prob_mismapped
                );
                assert!(!total.is_nan());
                total
            },
            (None, purity) => {
                // no AF for control sample given
                if let Some(purity) = purity {
                    assert!(
                        purity == LogProb::ln_one(),
                        "no control allele frequency given but purity is not 1.0"
                    );
                }

                // read comes from case sample and is correctly mapped
                let prob_case = (prob_sample_alt_case + observation.prob_alt).ln_add_exp(
                    prob_sample_alt_case.ln_one_minus_exp() +
                    observation.prob_ref
                );
                assert!(!prob_case.is_nan());
                // total probability
                let total = (observation.prob_mapping + prob_case).ln_add_exp(
                    observation.prob_mapping.ln_one_minus_exp() +
                    observation.prob_mismapped
                );
                assert!(!total.is_nan());
                total
            }
            (Some(_), None) => {
                panic!("control allele frequency given but purity not defined")
            }
        }
    }

    /// Likelihood to observe a pileup given allele frequencies for case and control.
    pub fn likelihood_pileup(&self,
                             pileup: &[Observation],
                             allele_freq_case: AlleleFreq,
                             allele_freq_control: Option<AlleleFreq>,
                             prob_sample_alt: LogProb
    ) -> LogProb {
        // calculate product of per-read likelihoods in log space
        let likelihood = pileup.iter().fold(
            LogProb::ln_one(),
            |prob, obs| {
                prob + self.likelihood_observation(
                    obs, allele_freq_case, allele_freq_control, prob_sample_alt
                )
            }
        );
        assert!(!likelihood.is_nan());
        likelihood
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::evidence::{Observation, Evidence};
    use itertools_num::linspace;
    use bio::stats::LogProb;

    #[test]
    fn test_likelihood_observation() {
        let model = LatentVariableModel::new(1.0);
        let observation = Observation{
            prob_mapping: LogProb::ln_one(),
            prob_alt: LogProb::ln_one(),
            prob_ref: LogProb::ln_zero(),
            prob_sample_alt: LogProb::ln_one(),
            prob_mismapped: LogProb::ln_one(),
            evidence: Evidence::dummy_alignment()
        };

        let lh = model.likelihood_observation(&observation, LogProb::ln_one(), Some(LogProb::ln_zero()));
        assert_relative_eq!(*lh, *LogProb::ln_one());

        let lh = model.likelihood_observation(&observation, LogProb::ln_zero(), Some(LogProb::ln_zero()));
        assert_relative_eq!(*lh, *LogProb::ln_zero());

        let lh = model.likelihood_observation(&observation, LogProb(0.5f64.ln()), Some(LogProb::ln_zero()));
        assert_relative_eq!(*lh, 0.5f64.ln());

        let lh = model.likelihood_observation(&observation, LogProb(0.5f64.ln()), Some(LogProb(0.5f64.ln())));
        assert_relative_eq!(*lh, 0.5f64.ln());

        let lh = model.likelihood_observation(&observation, LogProb(0.1f64.ln()), Some(LogProb::ln_zero()));
        assert_relative_eq!(*lh, 0.1f64.ln());

        // test with 50% purity
        let model = LatentVariableModel::new(0.5);

        let lh = model.likelihood_observation(&observation, LogProb::ln_zero(), Some(LogProb::ln_one()));
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
                prob_sample_alt: LogProb::ln_one(),
                prob_mismapped: LogProb::ln_one(),
                evidence: Evidence::dummy_alignment()
            });
        }
        for _ in 0..5 {
            observations.push(Observation{
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_zero(),
                prob_ref: LogProb::ln_one(),
                prob_sample_alt: LogProb::ln_one(),
                prob_mismapped: LogProb::ln_one(),
                evidence: Evidence::dummy_alignment()
            });
        }
        let lh = model.likelihood_pileup(&observations, 0.5, Some(0.0));
        for af in linspace(0.0, 1.0, 10) {
            if af != 0.5 {
                let l = model.likelihood_pileup(&observations, af, Some(0.0));
                assert!(lh > l);
            }
        }
    }
}
