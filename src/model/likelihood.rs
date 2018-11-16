use bio::stats::LogProb;

use model::evidence::Observation;
use model::AlleleFreq;

/// Variant calling model, taking purity and allele frequencies into account.
#[derive(Clone, Copy, Debug)]
pub struct LatentVariableModel {
    /// Purity of the case sample.
    purity: Option<LogProb>,
    impurity: Option<LogProb>,
}

impl LatentVariableModel {
    /// Create new model.
    pub fn new(purity: f64) -> Self {
        assert!(purity > 0.0 && purity <= 1.0);
        LatentVariableModel {
            purity: Some(LogProb(purity.ln())),
            impurity: Some(LogProb(purity.ln()).ln_one_minus_exp()),
        }
    }

    pub fn with_single_sample() -> Self {
        LatentVariableModel {
            purity: None,
            impurity: None,
        }
    }

    /// Likelihood to observe a read given allele frequencies for case and control.
    fn likelihood_observation_case_control(
        &self,
        observation: &Observation,
        allele_freq_case: LogProb,
        allele_freq_control: LogProb,
    ) -> LogProb {
        // Step 1: probability to sample observation: AF * placement induced probability
        let prob_sample_alt_case = allele_freq_case + observation.prob_sample_alt;
        let prob_sample_alt_control = allele_freq_control + observation.prob_sample_alt;

        // Step 2: read comes from control sample and is correctly mapped
        let prob_control = self.impurity.unwrap()
            + (prob_sample_alt_control + observation.prob_alt)
                .ln_add_exp(prob_sample_alt_control.ln_one_minus_exp() + observation.prob_ref);
        assert!(!prob_control.is_nan());

        // Step 3: read comes from case sample and is correctly mapped
        println!(
            "{} {} {}",
            *prob_sample_alt_case, *observation.prob_sample_alt, *allele_freq_case
        );
        let prob_case = self.purity.unwrap()
            + (prob_sample_alt_case + observation.prob_alt)
                .ln_add_exp(prob_sample_alt_case.ln_one_minus_exp() + observation.prob_ref);
        assert!(!prob_case.is_nan());

        // Step 4: total probability
        let total = (observation.prob_mapping + prob_control.ln_add_exp(prob_case))
            .ln_add_exp(observation.prob_mismapping);
        //println!("afs {}:{}, case {:?} control {:?} total {:?}", allele_freq_case, allele_freq_control, prob_case, prob_control, total);
        assert!(!total.is_nan());
        total
    }

    /// Likelihood to observe a read given allele frequency for a single sample.
    fn likelihood_observation_single_sample(
        observation: &Observation,
        allele_freq_case: LogProb,
    ) -> LogProb {
        // Step 1: calculate probability to sample from alt allele
        let prob_sample_alt = allele_freq_case + observation.prob_sample_alt;

        // Step 2: read comes from case sample and is correctly mapped
        let prob_case = (prob_sample_alt + observation.prob_alt)
            .ln_add_exp(prob_sample_alt.ln_one_minus_exp() + observation.prob_ref);
        assert!(!prob_case.is_nan());

        // Step 3: total probability
        let total = (observation.prob_mapping + prob_case).ln_add_exp(observation.prob_mismapping);
        assert!(!total.is_nan());
        total
    }

    /// Likelihood to observe a pileup given allele frequencies for case and control.
    pub fn likelihood_pileup(
        &self,
        pileup: &[Observation],
        allele_freq_case: AlleleFreq,
        allele_freq_control: Option<AlleleFreq>,
    ) -> LogProb {
        let likelihood;
        match allele_freq_control {
            Some(allele_freq_control) => {
                let ln_af_case = LogProb(allele_freq_case.ln());
                let ln_af_control = LogProb(allele_freq_control.ln());
                // calculate product of per-read likelihoods in log space
                likelihood = pileup.iter().fold(LogProb::ln_one(), |prob, obs| {
                    let lh =
                        self.likelihood_observation_case_control(obs, ln_af_case, ln_af_control);
                    prob + lh
                });
            }
            None => {
                // no AF for control sample given
                if let Some(purity) = self.purity {
                    assert!(
                        purity == LogProb::ln_one(),
                        "no control allele frequency given but purity is not 1.0"
                    );
                }
                let ln_af_case = LogProb(allele_freq_case.ln());
                // calculate product of per-read likelihoods in log space
                likelihood = pileup.iter().fold(LogProb::ln_one(), |prob, obs| {
                    let lh = Self::likelihood_observation_single_sample(obs, ln_af_case);
                    prob + lh
                });
            }
        }
        assert!(!likelihood.is_nan());
        likelihood
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::stats::LogProb;
    use itertools_num::linspace;
    use model::tests::observation;

    #[test]
    fn test_likelihood_observation_absent_single() {
        let observation = observation(LogProb::ln_one(), LogProb::ln_zero(), LogProb::ln_one());

        let lh = LatentVariableModel::likelihood_observation_single_sample(
            &observation,
            LogProb(AlleleFreq(0.0).ln()),
        );
        assert_relative_eq!(*lh, *LogProb::ln_one());
    }

    #[test]
    fn test_likelihood_observation_absent() {
        let model = LatentVariableModel::new(1.0);
        let observation = observation(LogProb::ln_one(), LogProb::ln_zero(), LogProb::ln_one());

        let lh = model.likelihood_observation_case_control(
            &observation,
            LogProb(AlleleFreq(0.0).ln()),
            LogProb(AlleleFreq(0.0).ln()),
        );
        assert_relative_eq!(*lh, *LogProb::ln_one());
    }

    #[test]
    fn test_likelihood_pileup_absent() {
        let model = LatentVariableModel::new(1.0);
        let mut observations = Vec::new();
        for _ in 0..10 {
            observations.push(observation(
                LogProb::ln_one(),
                LogProb::ln_zero(),
                LogProb::ln_one(),
            ));
        }

        let lh = model.likelihood_pileup(&observations, AlleleFreq(0.0), Some(AlleleFreq(0.0)));
        assert_relative_eq!(*lh, *LogProb::ln_one());
    }

    #[test]
    fn test_likelihood_pileup_absent_single() {
        let model = LatentVariableModel::new(1.0);
        let mut observations = Vec::new();
        for _ in 0..10 {
            observations.push(observation(
                LogProb::ln_one(),
                LogProb::ln_zero(),
                LogProb::ln_one(),
            ));
        }

        let lh = model.likelihood_pileup(&observations, AlleleFreq(0.0), None);
        assert_relative_eq!(*lh, *LogProb::ln_one());
    }

    #[test]
    fn test_likelihood_observation_case_control() {
        let model = LatentVariableModel::new(1.0);
        let observation = observation(LogProb::ln_one(), LogProb::ln_one(), LogProb::ln_zero());

        let lh = model.likelihood_observation_case_control(
            &observation,
            LogProb(AlleleFreq(1.0).ln()),
            LogProb(AlleleFreq(0.0).ln()),
        );
        assert_relative_eq!(*lh, *LogProb::ln_one());

        let lh = model.likelihood_observation_case_control(
            &observation,
            LogProb(AlleleFreq(0.0).ln()),
            LogProb(AlleleFreq(0.0).ln()),
        );
        assert_relative_eq!(*lh, *LogProb::ln_zero());

        let lh = model.likelihood_observation_case_control(
            &observation,
            LogProb(AlleleFreq(0.5).ln()),
            LogProb(AlleleFreq(0.0).ln()),
        );
        assert_relative_eq!(*lh, 0.5f64.ln());

        let lh = model.likelihood_observation_case_control(
            &observation,
            LogProb(AlleleFreq(0.5).ln()),
            LogProb(AlleleFreq(0.5).ln()),
        );
        assert_relative_eq!(*lh, 0.5f64.ln());

        let lh = model.likelihood_observation_case_control(
            &observation,
            LogProb(AlleleFreq(0.1).ln()),
            LogProb(AlleleFreq(0.0).ln()),
        );
        assert_relative_eq!(*lh, 0.1f64.ln());

        // test with 50% purity
        let model = LatentVariableModel::new(0.5);

        let lh = model.likelihood_observation_case_control(
            &observation,
            LogProb(AlleleFreq(0.0).ln()),
            LogProb(AlleleFreq(1.0).ln()),
        );
        assert_relative_eq!(*lh, 0.5f64.ln(), epsilon = 0.0000000001);
    }

    #[test]
    fn test_likelihood_observation_single_sample() {
        let observation = observation(
            // prob_mapping
            LogProb::ln_one(),
            // prob_alt
            LogProb::ln_one(),
            // prob_ref
            LogProb::ln_zero(),
        );

        let lh = LatentVariableModel::likelihood_observation_single_sample(
            &observation,
            LogProb(AlleleFreq(1.0).ln()),
        );
        assert_relative_eq!(*lh, *LogProb::ln_one());

        let lh = LatentVariableModel::likelihood_observation_single_sample(
            &observation,
            LogProb(AlleleFreq(0.0).ln()),
        );
        assert_relative_eq!(*lh, *LogProb::ln_zero());

        let lh = LatentVariableModel::likelihood_observation_single_sample(
            &observation,
            LogProb(AlleleFreq(0.5).ln()),
        );
        assert_relative_eq!(*lh, 0.5f64.ln());

        let lh = LatentVariableModel::likelihood_observation_single_sample(
            &observation,
            LogProb(AlleleFreq(0.1).ln()),
        );
        assert_relative_eq!(*lh, 0.1f64.ln());
    }

    #[test]
    fn test_likelihood_pileup() {
        let model = LatentVariableModel::new(1.0);
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
        let lh = model.likelihood_pileup(&observations, AlleleFreq(0.5), Some(AlleleFreq(0.0)));
        for af in linspace(0.0, 1.0, 10) {
            if af != 0.5 {
                let l =
                    model.likelihood_pileup(&observations, AlleleFreq(af), Some(AlleleFreq(0.0)));
                assert!(lh > l);
            }
        }
    }
}
