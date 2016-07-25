use bio::stats::logprobs;
use bio::stats::LogProb;

use model::observations::Observation;


pub type AlleleFreq = f64;


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
