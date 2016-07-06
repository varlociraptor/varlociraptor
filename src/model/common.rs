use bio::stats::logprobs;
use bio::stats::LogProb;


pub type AlleleFreq = f64;


pub trait SequenceRead {
    /// Posterior probability that the read has been mapped correctly (1 - MAPQ).
    fn prob_mapping(&self) -> LogProb;
    /// Probability that the read comes from the alternative allele.
    fn prob_alt(&self) -> LogProb;
    /// Probability that the read comes from the reference allele.
    fn prob_ref(&self) -> LogProb;
    /// Probability of the read given that it has been mismapped.
    fn prob_mismapped(&self) -> LogProb;
}


pub trait LatentVariableModel {
    /// Purity of the case sample.
    fn purity(&self) -> LogProb;

    /// Impurity of the case sample (fraction of control cells in the case sample).
    fn impurity(&self) -> LogProb {
        logprobs::ln_1m_exp(self.purity())
    }

    /// Likelihood to observe a read given allele frequencies for case and control.
    fn likelihood_read(&self,
                       read: &Box<SequenceRead>,
                       allele_freq_case: AlleleFreq,
                       allele_freq_control: AlleleFreq) -> LogProb {
        // read comes from control sample and is correctly mapped
        let prob_control = self.impurity() +
                           logprobs::add(allele_freq_control.ln() + read.prob_alt(),
                                         (1.0 - allele_freq_control).ln() + read.prob_ref());
        // read comes from case sample and is correctly mapped
        let prob_case = self.purity() +
                        logprobs::add(allele_freq_case.ln() + read.prob_alt(),
                                      (1.0 - allele_freq_case).ln() + read.prob_ref());
        // total probability
        let total = logprobs::add(read.prob_mapping() + logprobs::add(prob_control, prob_case),
                                  logprobs::ln_1m_exp(read.prob_mapping()) + read.prob_mismapped());
        total
    }

    /// Likelihood to observe a pileup given allele frequencies for case and control.
    fn likelihood_pileup<'a, I: IntoIterator<Item=&'a Box<SequenceRead>>>(&self,
                                                                 pileup: I,
                                                                 allele_freq_case: AlleleFreq,
                                                                 allele_freq_control: AlleleFreq) -> LogProb {
        // calculate product of per-read likelihoods in log space
        let likelihood = pileup.into_iter().fold(0.0,
            |prob, read| prob + self.likelihood_read(read, allele_freq_case, allele_freq_control));
        likelihood
    }

    fn prior_prob(&self, allele_freq_case: LogProb, allele_freq_control: LogProb);
}
