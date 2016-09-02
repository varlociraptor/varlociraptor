use std::f64;

use itertools::Itertools;
use bio::stats::LogProb;

use model::{Variant, DiscreteAlleleFreqs, AlleleFreq};
use priors::Model;


/// The classical population genetic model used for variant calling in e.g. GATK and Samtools.
pub struct InfiniteSitesNeutralVariationModel {
    ploidy: u32,
    heterozygosity: LogProb,
    zero_prob: LogProb,
    allele_freqs: DiscreteAlleleFreqs
}


impl InfiniteSitesNeutralVariationModel {
    /// Create new model for given ploidy and heterozygosity.
    pub fn new(ploidy: u32, heterozygosity: f64) -> Self {
        let heterozygosity = LogProb(heterozygosity.ln());
        let zero_prob = LogProb(*heterozygosity +
            (1..ploidy + 1).fold(0.0, |s, m| s + 1.0 / m as f64).ln()
        ).ln_one_minus_exp();

        let allele_freqs = (0..ploidy + 1).map(|m| AlleleFreq(m as f64 / ploidy as f64)).collect_vec();

        InfiniteSitesNeutralVariationModel {
            ploidy: ploidy,
            heterozygosity: heterozygosity,
            zero_prob: zero_prob,
            allele_freqs: allele_freqs
        }
    }
}


impl Model<DiscreteAlleleFreqs> for InfiniteSitesNeutralVariationModel {
    fn prior_prob(&self, af: AlleleFreq, _: Variant) -> LogProb {
        if *af > 0.0 {
            let m = *af * self.ploidy as f64;
            if relative_eq!(m % 1.0, 0.0) {
                // if m is discrete
                LogProb(*self.heterozygosity - m.ln())
            } else {
                // invalid allele frequency
                LogProb::ln_zero()
            }
        } else {
            self.zero_prob
        }
    }

    fn joint_prob<L>(&self, afs: &DiscreteAlleleFreqs, likelihood: &L, variant: Variant) -> LogProb where
        L: Fn(AlleleFreq) -> LogProb
    {
        let summands = afs.iter().map(|af| self.prior_prob(*af, variant) + likelihood(*af)).collect_vec();
        LogProb::ln_sum_exp(&summands)
    }

    fn marginal_prob<L>(&self, likelihood: &L, variant: Variant) -> LogProb where
        L: Fn(AlleleFreq) -> LogProb
    {
        self.joint_prob(self.allele_freqs(), likelihood, variant)
    }

    fn allele_freqs(&self) -> &DiscreteAlleleFreqs {
        &self.allele_freqs
    }
}
