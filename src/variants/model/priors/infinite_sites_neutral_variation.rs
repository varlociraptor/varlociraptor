use std::f64;

use bio::stats::{LogProb, Prob};

/// The classical population genetic model used for variant calling in e.g. GATK and Samtools.
pub struct InfiniteSitesNeutralVariationModel {
    heterozygosity: LogProb,
    zero_prob: LogProb,
    max_m: u32,
}

impl InfiniteSitesNeutralVariationModel {
    pub fn new(n_samples: u32, ploidy: u32, heterozygosity: Prob) -> Self {
        let heterozygosity = LogProb::from(heterozygosity);

        let zero_prob = {
            let allele_freq_sum = (1..n_samples * ploidy + 1).fold(0.0, |s, m| s + 1.0 / m as f64);
            LogProb(*heterozygosity + allele_freq_sum.ln()).ln_one_minus_exp()
        };

        InfiniteSitesNeutralVariationModel {
            heterozygosity: heterozygosity,
            zero_prob: zero_prob,
            max_m: n_samples * ploidy,
        }
    }

    /// Prior probability for m alternative alleles.
    pub fn prior_prob(&self, m: u32) -> LogProb {
        if m == 0 {
            self.zero_prob
        } else {
            assert!(m <= self.max_m, "m too large (at most ploidy * n_samples)");
            LogProb(*self.heterozygosity - (m as f64).ln())
        }
    }
}
