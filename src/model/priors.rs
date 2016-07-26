
use bio::stats::{LogProb, logprobs};


/// A prior model of the allele frequency spectrum.
pub trait Model {
    /// Calculate prior probability of given allele frequency.
    fn prior_prob(&self, af: f64) -> LogProb;
}


pub trait DiscreteModel: Model {}


pub trait ContinuousModel: Model {}


/// The classical population genetic model used for variant calling in e.g. GATK and Samtools.
pub struct InfiniteSitesNeutralVariationModel {
    ploidy: u32,
    heterozygosity: f64,
    zero_prob: LogProb
}


impl InfiniteSitesNeutralVariationModel {
    /// Create new model for given ploidy and heterozygosity.
    pub fn new(ploidy: u32, heterozygosity: f64) -> Self {
        let zero_prob = logprobs::ln_1m_exp(
            heterozygosity.ln() +
            (0..ploidy * 2).fold(0.0, |s, m| s + 1.0 / m as f64).ln()
        );
        InfiniteSitesNeutralVariationModel {
            ploidy: ploidy,
            heterozygosity: heterozygosity,
            zero_prob: zero_prob
        }
    }
}


impl Model for InfiniteSitesNeutralVariationModel {
    fn prior_prob(&self, af: f64) -> LogProb {
        if af > 0.0 {
            // TODO fail for non-discrete m
            let m = af * self.ploidy as f64;
            self.heterozygosity.ln() - m.ln()
        } else {
            self.zero_prob
        }
    }
}


impl DiscreteModel for InfiniteSitesNeutralVariationModel {}


/// Flat priors.
pub struct FlatModel;


impl FlatModel {
    pub fn new() -> Self {
        FlatModel
    }
}


impl Model for FlatModel {
    fn prior_prob(&self, _: f64) -> LogProb {
        0.0
    }
}
