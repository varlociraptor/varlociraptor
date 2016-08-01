use itertools::Itertools;

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
            (1..ploidy + 1).fold(0.0, |s, m| s + 1.0 / m as f64).ln()
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
            let m = af * self.ploidy as f64;
            if relative_eq!(m % 1.0, 0.0) {
                // if m is discrete
                self.heterozygosity.ln() - m.ln()
            } else {
                // invalid allele frequency
                0.0f64.ln()
            }
        } else {
            self.zero_prob
        }
    }
}


impl DiscreteModel for InfiniteSitesNeutralVariationModel {}


/// Mixture of tumor and normal priors.
/// The tumor model uses a published model of neutral mutation in tumor cell populations
/// described in Williams et al. Nature Genetics 2016.
pub struct TumorModel {
    effective_mutation_rate: f64,
    genome_size: u64,
    purity: f64,
    normal_model: InfiniteSitesNeutralVariationModel
}


impl TumorModel {
    /// Create new model for given ploidy, heterozygosity (in normal tissue) and tumor mutation rate
    /// per effective cell division.
    /// The latter is the quotient mu/beta, with mu being the mutation rate and beta being the fraction
    /// of effective cell divisions (both lineages survive). Alone, the parameters are not observable.
    /// However, mu/beta can be estimated from e.g. SNV calls. It is the slope of the linear model
    /// `y = mu/beta * (x -  1 / fmax)``, with `x` being the reciprocal of the observed allele frequencies
    /// and y being the number of observed mutations corresponding to each frequency
    /// (see Williams et al. Nature Genetics 2016).
    ///
    /// Based on the Williams model, the prior probability of a subclonal allele frequency F = f is
    /// `p_sub = M(f) / n = mu/beta (1 / f - 1 / fmax) / n`
    /// with `n` being the size of the genome and `fmax` is the expected allele frequency of clonal variants
    /// at the beginning of tumor evolution. The model assumes that the tumor comes from a single cell.
    /// Hence, `fmax` is set to `fmax = 1 / ploidy`.
    ///
    /// The prior probability for a clonal allele frequency f (e.g. 0.0, 0.5, 1.0) in the tumor is
    /// calculated with an `InfiniteSitesNeutralVariationModel`. This is valid since clonal variants
    /// come from the underlying normal tissue and Williams model assumes that allele frequencies
    /// do not change during tumor evolution (no genetic drift, no selection).
    ///
    /// For the final prior, we consider a given tumor purity and calculate the combined prior
    /// for all possible allele frequency combinations satisfying `af = purity * af_tumor + (1-purity) * af_normal`.
    ///
    /// # Arguments
    ///
    /// * `ploidy` - the ploidy in the corresponding normal sample (e.g. 2 for diploid)
    /// * `effective_mutation_rate` - the mutation rate per effective cell division in the tumor
    /// * `genome_size` - the size of the genome
    /// * `purity` - tumor purity
    /// * `heterozygosity` - expected heterozygosity in the corresponding normal
    pub fn new(
        ploidy: u32,
        effective_mutation_rate: f64,
        genome_size: u64,
        purity: f64,
        heterozygosity: f64) -> Self {
        let normal_model = InfiniteSitesNeutralVariationModel::new(ploidy, heterozygosity);
        TumorModel {
            effective_mutation_rate: effective_mutation_rate,
            genome_size: genome_size,
            purity: purity,
            normal_model: normal_model
        }
    }

    fn prior_prob_tumor(&self, af_tumor: f64) -> LogProb {
        if af_tumor > 0.0 && af_tumor < 1.0 / self.normal_model.ploidy as f64 {
            // subclonal event
            let x = 1.0 / af_tumor - self.normal_model.ploidy as f64;
            (self.effective_mutation_rate * x).ln() - (self.genome_size as f64).ln()
        } else {
            // clonal event
            self.normal_model.prior_prob(af_tumor)
        }
    }
}


impl Model for TumorModel {
    fn prior_prob(&self, af: f64) -> LogProb {
        // af = purity * af_tumor + (1-purity) * af_normal

        // sum over the different possibilities to obtain af
        let probs = (0..self.normal_model.ploidy + 1).filter_map(|m| {
            let af_normal = m as f64 / self.normal_model.ploidy as f64;
            if af >= af_normal {
                let af_tumor = (af - (1.0 - self.purity) * af_normal) / self.purity;

                let p_tumor = self.prior_prob_tumor(af_tumor);
                let p_normal = self.normal_model.prior_prob(af_normal);
                Some(p_tumor + p_normal)
            } else {
                None
            }
        }).collect_vec();

        logprobs::sum(&probs) - (probs.len() as f64).ln()
    }
}


impl ContinuousModel for TumorModel {}


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


#[cfg(test)]
mod tests {
    use super::*;
    use itertools::linspace;
    use bio::stats::logprobs;

    #[test]
    fn test_flat() {
        let model = FlatModel::new();
        assert_relative_eq!(model.prior_prob(0.1), 1.0f64.ln());
    }

    #[test]
    fn test_infinite_sites_neutral_variation() {
        let ploidy = 2;
        let het = 0.001;
        let model = InfiniteSitesNeutralVariationModel::new(ploidy, het);
        assert_relative_eq!(model.prior_prob(0.5).exp(), 0.001);
        assert_relative_eq!(model.prior_prob(1.0).exp(), 0.0005);
        assert_relative_eq!(model.prior_prob(0.0).exp(), 0.9985);
    }

    #[test]
    fn test_tumor() {
        //let mut model = TumorModel::new(2, 30.0, 3e9 as u64, 0.5, 0.001);
        //model.prior_prob(0.0);
        //assert!(false);
        for purity in linspace(0.5, 1.0, 5) {
            let model = TumorModel::new(2, 30.0, 3e9 as u64, purity, 0.001);
            let density = |af| model.prior_prob(af);
            let total = logprobs::integrate(&density, 0.0, 1.0, 2000);
            println!("purity={}", purity);
            for af in linspace(0.0, 1.0, 10) {
                println!("af={}, p={}", af, density(af).exp());
            }

            println!("total subclonal {}", logprobs::integrate(&density, 0.0, 0.5, 2000).exp());
            assert_relative_eq!(total.exp() + density(0.0).exp(), 1.0, epsilon=0.01);
        }
    }
}
