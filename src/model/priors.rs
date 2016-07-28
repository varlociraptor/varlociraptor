
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


/// Model of neutral mutation in tumor cell populations described in
/// Williams et al. Nature Genetics 2016.
pub struct WilliamsTumorModel {
    effective_mutation_rate: f64,
    average_ploidy: u32,
    genome_size: u64
}


impl WilliamsTumorModel {
    /// Create new model for given average ploidy and mutation rate per effective cell division.
    /// The latter, is the quotient mu/beta, with mu being the mutation rate and beta being the fraction
    /// of effective cell divisions (both lineages survive). Alone, the parameters are not observable.
    /// However, mu/beta can be estimated from e.g. SNV calls. It is the slope of the linear model
    /// y = mu/beta * (x -  1 / fmax), with x being the reciprocal of the observed allele frequencies
    /// and y being the number of observed mutations corresponding to each frequency
    /// (see Williams et al. Nature Genetics 2016).
    ///
    /// Based on the model, the prior probability of an allele frequency F = f is
    /// Pr(F = f) = M(f) / n = mu/beta (1 / f - 1 / fmax) / n
    /// with n being the size of the genome and fmax is the expected allele frequency of clonal variants
    /// at the beginning of tumor evolution. Hence, fmax is set to fmax = 1 / ploidy.
    ///
    /// # Arguments
    ///
    /// * `average_ploidy` - the average ploidy of the tumor (e.g. 2 for diploid)
    /// * `effective_mutation_rate` - the mutation rate per effective cell division
    pub fn new(average_ploidy: u32, effective_mutation_rate: f64, genome_size: u64) -> Self {
        WilliamsTumorModel {
            average_ploidy: average_ploidy,
            effective_mutation_rate: effective_mutation_rate,
            genome_size: genome_size
        }
    }
}


impl Model for WilliamsTumorModel {
    fn prior_prob(&self, af: f64) -> LogProb {
        let x = 1.0 / af - self.average_ploidy as f64;
        self.effective_mutation_rate.ln() + x.ln() - (self.genome_size as f64).ln()
    }
}


impl ContinuousModel for WilliamsTumorModel {}


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
