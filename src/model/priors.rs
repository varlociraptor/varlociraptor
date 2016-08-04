use std::f64;
use std::ops::Range;
use std::cmp;

use itertools::Itertools;
use ordered_float::NotNaN;
use bio::stats::{LogProb, logprobs};


pub type DiscreteAlleleFreq = Vec<f64>;
pub type ContinousAlleleFreq = Range<f64>;


pub trait AlleleFreq {}


impl AlleleFreq for DiscreteAlleleFreq {}


impl AlleleFreq for ContinousAlleleFreq {}


/// A prior model of the allele frequency spectrum.
pub trait Model<A: AlleleFreq> {
    fn prior_prob(&self, af: f64) -> LogProb;

    /// Calculate prior probability of given allele frequency.
    fn joint_prob<E: Fn(f64) -> LogProb>(&self, af: &A, event_prob: &E) -> LogProb;

    fn allele_freqs(&self) -> &A;
}


/// The classical population genetic model used for variant calling in e.g. GATK and Samtools.
pub struct InfiniteSitesNeutralVariationModel {
    ploidy: u32,
    heterozygosity: f64,
    zero_prob: LogProb,
    allele_freqs: Vec<f64>
}


impl InfiniteSitesNeutralVariationModel {
    /// Create new model for given ploidy and heterozygosity.
    pub fn new(ploidy: u32, heterozygosity: f64) -> Self {
        let zero_prob = logprobs::ln_1m_exp(
            heterozygosity.ln() +
            (1..ploidy + 1).fold(0.0, |s, m| s + 1.0 / m as f64).ln()
        );

        let allele_freqs = (0..ploidy + 1).map(|m| m as f64 / ploidy as f64).collect_vec();

        InfiniteSitesNeutralVariationModel {
            ploidy: ploidy,
            heterozygosity: heterozygosity,
            zero_prob: zero_prob,
            allele_freqs: allele_freqs
        }
    }
}


impl Model<DiscreteAlleleFreq> for InfiniteSitesNeutralVariationModel {
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

    fn allele_freqs(&self) -> &DiscreteAlleleFreq {
        &self.allele_freqs
    }

    fn joint_prob<E: Fn(f64) -> LogProb>(&self, af: &DiscreteAlleleFreq, event_prob: &E) -> LogProb {
        logprobs::sum(&af.iter().map(|&af| self.prior_prob(af) + event_prob(af)).collect_vec())
    }
}


/// Mixture of tumor and normal priors.
/// The tumor model uses a published model of neutral mutation in tumor cell populations
/// described in Williams et al. Nature Genetics 2016.
pub struct TumorModel {
    effective_mutation_rate: f64,
    genome_size: u64,
    purity: f64,
    normal_model: InfiniteSitesNeutralVariationModel,
    grid_points: usize,
    panels: Vec<NotNaN<f64>>,
    allele_freqs: ContinousAlleleFreq
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
    /// Based on the Williams model, the tail probability of a somatic allele frequency F > f can be expressed
    /// as
    /// `Pr(F > f) = M(f) / n = mu/beta (1 / f - 1 / fmax) / n`
    /// with `n` being the size of the genome and `fmax` is the expected allele frequency of clonal variants
    /// at the beginning of tumor evolution.
    /// From this, we can obtain the cumulative distribution function as `Pr(F <= f) = 1 - Pr(F > f)`.
    /// Consequently, the density becomes the first derivate, i.e. `Pr(F = f) = - M(f)' / n = mu/beta * 1 / (f²n)`.
    ///
    /// The prior probability for a germline allele frequency f (e.g. 0.0, 0.5, 1.0) in the tumor is
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
        assert!(purity <= 1.0 && purity >= 0.0);
        let normal_model = InfiniteSitesNeutralVariationModel::new(ploidy, heterozygosity);

        // calculate outer panels for integration
        // we choose these such that we always start at a new peak in the density (e.g. at 0.0, 0.5)
        let mut panels = Vec::new();
        panels.extend(normal_model.allele_freqs().iter().map(|&af| NotNaN::new(af).unwrap()));
        if purity < 1.0 {
            // handle additional peaks
            panels.extend(normal_model.allele_freqs().iter().filter_map(|af| {
                let f = af / (1.0 - purity);
                if f > 0.0 && f <= 1.0 {
                    Some(NotNaN::new(f).unwrap())
                } else {
                    None
                }
            }));
        }
        panels = panels.into_iter().unique().collect_vec();

        let model = TumorModel {
            effective_mutation_rate: effective_mutation_rate,
            genome_size: genome_size,
            purity: purity,
            normal_model: normal_model,
            grid_points: 200,
            panels: panels,
            allele_freqs: 0.0..1.0
        };

        model
    }

    fn somatic_density(&self, af: f64) -> LogProb {
        // mu/beta * 1 / (af**2 * n)
        if af == 0.0 {
            // undefined for 0, but returning 0 for convenience
            return 0.0f64.ln()
        }
        self.effective_mutation_rate.ln() - (2.0 * af.ln() + (self.genome_size as f64).ln())
    }

    fn germline_density(&self, af: f64) -> LogProb {
        self.normal_model.prior_prob(af)
    }

    fn tumor_density(&self, af: f64) -> LogProb {
        let probs = self.normal_model.allele_freqs().iter().filter_map(|&af_germline| {
            if af >= af_germline {
                let af_somatic = af - af_germline;

                let p_germline = self.germline_density(af_germline);

                let p_somatic = if af_somatic != 0.0 {
                    self.somatic_density(af_somatic)
                } else {
                    // do not consider somatic density for allele freq = 0
                    0.0
                };
                Some(p_germline + p_somatic)
            } else {
                None
            }
        }).collect_vec();

        // summands are disjoint because we go over disjoint germline events
        logprobs::sum(&probs)
    }
}


impl Model<ContinousAlleleFreq> for TumorModel {
    fn prior_prob(&self, af: f64) -> LogProb {
        // sum over the different possibilities to obtain af
        let probs = self.normal_model.allele_freqs().iter().filter_map(|&af_normal| {
            let f = (1.0 - self.purity) * af_normal;
            if af >= f {
                // af = purity * af_tumor + (1-purity) * af_normal
                let af_tumor = (af - f) / self.purity;

                let p_tumor = self.tumor_density(af_tumor);
                let p_normal = self.normal_model.prior_prob(af_normal);
                Some(p_tumor + p_normal)
            } else {
                None
            }
        }).collect_vec();

        // summands are disjoint because we go over disjoint events on the normal sample
        let p = logprobs::sum(&probs);
        p
    }

    fn joint_prob<E: Fn(f64) -> LogProb>(&self, af: &ContinousAlleleFreq, event_prob: &E) -> LogProb {
        let density = |af| self.prior_prob(af) + event_prob(af);
        let mut summands = Vec::with_capacity(self.panels.len());
        let af_min = NotNaN::new(af.start).expect("Allele frequency range may not be NaN.");
        let af_max = NotNaN::new(af.end).expect("Allele frequency range may not be NaN.");

        for i in 0..self.panels.len() - 1 {
            let (fmin, fmax) = (self.panels[i], self.panels[i + 1]);
            if fmin >= af_max {
                break;
            }
            let fmax = cmp::min(fmax, af_max);
            if fmin >= af_min {
                // add density(fmin) because it contains the discrete peak
                summands.push(density(*fmin));
            }
            let fmin = cmp::max(fmin, af_min);
            summands.push(logprobs::integrate(&density, *fmin, *fmax, self.grid_points));
        }
        logprobs::sum(&summands)
    }

    fn allele_freqs(&self) -> &ContinousAlleleFreq {
        &self.allele_freqs
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use itertools::linspace;

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
        for purity in linspace(0.5, 1.0, 5) {
            println!("purity {}", purity);
            let model = TumorModel::new(2, 300.0, 3e9 as u64, purity, 0.001);
            println!("af=0.0 -> {}", model.prior_prob(0.0));
            let total = model.integrate(&(0.0..1.0), &|_| 0.0);
            println!("total {}", total);

            for af in linspace(0.0, 1.0, 20) {
                println!("af={} p={}", af, model.prior_prob(af));
            }
            assert_relative_eq!(total.exp(), 1.0, epsilon=0.01);
        }
    }
}
