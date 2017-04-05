use itertools::{Itertools};
use ordered_float::NotNaN;
use bio::stats::LogProb;
use statrs::function::beta::ln_beta;
use statrs::function::factorial::ln_binomial;

use model::{Variant, ContinuousAlleleFreqs, DiscreteAlleleFreqs, AlleleFreq};

use priors::PairModel;

/// Prior model for a Single Cell against a Bulk background from the same individual (optimally the
/// same cell type). It uses the ploidy of the organism as well as a WGA method specific single cell
/// model (so far only for MDA) to account for differential allelic amplification.
/// TODO: * use the general level of heterozygosity through the InfiniteSitesNeutralEvolutionModel as
///         a prior? "The prior probability for a germline allele frequency theta_g (e.g. 0.0, 0.5 or 1.0 for the diploid case) in the bulk background can be calculated with an `InfiniteSitesNeutralVariationModel`. This is valid since clonal variants come from the last common ancestor and analogously to tumor evolution in the Williams model, we can assume neutral mutations (no genetic drift, no selection) and thus no change of allele frequencies in cell divisions that do not introduce new mutations. The `InfiniteSitesNeutralVariationModel` requires the ploidy and the level of heterozygosity."
///       * use the somatic mutation rate per effective cell division? ("The somatic mutation rate per effective cell division in the bulk is the quotient mu/beta, with mu being the somatic mutation rate and beta being the fraction of effective cell divisions (i.e. where both daugther cells survive and form a lineage). Alone, these parameters are not easily obtained. However, assuming mostly neutral mutations, mu/beta can be estimated from SNV calls with a low frequency in the bulk sample, analogous to the tumour sample in Williams et al. (2016). It is the slope of the linear model `y = mu/beta * (x -  1 / fmax)`, with `x` being the reciprocal of the observed allele frequencies and y being the number of observed mutations corresponding to each frequency (see: Williams MJ, Werner B, Barnes CP, Graham TA, Sottoriva A. Identification of neutral tumor evolution across cancer types. Nat Genet. 2016;48: 238–244. doi:10.1038/ng.3489). Based on the Williams model, the tail probability of a somatic allele frequency F > f can be expressed as `Pr(F > f) = M(f) / n = mu/beta (1 / f - 1 / fmax) / n`, with `n` being the size of the genome and `fmax` the expected allele frequency of clonal variants at the beginning of tumor history, overall somatic history in our case. From this, we can obtain the cumulative distribution function as `Pr(F <= f) = 1 - Pr(F > f)`. Consequently, the density becomes the first derivative, i.e. `Pr(F = f) = - M(f)' / n = mu/beta * 1/n * 1/f²` for `f>=fmin`, with `fmin = sqrt(mu/beta * 1/n)`."
pub struct SingleCellBulkModel {
    allele_freqs_single: DiscreteAlleleFreqs,
    allele_freqs_bulk: ContinuousAlleleFreqs
}

impl SingleCellBulkModel {
    /// Create new model.
    ///
    /// # Arguments
    ///
    /// * `ploidy` - the ploidy in the single cell sample (e.g. 2 for diploid)
    pub fn new(ploidy: u32) -> Self {
        let allele_freqs = (0..ploidy + 1).map(|m| AlleleFreq(m as f64 / ploidy as f64)).collect_vec();
        SingleCellBulkModel {
            allele_freqs_single: allele_freqs,
            allele_freqs_bulk: AlleleFreq(0.0)..AlleleFreq(1.0),
        }
    }

    // Lodato et al. 2015, Science, Supplementary Information, pages 8f and Fig. S5 (A, C, E)
    // TODO: allow for non-default Lodato model parameters, e.g. learned from the data at hand
    // TODO: allow for non-Lodato models
    fn prob_rho(
        &self,
        af_single_underlying: &f64,
        n_obs: &usize,
        k: &usize
    ) -> LogProb
    {
        let binomial_coeff = ln_binomial(*n_obs as u64,*k as u64);

        match *af_single_underlying {
            // model for hom ref sites
            0.0 => {
                let alpha = |cov| {
                    -0.000027183 * cov as f64 + 0.068567471
                };
                let beta = |cov| {
                    0.007454388 * cov as f64 + 2.367486659
                };
                let a = alpha(*n_obs);
                let b = beta(*n_obs);

                LogProb(
                    binomial_coeff +
                        ln_beta((*k as f64) + a,  (*n_obs) as f64 - (*k as f64) + b) - ln_beta(a,b)
                )
            },
            // model for heterozygous sites
            0.5 => {
                let weight = |cov| {
                    0.000548761 * cov as f64 + 0.540396786
                };
                let alpha_1 = |cov| {
                    0.057378844 * cov as f64 + 0.669733191
                };
                let alpha_2 = |cov| {
                    0.003233912 * cov as f64 + 0.399261625
                };

                let a1 = alpha_1(*n_obs);
                let a2 = alpha_2(*n_obs);
                let w = weight(*n_obs);

                // Lodato et al. 2015, Science, Supplementary Information, pages 8f and Fig. S5 (A, C, E)
                LogProb(
                    w.ln() + binomial_coeff +
                        ln_beta((*k as f64) + a1,  (*n_obs) as f64 - (*k as f64) + a1) - ln_beta(a1,a1)
                ).ln_add_exp(
                    LogProb(
                        ((-w).ln_1p()) + binomial_coeff +
                            ln_beta((*k as f64) + a2, (*n_obs) as f64 - (*k as f64) + a2) - ln_beta(a2,a2) )
                )
            },
            // model for hom alt sites (hom ref density mirrored)
            1.0 => {
                let alpha = |cov| {
                    0.007454388 * cov as f64 + 2.367486659
                };
                let beta = |cov| {
                    -0.000027183 * cov as f64 + 0.068567471
                };
                let a = alpha(*n_obs);
                let b = beta(*n_obs);

                LogProb(
                    binomial_coeff +
                        ln_beta((*k as f64) + a,  (*n_obs) as f64 - (*k as f64) + b) - ln_beta(a,b)
                )
            },
            _ => panic!("SingleCellBulkModel is currently only implemented for the diploid case with allele frequencies 0.0, 0.5 and 1.0.")
        }
    }
}


impl PairModel<DiscreteAlleleFreqs, ContinuousAlleleFreqs> for SingleCellBulkModel {

    fn prior_prob(&self, _: AlleleFreq, _: AlleleFreq, _: Variant) -> LogProb {
        // TODO: stick in the InfiniteSitesNeutralVariationModel here?
        LogProb::ln_one()
    }

    fn joint_prob<L, O>(
        &self,
        af_single: &DiscreteAlleleFreqs,
        af_bulk: &ContinuousAlleleFreqs,
        likelihood_single_distorted: &L,
        likelihood_bulk: &O,
        _: Variant,
        n_obs_single: usize,
        n_obs_bulk: usize
    ) -> LogProb where
        L: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb,
        O: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb
    {
        // cap the use of the single cell amplification bias model at a coverage of 100,
        // as the Lodato et al. model was fit with coverages capped at 60 and starts behaving
        // weirdly above 100
        // TODO: make this optional and dependent on the usage of the Lodato model with their params
        let n_obs_s = if n_obs_single > 100 {
            100
        } else {
            n_obs_single
        };
        let k_single = 0..n_obs_s + 1;

        let k_start = if *af_bulk.start == 0.0 { // 0 as a bulk range start is always inclusive
            0 as u64
        } else {
            // any non-zero bulk range start is always exclusive
            (*af_bulk.start * n_obs_bulk as f64).floor() as u64 + 1
        };
        // any bulk range end is always inclusive
        let k_end = (*af_bulk.end * n_obs_bulk as f64).floor() as u64 + 1;
        let k_bulk = k_start..k_end;

        // sum up all possible discrete bulk allele frequencies with current number of observations
        let p_bulk = LogProb::ln_sum_exp(&k_bulk.map(|k_b| {
            let af_bulk = AlleleFreq(k_b as f64/n_obs_bulk as f64);
            likelihood_bulk(af_bulk, None)
        }).collect_vec() );

        // go through all possible underlying single cell allele frequencies
        let prob = LogProb::ln_sum_exp(&af_single.iter().map(|&af_single| {
            let p_single =
                    LogProb::ln_sum_exp(&k_single.clone().map(|k_s| { // sum up all possible discrete single cell allele frequencies with current number of observations
                        let af_single_distorted = AlleleFreq(k_s as f64/n_obs_s as f64);
                        likelihood_single_distorted(af_single_distorted, None) + self.prob_rho(&af_single, &n_obs_s, &k_s)
                    }).collect_vec());
            let prob = p_bulk + p_single;

            prob
        }).collect_vec());

        prob
    }

    fn marginal_prob<L, O>(
        &self,
        likelihood_single_distorted: &L,
        likelihood_bulk: &O,
        variant: Variant,
        n_obs_single: usize,
        n_obs_bulk: usize
    ) -> LogProb where
        L: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb,
        O: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb
    {
        self.joint_prob(
            self.allele_freqs().0,
            self.allele_freqs().1,
            likelihood_single_distorted,
            likelihood_bulk,
            variant,
            n_obs_single,
            n_obs_bulk
        )
    }

    fn map<L, O>(
        &self,
        likelihood_single_distorted: &L,
        likelihood_bulk: &O,
        _: Variant,
        n_obs_single: usize,
        n_obs_bulk: usize
    ) -> (AlleleFreq, AlleleFreq) where
        L: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb,
        O: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb
    {
        // cap the use of the single cell amplification bias model at a coverage of 100,
        // as the Lodato et al. model was fit with coverages capped at 60 and starts behaving
        // weirdly above 100
        let n_obs_s = if n_obs_single > 100 {
            100
        } else {
            n_obs_single
        };
        let k_single = 0..n_obs_s + 1;
        let (_, map_single) = self.allele_freqs().0.iter().minmax_by_key(
            |&&af_single| {
                let p_single =
                    LogProb::ln_sum_exp(&k_single.clone().map(|k_s| { // sum up all possible discrete single cell allele frequencies with current number of observations
                        let af_single_distorted = AlleleFreq(k_s as f64/n_obs_s as f64);
                        likelihood_single_distorted(af_single_distorted, None) + self.prob_rho(&af_single, &n_obs_s, &k_s)
                    }).collect_vec());
                NotNaN::new(*p_single).expect("posterior probability is NaN")
            }
        ).into_option().expect("prior has empty allele frequency spectrum");

        let k_bulk = 0..n_obs_bulk + 1;
        let afs_bulk = k_bulk.map(|k| AlleleFreq(k as f64/n_obs_bulk as f64));
        let (_, map_bulk) = afs_bulk.minmax_by_key(
            |af_bulk| {
                let p_bulk = likelihood_bulk(*af_bulk, None);
                NotNaN::new(*p_bulk).expect("posterior probability is NaN")
            }
        ).into_option().expect("prior has empty allele frequency spectrum");

        (*map_single, map_bulk)
    }

    fn allele_freqs(&self) -> (&DiscreteAlleleFreqs, &ContinuousAlleleFreqs) {
        (&self.allele_freqs_single, &self.allele_freqs_bulk)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::stats::LogProb;
    use model::priors::PairModel;
    use model::AlleleFreq;

    #[test]
    fn test_prob_rho() {
        let model = SingleCellBulkModel::new(2);
        // expected results calculated with implementation in R version 3.3.2, using dbetabinom.ab()
        // from the R package VGAM_1.0-2
        let results_5 = [0.2158109, 0.1475298, 0.1366592, 0.1366592, 0.1475298, 0.2158109];
        for k in 0..5+1 {
            assert_relative_eq!( model.prob_rho(&5,&(k as usize)).exp(), results_5[k] as f64,
                                    max_relative = 1.0, epsilon = 0.0000001 );
        }
        // expected results calculated with implementation in R version 3.3.2, using dbetabinom.ab()
        // from the R package VGAM_1.0-2
        let results_60 = [  0.022964692, 0.013791217, 0.011203809, 0.009976873, 0.009344533,
                            0.009066574, 0.009039972, 0.009210743, 0.009545572, 0.010020673,
                            0.010616821, 0.011316971, 0.012105068, 0.012965467, 0.013882668,
                            0.014841236, 0.015825822, 0.016821245, 0.017812600, 0.018785396,
                            0.019725687, 0.020620215, 0.021456542, 0.022223176, 0.022909689,
                            0.023506823, 0.024006581, 0.024402309, 0.024688762, 0.024862159,
                            0.024920212, 0.024862159, 0.024688762, 0.024402309, 0.024006581,
                            0.023506823, 0.022909689, 0.022223176, 0.021456542, 0.020620215,
                            0.019725687, 0.018785396, 0.017812600, 0.016821245, 0.015825822,
                            0.014841236, 0.013882668, 0.012965467, 0.012105068, 0.011316971,
                            0.010616821, 0.010020673, 0.009545572, 0.009210743, 0.009039972,
                            0.009066574, 0.009344533, 0.009976873, 0.011203809, 0.013791217, 0.022964692];
        for k in 0..60+1 {
            assert_relative_eq!( model.prob_rho(&60,&(k as usize)).exp(), results_60[k] as f64,
                                    max_relative = 1.0, epsilon = 0.0000001 );
        }
    }
}
