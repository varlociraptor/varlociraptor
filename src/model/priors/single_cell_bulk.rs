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
            allele_freqs_bulk: AlleleFreq(0.0)..AlleleFreq(1.00000001) // TODO: make use of inclusive ranges as soon as they arrive in stable rust
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
        // TODO: stick in meaningful prior for the bulk sample, derivative of InfiniteSitesNeutralVariationModel?
        // TODO: create a
        panic!("Priors are currently unsupported in the SingleCellBulkModel, feel free to implement them.")
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

        // stable rust ranges currently only allow inclusive start points, for exclusive add
        // 0.00000001 in the Event definition
        let k_start = (*af_bulk.start * n_obs_bulk as f64).ceil() as u64;
        // stable rust ranges currently only allow exclusive end points, for including an end point,
        // add 0.00000001 in the Event definition
        let k_end = (*af_bulk.end * n_obs_bulk as f64).ceil() as u64;
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
    use model::{AlleleFreq, AlleleFreqs, likelihood, Variant, PairPileup, priors};
    use model::sample::{Observation, Evidence};

    #[test]
    fn test_prob_rho() {
        let model = SingleCellBulkModel::new(2);
        // all expected results calculated with implementation in R version 3.3.3, using
        // dbetabinom.ab() from the R package VGAM_1.0-2

        // HOM REF model
        let results_5_hom_ref = [0.917057691405, 0.048991452054, 0.019369602690, 0.009095775973, 0.004098640948, 0.001386836930];
        let results_60_hom_ref = [  0.803320245461, 0.052192823008, 0.027012354914, 0.018046293954, 0.013409772796, 0.010564967885,
                                    0.008637021122, 0.007242360970, 0.006185855926, 0.005357571136, 0.004690762755, 0.004142520913,
                                    0.003683979192, 0.003294990190, 0.002961056264, 0.002671473498, 0.002418163441, 0.002194912251,
                                    0.001996860416, 0.001820151668, 0.001661685948, 0.001518942088, 0.001389848262, 0.001272685809,
                                    0.001166016809, 0.001068628819, 0.000979492206, 0.000897726836, 0.000822575822, 0.000753384624,
                                    0.000689584283, 0.000630677854, 0.000576229360, 0.000525854718, 0.000479214255, 0.000436006489,
                                    0.000395962936, 0.000358843750, 0.000324434041, 0.000292540765, 0.000262990063, 0.000235625005,
                                    0.000210303643, 0.000186897346, 0.000165289364, 0.000145373587, 0.000127053478, 0.000110241151,
                                    0.000094856580, 0.000080826922, 0.000068085955, 0.000056573603, 0.000046235582, 0.000037023138,
                                    0.000028892926, 0.000021807042, 0.000015733297, 0.000010645858, 0.000006526566, 0.000003367639, 0.000001177820];
        // HET model
        let results_5_het = [0.215810926125, 0.147529824097, 0.136659249777, 0.136659249777, 0.147529824097, 0.215810926125];
        let results_60_het = [  0.022964691724, 0.013791217110, 0.011203808536, 0.009976873288, 0.009344532591, 0.009066573603,
                                0.009039972286, 0.009210742960, 0.009545572404, 0.010020673323, 0.010616821444, 0.011316971031,
                                0.012105068132, 0.012965466790, 0.013882667547, 0.014841235701, 0.015825822418, 0.016821244988,
                                0.017812600331, 0.018785395739, 0.019725686644, 0.020620214709, 0.021456541698, 0.022223175971,
                                0.022909689396, 0.023506823069, 0.024006580679, 0.024402308658, 0.024688762485, 0.024862158676,
                                0.024920212140, 0.024862158676, 0.024688762485, 0.024402308658, 0.024006580679, 0.023506823069,
                                0.022909689396, 0.022223175971, 0.021456541698, 0.020620214709, 0.019725686644, 0.018785395739,
                                0.017812600331, 0.016821244988, 0.015825822418, 0.014841235701, 0.013882667547, 0.012965466790,
                                0.012105068132, 0.011316971031, 0.010616821444, 0.010020673323, 0.009545572404, 0.009210742960,
                                0.009039972286, 0.009066573603, 0.009344532591, 0.009976873288, 0.011203808536, 0.013791217110, 0.022964691724];
        // HOM ALT model
        let results_5_hom_alt = [0.001386836930, 0.004098640948, 0.009095775973, 0.019369602690, 0.048991452054, 0.917057691405];
        let results_60_hom_alt = [  0.000001177820, 0.000003367639, 0.000006526566, 0.000010645858, 0.000015733297, 0.000021807042,
                                    0.000028892926, 0.000037023138, 0.000046235582, 0.000056573603, 0.000068085955, 0.000080826922,
                                    0.000094856580, 0.000110241151, 0.000127053478, 0.000145373587, 0.000165289364, 0.000186897346,
                                    0.000210303643, 0.000235625005, 0.000262990063, 0.000292540765, 0.000324434041, 0.000358843750,
                                    0.000395962936, 0.000436006489, 0.000479214255, 0.000525854718, 0.000576229360, 0.000630677854,
                                    0.000689584283, 0.000753384624, 0.000822575822, 0.000897726836, 0.000979492206, 0.001068628819,
                                    0.001166016809, 0.001272685809, 0.001389848262, 0.001518942088, 0.001661685948, 0.001820151668,
                                    0.001996860416, 0.002194912251, 0.002418163441, 0.002671473498, 0.002961056264, 0.003294990190,
                                    0.003683979192, 0.004142520913, 0.004690762755, 0.005357571136, 0.006185855926, 0.007242360970,
                                    0.008637021122, 0.010564967885, 0.013409772796, 0.018046293954, 0.027012354914, 0.052192823008, 0.803320245461];
        // test all models
        for k in 0..5+1 {
            assert_relative_eq!( model.prob_rho(&AlleleFreq(0.0), &5, &(k as usize)).exp(),
                                    results_5_hom_ref[k] as f64, max_relative = 1.0,
                                    epsilon = 0.000000000001);
            assert_relative_eq!( model.prob_rho(&AlleleFreq(0.5), &5, &(k as usize)).exp(),
                                    results_5_het[k] as f64, max_relative = 1.0,
                                    epsilon = 0.000000000001);
            assert_relative_eq!( model.prob_rho(&AlleleFreq(1.0), &5, &(k as usize)).exp(),
                                    results_5_hom_alt[k] as f64, max_relative = 1.0,
                                    epsilon = 0.000000000001);
        }
        for k in 0..60+1 {
            assert_relative_eq!( model.prob_rho(&AlleleFreq(0.0), &60, &(k as usize)).exp(),
                                    results_60_hom_ref[k] as f64, max_relative = 1.0,
                                    epsilon = 0.000000000001);
            assert_relative_eq!( model.prob_rho(&AlleleFreq(0.5), &60, &(k as usize)).exp(),
                                    results_60_het[k] as f64, max_relative = 1.0,
                                    epsilon = 0.000000000001);
            assert_relative_eq!( model.prob_rho(&AlleleFreq(1.0), &60, &(k as usize)).exp(),
                                    results_60_hom_alt[k] as f64, max_relative = 1.0,
                                    epsilon = 0.000000000001);
        }
    }

    fn create_obs_vector(
        n_obs_ref: usize,
        n_obs_alt: usize
    ) -> Vec<Observation> {
        let obs_ref_abs = Observation {
            prob_mapping: LogProb::ln_one(),
            //prob_mapping: LogProb(0.9f64.ln()),
            prob_alt: LogProb::ln_zero(),
            prob_ref: LogProb::ln_one(),
            prob_mismapped: LogProb::ln_one(),
            evidence: Evidence::Alignment
        };
        let obs_alt_abs = Observation {
            prob_mapping: LogProb::ln_one(),
            //prob_mapping: LogProb(0.9f64.ln()),
            prob_alt: LogProb::ln_one(),
            prob_ref: LogProb::ln_zero(),
            prob_mismapped: LogProb::ln_one(),
            evidence: Evidence::Alignment
        };

        let mut obs = Vec::new();
        for _ in 0..n_obs_ref {
            obs.push(obs_ref_abs.clone());
        }
        for _ in 0..n_obs_alt {
            obs.push(obs_alt_abs.clone());
        }
        obs
    }

    fn create_test_snv_pileup<'a, A: AlleleFreqs, B: AlleleFreqs, P: priors::PairModel<A, B>>(
        prior_model: &'a P,
        s_obs: &Vec<Observation>,
        b_obs: &Vec<Observation>
    ) -> PairPileup<'a, A, B, P> {
        let variant = Variant::SNV(b'T');
        let single_sample_model = likelihood::LatentVariableModel::with_single_sample();
        let bulk_sample_model = likelihood::LatentVariableModel::with_single_sample();
        PairPileup::new(
            s_obs.clone(),
            b_obs.clone(),
            variant,
            prior_model,
            single_sample_model,
            bulk_sample_model,
        )
    }

    // tests that bulk range discretization inludes and excludes the right discrete values per
    // bulk frequency range (indirect through probabilities, as these bulk ranges aren't public)
    #[test]
    fn test_bulk_range_discretization() {
        let model = SingleCellBulkModel::new(2);

        let af_single_all = vec![AlleleFreq(0.0), AlleleFreq(0.5), AlleleFreq(1.0)];
        let af_bulk_zero = AlleleFreq(0.0)..AlleleFreq(0.00000001);
        let af_bulk_not_zero = AlleleFreq(0.00000001)..AlleleFreq(1.00000001);
        let af_bulk_quarter = AlleleFreq(0.25)..AlleleFreq(0.25000001);
        let af_bulk_half = AlleleFreq(0.5)..AlleleFreq(0.50000001);
        let af_bulk_not_one = AlleleFreq(0.0)..AlleleFreq(1.0);
        let af_bulk_one = AlleleFreq(1.0)..AlleleFreq(1.00000001);

        let s_dummy_obs = create_obs_vector(1, 1);

        let b_1_ref = create_obs_vector(1, 0);
        let mut pileup = create_test_snv_pileup(&model, &s_dummy_obs, &b_1_ref);
        println!("SCBM bulk ranges: pileup.marginal_prob(): {}", pileup.marginal_prob().exp() );
        println!("SCBM bulk ranges: pileup.joint_prob(af_single_all, af_bulk_zero): {}",
                    pileup.joint_prob(&af_single_all, &af_bulk_zero).exp() );
        println!("SCBM bulk ranges: pileup.joint_prob(af_single_all, af_bulk_not_zero): {}",
                    pileup.joint_prob(&af_single_all, &af_bulk_not_zero).exp() );
        assert_eq!(pileup.joint_prob(&af_single_all, &af_bulk_not_zero).exp(), 0.0);
        assert_eq!( (pileup.joint_prob(&af_single_all, &af_bulk_zero) - pileup.marginal_prob() ).exp(), 1.0);

        let b_1_alt = create_obs_vector(0, 1);
        pileup = create_test_snv_pileup(&model, &s_dummy_obs, &b_1_alt);
        println!("SCBM bulk ranges: pileup.marginal_prob(): {}", pileup.marginal_prob().exp() );
        println!("SCBM bulk ranges: pileup.joint_prob(af_single_all, af_bulk_one): {}",
                    pileup.joint_prob(&af_single_all, &af_bulk_one).exp() );
        println!("SCBM bulk ranges: pileup.joint_prob(af_single_all, af_bulk_not_one): {}",
                    pileup.joint_prob(&af_single_all, &af_bulk_not_one).exp() );
        assert_eq!(pileup.joint_prob(&af_single_all, &af_bulk_not_one).exp(), 0.0);
        assert_eq!( (pileup.joint_prob(&af_single_all, &af_bulk_one) - pileup.marginal_prob() ).exp(), 1.0);

        let b_2_ref = create_obs_vector(2, 0);
        pileup = create_test_snv_pileup(&model, &s_dummy_obs, &b_2_ref);
        println!("SCBM bulk ranges: pileup.marginal_prob(): {}", pileup.marginal_prob().exp() );
        println!("SCBM bulk ranges: pileup.joint_prob(af_single_all, af_bulk_zero): {}",
                    pileup.joint_prob(&af_single_all, &af_bulk_zero).exp() );
        println!("SCBM bulk ranges: pileup.joint_prob(af_single_all, af_bulk_half): {}",
                    pileup.joint_prob(&af_single_all, &af_bulk_half).exp() );
        println!("SCBM bulk ranges: pileup.joint_prob(af_single_all, af_bulk_one): {}",
                    pileup.joint_prob(&af_single_all, &af_bulk_one).exp() );
        assert_eq!( ( pileup.joint_prob(&af_single_all, &af_bulk_zero).ln_add_exp(
                        pileup.joint_prob(&af_single_all, &af_bulk_half)).ln_add_exp(
                            pileup.joint_prob(&af_single_all, &af_bulk_one)) -
                    pileup.marginal_prob() ).exp(), 1.0);
        println!("SCBM bulk ranges: pileup.joint_prob(af_single_all, af_bulk_quarter): {}",
                    pileup.joint_prob(&af_single_all, &af_bulk_quarter).exp() );
        assert_eq!(pileup.joint_prob(&af_single_all, &af_bulk_quarter).exp(), 0.0);

        let b_1_ref_1_alt = create_obs_vector(1, 1);
        pileup = create_test_snv_pileup(&model, &s_dummy_obs, &b_1_ref_1_alt);
        println!("SCBM bulk ranges: pileup.marginal_prob(): {}", pileup.marginal_prob().exp() );
        println!("SCBM bulk ranges: pileup.joint_prob(af_single_all, af_bulk_zero): {}",
                    pileup.joint_prob(&af_single_all, &af_bulk_zero).exp() );
        println!("SCBM bulk ranges: pileup.joint_prob(af_single_all, af_bulk_half): {}",
                    pileup.joint_prob(&af_single_all, &af_bulk_half).exp() );
        println!("SCBM bulk ranges: pileup.joint_prob(af_single_all, af_bulk_one): {}",
                    pileup.joint_prob(&af_single_all, &af_bulk_one).exp() );
        assert_eq!( ( pileup.joint_prob(&af_single_all, &af_bulk_zero).ln_add_exp(
                        pileup.joint_prob(&af_single_all, &af_bulk_half)).ln_add_exp(
                            pileup.joint_prob(&af_single_all, &af_bulk_one)) -
                    pileup.marginal_prob() ).exp(), 1.0);
        println!("SCBM bulk ranges: pileup.joint_prob(af_single_all, af_bulk_quarter): {}",
                    pileup.joint_prob(&af_single_all, &af_bulk_quarter).exp() );
        assert_eq!(pileup.joint_prob(&af_single_all, &af_bulk_quarter).exp(), 0.0);

        let b_2_alt = create_obs_vector(0, 2);
        pileup = create_test_snv_pileup(&model, &s_dummy_obs, &b_2_alt);
        println!("SCBM bulk ranges: pileup.marginal_prob(): {}", pileup.marginal_prob().exp() );
        println!("SCBM bulk ranges: pileup.joint_prob(af_single_all, af_bulk_zero): {}",
                    pileup.joint_prob(&af_single_all, &af_bulk_zero).exp() );
        println!("SCBM bulk ranges: pileup.joint_prob(af_single_all, af_bulk_half): {}",
                    pileup.joint_prob(&af_single_all, &af_bulk_half).exp() );
        println!("SCBM bulk ranges: pileup.joint_prob(af_single_all, af_bulk_one): {}",
                    pileup.joint_prob(&af_single_all, &af_bulk_one).exp() );
        assert_eq!( ( pileup.joint_prob(&af_single_all, &af_bulk_zero).ln_add_exp(
                        pileup.joint_prob(&af_single_all, &af_bulk_half)).ln_add_exp(
                            pileup.joint_prob(&af_single_all, &af_bulk_one)) -
                    pileup.marginal_prob() ).exp(), 1.0);
        println!("SCBM bulk ranges: pileup.joint_prob(af_single_all, af_bulk_quarter): {}",
                    pileup.joint_prob(&af_single_all, &af_bulk_quarter).exp() );
        assert_eq!(pileup.joint_prob(&af_single_all, &af_bulk_quarter).exp(), 0.0);

    }

    #[test]
    fn test_scbm_sc_het_germ() {

        let model = SingleCellBulkModel::new(2);

        // single cell is het against het germline
        let af_single = vec![AlleleFreq(0.5)];
        let af_bulk = AlleleFreq(0.25000001)..AlleleFreq(0.75000001);

        let single_obs = create_obs_vector(3, 3);
        let bulk_obs = create_obs_vector(3, 3);

        let pileup = create_test_snv_pileup(&model, &single_obs, &bulk_obs);

//        println!("SCBM.prior_prob(af_single[0] = {}, af_bulk.start = {}): {}", af_single[0], af_bulk.start, model.prior_prob(af_single[0], af_bulk.start, variant).exp() );
//        assert_eq!( model.prior_prob(af_single[0], af_bulk.start, variant).exp(), 1.0 );
        println!("SCBM het: pileup.joint_prob(af_single, af_bulk): {}", pileup.joint_prob(&af_single, &af_bulk).exp() );
        assert_relative_eq!( pileup.joint_prob(&af_single, &af_bulk).exp(), 0.00019473983947767667, epsilon = 0.000000000000000001 );
        println!("SCBM het: pileup.marginal_prob: {}", pileup.marginal_prob().exp() );
        assert_relative_eq!( pileup.marginal_prob().exp(), 0.00027403381880824275, epsilon = 0.000000000000000001 );
        println!("SCBM het: pileup.posterior_prob: {}", pileup.posterior_prob(&af_single, &af_bulk).exp() );
        assert_relative_eq!( pileup.posterior_prob(&af_single, &af_bulk).exp(), 0.710641629287177, epsilon = 0.000000000000001 );
        let (sc, blk) = pileup.map_allele_freqs();
        println!("SCBM het: pileup.map_allele_freqs: ({} {})", sc, blk );
        assert_eq!( pileup.map_allele_freqs(), ( AlleleFreq(0.5), AlleleFreq(0.5) ) );

        // test maximum a posteriori allele frequency estimates with uneven allele observations
        let single_obs = create_obs_vector(3, 7);
        let bulk_obs = create_obs_vector(3, 7);

        let pileup = create_test_snv_pileup(&model, &single_obs, &bulk_obs);
        let (sc, blk) = pileup.map_allele_freqs();
        println!("SCBM het: pileup.map_allele_freqs: ({} {})", sc, blk );
        assert_eq!( pileup.map_allele_freqs(), ( AlleleFreq(0.5), AlleleFreq(0.7) ) );
    }

    #[test]
    fn test_scbm_sc_hom_ref_germ() {

        let model = SingleCellBulkModel::new(2);

        // single cell is hom ref against hom ref germline
        let af_single = vec![AlleleFreq(0.0)];
        let af_bulk = AlleleFreq(0.0)..AlleleFreq(0.25000001);

        let single_obs = create_obs_vector(5, 0);
        let bulk_obs = create_obs_vector(5, 0);

        let pileup = create_test_snv_pileup(&model, &single_obs, &bulk_obs);

//        println!("SCBM.prior_prob(af_single[0] = {}, af_bulk.start = {}): {}", af_single[0], af_bulk.start, model.prior_prob(af_single[0], af_bulk.start, variant).exp() );
//        assert_eq!( model.prior_prob(af_single[0], af_bulk.start, variant).exp(), 1.0 );
        println!("SCBM hom ref: pileup.joint_prob(af_single, af_bulk): {}", pileup.joint_prob(&af_single, &af_bulk).exp() );
        assert_relative_eq!( pileup.joint_prob(&af_single, &af_bulk).exp(), 1.2409982197541034, epsilon = 0.0000000000000001 );
        println!("SCBM hom ref: pileup.marginal_prob: {}", pileup.marginal_prob().exp() );
        assert_relative_eq!( pileup.marginal_prob().exp(), 1.719859092473503, epsilon = 0.000000000000001 );
        println!("SCBM hom ref: pileup.posterior_prob: {}", pileup.posterior_prob(&af_single, &af_bulk).exp() );
        assert_relative_eq!( pileup.posterior_prob(&af_single, &af_bulk).exp(), 0.7215697060212639, epsilon = 0.0000000000000001 );
        let (sc, blk) = pileup.map_allele_freqs();
        println!("SCBM hom ref: pileup.map_allele_freqs: ({} {})", sc, blk );
        assert_eq!( pileup.map_allele_freqs(), ( AlleleFreq(0.0), AlleleFreq(0.0) ) );

        // test maximum a posteriori allele frequency estimates with alt reads
        let single_obs = create_obs_vector(57, 3);
        let bulk_obs = create_obs_vector(27, 3);

        let pileup = create_test_snv_pileup(&model, &single_obs, &bulk_obs);
        let (sc, blk) = pileup.map_allele_freqs();
        println!("SCBM hom ref: pileup.map_allele_freqs: ({} {})", sc, blk );
        assert_eq!( pileup.map_allele_freqs(), ( AlleleFreq(0.0), AlleleFreq(0.1) ) );
    }

    #[test]
    fn test_scbm_sc_hom_alt_germ() {

        let model = SingleCellBulkModel::new(2);

        // single cell is hom alt against hom alt germline
        let af_single = vec![AlleleFreq(1.0)];
        let af_bulk = AlleleFreq(0.75000001)..AlleleFreq(1.00000001);

        let single_obs = create_obs_vector(0, 5);
        let bulk_obs = create_obs_vector(0, 5);

        let pileup = create_test_snv_pileup(&model, &single_obs, &bulk_obs);

//        println!("SCBM.prior_prob(af_single[0] = {}, af_bulk.start = {}): {}", af_single[0], af_bulk.start, model.prior_prob(af_single[0], af_bulk.start, variant).exp() );
//        assert_eq!( model.prior_prob(af_single[0], af_bulk.start, variant).exp(), 1.0 );
        println!("SCBM hom alt: pileup.joint_prob(af_single, af_bulk): {}", pileup.joint_prob(&af_single, &af_bulk).exp() );
        assert_relative_eq!( pileup.joint_prob(&af_single, &af_bulk).exp(), 1.2409982197541034, epsilon = 0.0000000000000001 );
        println!("SCBM hom alt: pileup.marginal_prob: {}", pileup.marginal_prob().exp() );
        assert_relative_eq!( pileup.marginal_prob().exp(), 1.7198590924735035, epsilon = 0.0000000000000001 );
        println!("SCBM hom alt: pileup.posterior_prob: {}", pileup.posterior_prob(&af_single, &af_bulk).exp() );
        assert_relative_eq!( pileup.posterior_prob(&af_single, &af_bulk).exp(), 0.7215697060212638, epsilon = 0.0000000000000001 );
        let (sc, blk) = pileup.map_allele_freqs();
        println!("SCBM hom alt: pileup.map_allele_freqs: ({} {})", sc, blk );
        assert_eq!( pileup.map_allele_freqs(), ( AlleleFreq(1.0), AlleleFreq(1.0) ) );

        // test maximum a posteriori allele frequency estimates with alt reads
        let single_obs = create_obs_vector(3, 57);
        let bulk_obs = create_obs_vector(3, 27);

        let pileup = create_test_snv_pileup(&model, &single_obs, &bulk_obs);
        let (sc, blk) = pileup.map_allele_freqs();
        println!("SCBM hom alt: pileup.map_allele_freqs: ({} {})", sc, blk );
        assert_eq!( pileup.map_allele_freqs(), ( AlleleFreq(1.0), AlleleFreq(0.9) ) );
    }

    #[test]
    fn test_scbm_sc_hom_alt_som() {

        let model = SingleCellBulkModel::new(2);

        // single cell is hom alt against het germline
        let af_single = vec![AlleleFreq(1.0)];
        let af_bulk = AlleleFreq(0.25000001)..AlleleFreq(0.75000001);

        let single_obs = create_obs_vector(0, 4);
        let bulk_obs = create_obs_vector(3, 1);

        let pileup = create_test_snv_pileup(&model, &single_obs, &bulk_obs);

//        println!("SCBM.prior_prob(af_single[0] = {}, af_bulk.start = {}): {}", af_single[0], af_bulk.start, model.prior_prob(af_single[0], af_bulk.start, variant).exp() );
//        assert_eq!( model.prior_prob(af_single[0], af_bulk.start, variant).exp(), 1.0 );
        println!("SCBM hom alt: pileup.joint_prob(af_single, af_bulk): {}", pileup.joint_prob(&af_single, &af_bulk).exp() );
        assert_relative_eq!( pileup.joint_prob(&af_single, &af_bulk).exp(), 0.06996173990713889, epsilon = 0.0000000000000001 );
        println!("SCBM hom alt: pileup.marginal_prob: {}", pileup.marginal_prob().exp() );
        assert_relative_eq!( pileup.marginal_prob().exp(), 0.2267813569619291, epsilon = 0.0000000000000001 );
        println!("SCBM hom alt: pileup.posterior_prob: {}", pileup.posterior_prob(&af_single, &af_bulk).exp() );
        assert_relative_eq!( pileup.posterior_prob(&af_single, &af_bulk).exp(), 0.30849863870813554, epsilon = 0.00000000000000001 );
        let (sc, blk) = pileup.map_allele_freqs();
        println!("SCBM hom alt: pileup.map_allele_freqs: ({} {})", sc, blk );
        assert_eq!( pileup.map_allele_freqs(), ( AlleleFreq(1.0), AlleleFreq(0.25) ) );

        // test maximum a posteriori allele frequency estimates with alt reads
        let single_obs = create_obs_vector(4, 56);
        let bulk_obs = create_obs_vector(21, 9);

        let pileup = create_test_snv_pileup(&model, &single_obs, &bulk_obs);
        let (sc, blk) = pileup.map_allele_freqs();
        println!("SCBM hom alt: pileup.map_allele_freqs: ({} {})", sc, blk );
        assert_eq!( pileup.map_allele_freqs(), ( AlleleFreq(1.0), AlleleFreq(0.3) ) );
    }
}
