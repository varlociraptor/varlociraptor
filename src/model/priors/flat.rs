use itertools::Itertools;
use itertools_num::linspace;
use ordered_float::NotNaN;
use bio::stats::LogProb;

use model::{Variant, ContinuousAlleleFreqs, DiscreteAlleleFreqs, AlleleFreq};

use priors::PairModel;


pub struct FlatNormalNormalModel {
    allele_freqs: DiscreteAlleleFreqs
}


impl FlatNormalNormalModel {
    pub fn new(ploidy: u32) -> Self {
        let allele_freqs = (0..ploidy + 1).map(|m| AlleleFreq(m as f64 / ploidy as f64)).collect_vec();
        FlatNormalNormalModel {
            allele_freqs: allele_freqs
        }
    }
}


impl PairModel<DiscreteAlleleFreqs, DiscreteAlleleFreqs> for FlatNormalNormalModel {

    fn prior_prob(&self, _: AlleleFreq, _: AlleleFreq, _: Variant) -> LogProb {
        LogProb::ln_one()
    }

    fn joint_prob<L, O>(
        &self,
        af_first: &DiscreteAlleleFreqs,
        af_second: &DiscreteAlleleFreqs,
        likelihood_first: &L,
        likelihood_second: &O,
        _: Variant,
        _: usize,
        _: usize
    ) -> LogProb where
        L: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb,
        O: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb
    {
        let p_second = LogProb::ln_sum_exp(&af_second.iter().map(|&af_second| {
            likelihood_second(af_second, None)
        }).collect_vec());

        let prob = LogProb::ln_sum_exp(&af_first.iter().map(|&af_first| {
            let p_first = likelihood_first(af_first, None);
            let prob = p_first + p_second;

            prob
        }).collect_vec());

        prob
    }

    fn marginal_prob<L, O>(
        &self,
        likelihood_first: &L,
        likelihood_second: &O,
        variant: Variant,
        n_obs_first: usize,
        n_obs_second: usize
    ) -> LogProb where
        L: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb,
        O: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb
    {
        let p = self.joint_prob(
            self.allele_freqs().0,
            self.allele_freqs().1,
            likelihood_first,
            likelihood_second,
            variant,
            n_obs_first,
            n_obs_second
        );

        p
    }

    fn map<L, O>(
        &self,
        likelihood_first: &L,
        likelihood_second: &O,
        _: Variant,
        _: usize,
        _: usize
    ) -> (AlleleFreq, AlleleFreq) where
        L: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb,
        O: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb
    {
        fn calc_map<L: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb>(likelihood: &L, afs: &DiscreteAlleleFreqs) -> AlleleFreq {
            let (_, map) = afs.iter().minmax_by_key(|&af| {
                let p = likelihood(*af, None);
                NotNaN::new(*p).expect("probability is NaN")
            }).into_option().expect("prior has empty allele frequency spectrum");
            *map
        }

        let map_first = calc_map(likelihood_first, self.allele_freqs().0);
        let map_second = calc_map(likelihood_second, self.allele_freqs().1);

        (map_first, map_second)
    }

    fn allele_freqs(&self) -> (&DiscreteAlleleFreqs, &DiscreteAlleleFreqs) {
        (&self.allele_freqs, &self.allele_freqs)
    }
}


pub struct FlatTumorNormalModel {
    allele_freqs_tumor: ContinuousAlleleFreqs,
    allele_freqs_normal: DiscreteAlleleFreqs,
    grid_points: usize
}


impl FlatTumorNormalModel {
    pub fn new(ploidy: u32) -> Self {
        let allele_freqs = (0..ploidy + 1).map(|m| AlleleFreq(m as f64 / ploidy as f64)).collect_vec();
        FlatTumorNormalModel {
            allele_freqs_tumor: AlleleFreq(0.0)..AlleleFreq(1.0),
            allele_freqs_normal: allele_freqs,
            grid_points: 201
        }
    }
}


impl PairModel<ContinuousAlleleFreqs, DiscreteAlleleFreqs> for FlatTumorNormalModel {

    fn prior_prob(&self, _: AlleleFreq, _: AlleleFreq, _: Variant) -> LogProb {
        LogProb::ln_one()
    }

    fn joint_prob<L, O>(
        &self,
        af_tumor: &ContinuousAlleleFreqs,
        af_normal: &DiscreteAlleleFreqs,
        likelihood_tumor: &L,
        likelihood_normal: &O,
        _: Variant,
        _: usize,
        _: usize
    ) -> LogProb where
        L: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb,
        O: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb
    {
        let prob = LogProb::ln_sum_exp(&af_normal.iter().map(|&af_normal| {
            let density = |af_tumor| {
                let af_tumor = AlleleFreq(af_tumor);
                likelihood_tumor(af_tumor, Some(af_normal))
            };

            let p_tumor = if af_tumor.start == af_tumor.end {
                density(*af_tumor.start)
            } else {
                LogProb::ln_simpsons_integrate_exp(&density, *af_tumor.start, *af_tumor.end, self.grid_points)
            };
            let p_normal = likelihood_normal(af_normal, None);
            let prob = p_tumor + p_normal;

            prob
        }).collect_vec());

        prob
    }

    fn marginal_prob<L, O>(
        &self,
        likelihood_tumor: &L,
        likelihood_normal: &O,
        variant: Variant,
        n_obs_tumor: usize,
        n_obs_normal: usize
    ) -> LogProb where
        L: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb,
        O: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb
    {
        let p = self.joint_prob(
            self.allele_freqs().0,
            self.allele_freqs().1,
            likelihood_tumor,
            likelihood_normal,
            variant,
            n_obs_tumor,
            n_obs_normal
        ).ln_add_exp(
            // add prob for allele frequency zero (the density is non-continuous there)
            self.joint_prob(
                &(AlleleFreq(0.0)..AlleleFreq(0.0)),
                &vec![AlleleFreq(0.0)],
                likelihood_tumor,
                likelihood_normal,
                variant,
                n_obs_tumor,
                n_obs_normal
            )
        );
        p
    }

    fn map<L, O>(
        &self,
        likelihood_tumor: &L,
        likelihood_normal: &O,
        _: Variant,
        _: usize,
        _: usize
    ) -> (AlleleFreq, AlleleFreq) where
        L: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb,
        O: Fn(AlleleFreq, Option<AlleleFreq>) -> LogProb
    {
        let af_case = linspace(*self.allele_freqs_tumor.start, *self.allele_freqs_tumor.end, self.grid_points);
        let (_, (map_normal, map_tumor)) = self.allele_freqs().1.iter().cartesian_product(af_case).minmax_by_key(
            |&(&af_normal, af_tumor)| {
                let af_tumor = AlleleFreq(af_tumor);
                let p = likelihood_tumor(af_tumor, Some(af_normal)) +
                        likelihood_normal(af_normal, None);
                debug!("L(f_t={}, f_n={})={}", *af_tumor, *af_normal, *p);
                NotNaN::new(*p).expect("posterior probability is NaN")
            }
        ).into_option().expect("prior has empty allele frequency spectrum");

        (AlleleFreq(map_tumor), *map_normal)
    }

    fn allele_freqs(&self) -> (&ContinuousAlleleFreqs, &DiscreteAlleleFreqs) {
        (&self.allele_freqs_tumor, &self.allele_freqs_normal)
    }
}
