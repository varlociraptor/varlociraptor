use itertools::{Itertools, linspace};
use ordered_float::NotNaN;
use bio::stats::LogProb;

use model::{Variant, ContinuousAlleleFreqs, DiscreteAlleleFreqs, AlleleFreq};

use priors::PairModel;

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
            grid_points: 200
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
        _: Variant
    ) -> LogProb where
        L: Fn(AlleleFreq, AlleleFreq) -> LogProb,
        O: Fn(AlleleFreq, AlleleFreq) -> LogProb
    {
        let prob = LogProb::ln_sum_exp(&af_normal.iter().map(|&af_normal| {
            let density = |af_tumor| {
                let af_tumor = AlleleFreq(af_tumor);
                likelihood_tumor(af_tumor, af_normal)
            };

            let p_tumor = if af_tumor.start == af_tumor.end {
                density(*af_tumor.start)
            } else {
                LogProb::ln_simpsons_integrate_exp(&density, *af_tumor.start, *af_tumor.end, self.grid_points)
            };
            let p_normal = likelihood_normal(af_normal, AlleleFreq(0.0));
            let prob = p_tumor + p_normal;

            prob
        }).collect_vec());

        prob
    }

    fn marginal_prob<L, O>(
        &self,
        likelihood_tumor: &L,
        likelihood_normal: &O,
        variant: Variant
    ) -> LogProb where
        L: Fn(AlleleFreq, AlleleFreq) -> LogProb,
        O: Fn(AlleleFreq, AlleleFreq) -> LogProb
    {
        let p = self.joint_prob(
            self.allele_freqs().0,
            self.allele_freqs().1,
            likelihood_tumor,
            likelihood_normal,
            variant
        ).ln_add_exp(
            // add prob for allele frequency zero (the density is non-continuous there)
            self.joint_prob(
                &(AlleleFreq(0.0)..AlleleFreq(0.0)),
                &vec![AlleleFreq(0.0)],
                likelihood_tumor,
                likelihood_normal,
                variant
            )
        );
        p
    }

    fn map<L, O>(
        &self,
        likelihood_tumor: &L,
        likelihood_normal: &O,
        _: Variant
    ) -> (AlleleFreq, AlleleFreq) where
        L: Fn(AlleleFreq, AlleleFreq) -> LogProb,
        O: Fn(AlleleFreq, AlleleFreq) -> LogProb
    {
        let af_case = linspace(*self.allele_freqs_tumor.start, *self.allele_freqs_tumor.end, self.grid_points);
        let (_, (map_normal, map_tumor)) = self.allele_freqs().1.iter().cartesian_product(af_case).minmax_by_key(
            |&(&af_normal, af_tumor)| {
                let af_tumor = AlleleFreq(af_tumor);
                let p = likelihood_tumor(af_tumor, af_normal) +
                        likelihood_normal(af_normal, AlleleFreq(0.0));
                //println!("af {} vs {} = {} (prior={} tumor={} normal={})", *af_tumor, af_normal, *p, *self.prior_prob(af_tumor, af_normal, variant), *likelihood_tumor(af_tumor, af_normal), *likelihood_normal(af_normal, AlleleFreq(0.0)));
                NotNaN::new(*p).expect("posterior probability is NaN")
            }
        ).into_option().expect("prior has empty allele frequency spectrum");

        (AlleleFreq(map_tumor), *map_normal)
    }

    fn allele_freqs(&self) -> (&ContinuousAlleleFreqs, &DiscreteAlleleFreqs) {
        (&self.allele_freqs_tumor, &self.allele_freqs_normal)
    }
}
