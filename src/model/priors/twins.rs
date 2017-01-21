use itertools::{Itertools, linspace};
use ordered_float::NotNaN;
use bio::stats::LogProb;

use model::{Variant, ContinuousAlleleFreqs, DiscreteAlleleFreqs, AlleleFreq};

use priors::PairModel;


pub struct MonozygoticTwinsModel {
    allele_freqs: DiscreteAlleleFreqs,
    inner: InfiniteSitesNeutralVariationModel
}


impl MonozygoticTwinsModel {
    pub fn new(ploidy: u32, heterozygosity: Prob) -> Self {
        let allele_freqs = (0..ploidy + 1).map(|m| AlleleFreq(m as f64 / ploidy as f64)).collect_vec();
        MonozygoticTwinsModel {
            allele_freqs: allele_freqs,
            heterozygosity: heterozygosity
        }
    }
}


impl PairModel<DiscreteAlleleFreqs, DiscreteAlleleFreqs> for MonozygoticTwinsModel {

    fn prior_prob(&self, af_first: AlleleFreq, af_second: AlleleFreq, _: Variant) -> LogProb {
        if af_first == af_second {
            self.inner.prior_prob(af_first)
        } else {
            // TODO this is probably too extreme
            LogProb::ln_zero()
        }
    }

    fn joint_prob<L, O>(
        &self,
        af_first: &DiscreteAlleleFreqs,
        af_second: &DiscreteAlleleFreqs,
        likelihood_first: &L,
        likelihood_second: &O,
        _: Variant
    ) -> LogProb where
        L: Fn(AlleleFreq, AlleleFreq) -> LogProb,
        O: Fn(AlleleFreq, AlleleFreq) -> LogProb
    {
        // TODO
        let p_second = LogProb::ln_sum_exp(&af_second.iter().map(|&af_second| {
            likelihood_second(af_second, AlleleFreq(0.0))
        }).collect_vec());

        let prob = LogProb::ln_sum_exp(&af_first.iter().map(|&af_first| {
            let p_first = likelihood_first(af_first, AlleleFreq(0.0));
            let prob = p_first + p_second;

            prob
        }).collect_vec());

        prob
    }

    fn marginal_prob<L, O>(
        &self,
        likelihood_first: &L,
        likelihood_second: &O,
        variant: Variant
    ) -> LogProb where
        L: Fn(AlleleFreq, AlleleFreq) -> LogProb,
        O: Fn(AlleleFreq, AlleleFreq) -> LogProb
    {
        // TODO
        let p = self.joint_prob(
            self.allele_freqs().0,
            self.allele_freqs().1,
            likelihood_first,
            likelihood_second,
            variant
        );

        p
    }

    fn map<L, O>(
        &self,
        likelihood_first: &L,
        likelihood_second: &O,
        _: Variant
    ) -> (AlleleFreq, AlleleFreq) where
        L: Fn(AlleleFreq, AlleleFreq) -> LogProb,
        O: Fn(AlleleFreq, AlleleFreq) -> LogProb
    {
        // TODO
        fn calc_map<L: Fn(AlleleFreq, AlleleFreq) -> LogProb>(likelihood: &L, afs: &DiscreteAlleleFreqs) -> AlleleFreq {
            let (_, map) = afs.iter().minmax_by_key(|&af| {
                let p = likelihood(*af, AlleleFreq(0.0));
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
