// use std::f64;
//

// use ordered_float::NotNan;
// use bio::stats::{LogProb, Prob};
// use bio::stats::combinatorics::combinations;
//
// use model::{Variant, DiscreteAlleleFreqs, AlleleFreq};
// use priors::{Model, InfiniteSitesNeutralVariationModel, PairModel};

/*
pub(crate) struct NormalNormalModel {
    inner: InfiniteSitesNeutralVariationModel,
    ploidy: u32,
    allele_freqs: DiscreteAlleleFreqs
}


impl NormalNormalModel {
    pub(crate) fn new(ploidy: u32, heterozygosity: Prob) -> Self {
        NormalNormalModel {
            inner: InfiniteSitesNeutralVariationModel::new(2, ploidy, heterozygosity),
            ploidy: ploidy,
            allele_freqs: allele_freqs(ploidy)
        }
    }
}


impl PairModel<DiscreteAlleleFreqs, DiscreteAlleleFreqs> for NormalNormalModel {
    fn prior_prob(&self, af_first: AlleleFreq, af_second: AlleleFreq, _: &Variant) -> LogProb {
        let sample_prior = |af: AlleleFreq| {
            let m = *af * self.ploidy as f64;
            let valid_genotypes = combinations(self.ploidy, m);
            let all_genotypes = combinations()
        };
    }

    fn joint_prob<L, O>(
        &self,
        af_first: &DiscreteAlleleFreqs,
        af_second: &DiscreteAlleleFreqs,
        likelihood_first: &L,
        likelihood_second: &O,
        _: &Variant
    ) -> LogProb where
        L: Fn(AlleleFreq, AlleleFreq) -> LogProb,
        O: Fn(AlleleFreq, AlleleFreq) -> LogProb
    {
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
        variant: &Variant
    ) -> LogProb where
        L: Fn(AlleleFreq, AlleleFreq) -> LogProb,
        O: Fn(AlleleFreq, AlleleFreq) -> LogProb
    {
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
        _: &Variant
    ) -> (AlleleFreq, AlleleFreq) where
        L: Fn(AlleleFreq, AlleleFreq) -> LogProb,
        O: Fn(AlleleFreq, AlleleFreq) -> LogProb
    {
        fn calc_map<L: Fn(AlleleFreq, AlleleFreq) -> LogProb>(likelihood: &L, afs: &DiscreteAlleleFreqs) -> AlleleFreq {
            let (_, map) = afs.iter().minmax_by_key(|&af| {
                let p = likelihood(*af, AlleleFreq(0.0));
                NotNan::new(*p).expect("probability is NaN")
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
}*/
