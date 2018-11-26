use bio::stats::LogProb;
use itertools::Itertools;
use itertools_num::linspace;
use ordered_float::NotNan;

use model::{AlleleFreq, ContinuousAlleleFreqs, DiscreteAlleleFreqs, PairPileup, Variant};

use priors::PairModel;

pub struct FlatNormalNormalModel {
    allele_freqs: DiscreteAlleleFreqs,
}

impl FlatNormalNormalModel {
    pub fn new(ploidy: u32) -> Self {
        FlatNormalNormalModel {
            allele_freqs: DiscreteAlleleFreqs::feasible(ploidy),
        }
    }
}

impl PairModel<DiscreteAlleleFreqs, DiscreteAlleleFreqs> for FlatNormalNormalModel {
    fn prior_prob(&self, _: AlleleFreq, _: AlleleFreq, _: &Variant) -> LogProb {
        LogProb::ln_one()
    }

    fn joint_prob(
        &self,
        af_first: &DiscreteAlleleFreqs,
        af_second: &DiscreteAlleleFreqs,
        pileup: &mut PairPileup<DiscreteAlleleFreqs, DiscreteAlleleFreqs, Self>,
    ) -> LogProb {
        let p_second = LogProb::ln_sum_exp(
            &af_second
                .iter()
                .map(|&af_second| pileup.control_likelihood(af_second, None))
                .collect_vec(),
        );

        let prob = LogProb::ln_sum_exp(
            &af_first
                .iter()
                .map(|&af_first| {
                    let p_first = pileup.case_likelihood(af_first, None);
                    let prob = p_first + p_second;

                    prob
                }).collect_vec(),
        );

        prob
    }

    fn marginal_prob(
        &self,
        pileup: &mut PairPileup<DiscreteAlleleFreqs, DiscreteAlleleFreqs, Self>,
    ) -> LogProb {
        let p = self.joint_prob(self.allele_freqs().0, self.allele_freqs().1, pileup);

        p
    }

    fn map(
        &self,
        pileup: &mut PairPileup<DiscreteAlleleFreqs, DiscreteAlleleFreqs, Self>,
    ) -> (AlleleFreq, AlleleFreq) {
        fn calc_map<L: FnMut(AlleleFreq) -> LogProb>(
            likelihood: &mut L,
            afs: &DiscreteAlleleFreqs,
        ) -> AlleleFreq {
            let (_, map) = afs
                .iter()
                .minmax_by_key(|&af| {
                    let p = likelihood(*af);
                    NotNan::new(*p).expect("probability is NaN")
                }).into_option()
                .expect("prior has empty allele frequency spectrum");
            *map
        }

        let map_first: AlleleFreq;
        let map_second: AlleleFreq;

        {
            let mut first_likelihood =
                |af_first: AlleleFreq| pileup.case_likelihood(af_first, None);
            map_first = calc_map(&mut first_likelihood, self.allele_freqs().0);
        }
        {
            let mut second_likelihood =
                |af_second: AlleleFreq| pileup.control_likelihood(af_second, None);
            map_second = calc_map(&mut second_likelihood, self.allele_freqs().1);
        }

        (map_first, map_second)
    }

    fn allele_freqs(&self) -> (&DiscreteAlleleFreqs, &DiscreteAlleleFreqs) {
        (&self.allele_freqs, &self.allele_freqs)
    }
}

pub struct FlatTumorNormalModel {
    allele_freqs_tumor: ContinuousAlleleFreqs,
    allele_freqs_normal: DiscreteAlleleFreqs,
    grid_points: usize,
}

impl FlatTumorNormalModel {
    pub fn new(ploidy: u32) -> Self {
        FlatTumorNormalModel {
            allele_freqs_tumor: ContinuousAlleleFreqs::inclusive(0.0..1.0),
            allele_freqs_normal: DiscreteAlleleFreqs::feasible(ploidy),
            grid_points: 201,
        }
    }
}

impl PairModel<ContinuousAlleleFreqs, DiscreteAlleleFreqs> for FlatTumorNormalModel {
    fn prior_prob(&self, _: AlleleFreq, _: AlleleFreq, _: &Variant) -> LogProb {
        LogProb::ln_one()
    }

    fn joint_prob(
        &self,
        af_tumor: &ContinuousAlleleFreqs,
        af_normal: &DiscreteAlleleFreqs,
        pileup: &mut PairPileup<ContinuousAlleleFreqs, DiscreteAlleleFreqs, Self>,
    ) -> LogProb {
        let grid_points = self.grid_points;
        let prob = LogProb::ln_sum_exp(
            &af_normal
                .iter()
                .map(|&af_normal| {
                    let p_tumor: LogProb;
                    {
                        let mut density = |af_tumor| {
                            let af_tumor = AlleleFreq(af_tumor);
                            pileup.case_likelihood(af_tumor, Some(af_normal))
                        };

                        p_tumor = if af_tumor.start == af_tumor.end {
                            density(*af_tumor.start)
                        } else {
                            LogProb::ln_simpsons_integrate_exp(
                                density,
                                *af_tumor.start,
                                *af_tumor.end,
                                grid_points,
                            )
                        };
                    }
                    let p_normal = pileup.control_likelihood(af_normal, None);
                    let prob = p_tumor + p_normal;

                    prob
                }).collect_vec(),
        );

        prob
    }

    fn marginal_prob(
        &self,
        pileup: &mut PairPileup<ContinuousAlleleFreqs, DiscreteAlleleFreqs, Self>,
    ) -> LogProb {
        let p = self
            .joint_prob(self.allele_freqs().0, self.allele_freqs().1, pileup)
            .ln_add_exp(
                // add prob for allele frequency zero (the density is non-continuous there)
                self.joint_prob(
                    &ContinuousAlleleFreqs::inclusive(0.0..0.0),
                    &DiscreteAlleleFreqs::new(vec![AlleleFreq(0.0)]),
                    pileup,
                ),
            );
        p
    }

    fn map(
        &self,
        pileup: &mut PairPileup<ContinuousAlleleFreqs, DiscreteAlleleFreqs, Self>,
    ) -> (AlleleFreq, AlleleFreq) {
        let af_case = linspace(
            *self.allele_freqs_tumor.start,
            *self.allele_freqs_tumor.end,
            pileup.case.len() + 1,
        );

        let (_, (map_normal, map_tumor)) = self
            .allele_freqs()
            .1
            .iter()
            .cartesian_product(af_case)
            .minmax_by_key(|&(&af_normal, af_tumor)| {
                let af_tumor = AlleleFreq(af_tumor);
                let p = pileup.case_likelihood(af_tumor, Some(af_normal))
                    + pileup.control_likelihood(af_normal, None);
                NotNan::new(*p).expect("posterior probability is NaN")
            }).into_option()
            .expect("prior has empty allele frequency spectrum");

        (AlleleFreq(map_tumor), *map_normal)
    }

    fn allele_freqs(&self) -> (&ContinuousAlleleFreqs, &DiscreteAlleleFreqs) {
        (&self.allele_freqs_tumor, &self.allele_freqs_normal)
    }
}
