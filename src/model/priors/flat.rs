use std::cmp;

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
                })
                .collect_vec(),
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
                })
                .into_option()
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
    allele_freqs_normal_germline: DiscreteAlleleFreqs,
    allele_freqs_normal_somatic: ContinuousAlleleFreqs,
}

impl FlatTumorNormalModel {
    pub fn new(ploidy: u32) -> Self {
        let allele_freqs_normal_germline = DiscreteAlleleFreqs::feasible(ploidy);
        let allele_freqs_normal_somatic =
            ContinuousAlleleFreqs::exclusive(0.0..*allele_freqs_normal_germline[1]);
        FlatTumorNormalModel {
            // TODO how to sync this with zero_offset of other events
            allele_freqs_tumor: ContinuousAlleleFreqs::inclusive(0.0..1.0),
            allele_freqs_normal_germline: allele_freqs_normal_germline,
            allele_freqs_normal_somatic: allele_freqs_normal_somatic,
        }
    }

    fn grid_points_tumor(n_obs: usize) -> usize {
        let mut n = cmp::max(n_obs + 1, 5);
        if n % 2 == 0 {
            n += 1;
        }
        n
    }
}

impl PairModel<ContinuousAlleleFreqs, ContinuousAlleleFreqs> for FlatTumorNormalModel {
    fn prior_prob(&self, _: AlleleFreq, _: AlleleFreq, _: &Variant) -> LogProb {
        LogProb::ln_one()
    }

    fn joint_prob(
        &self,
        af_tumor: &ContinuousAlleleFreqs,
        af_normal: &ContinuousAlleleFreqs,
        pileup: &mut PairPileup<ContinuousAlleleFreqs, ContinuousAlleleFreqs, Self>,
    ) -> LogProb {
        let n_obs_tumor = pileup.case.len();
        let n_obs_normal = pileup.control.len();
        let grid_points_normal = 5;
        let grid_points_tumor = Self::grid_points_tumor(n_obs_tumor);

        let mut density = |af_normal| {
            let af_normal = AlleleFreq(af_normal);
            let p_tumor = {
                let mut tumor_density = |af_tumor| {
                    let af_tumor = AlleleFreq(af_tumor);
                    let p = pileup.case_likelihood(af_tumor, Some(af_normal));
                    p
                };

                if af_tumor.is_singleton() {
                    tumor_density(*af_tumor.start)
                } else {
                    LogProb::ln_simpsons_integrate_exp(
                        tumor_density,
                        *af_tumor.observable_min(n_obs_tumor),
                        *af_tumor.observable_max(n_obs_tumor),
                        grid_points_tumor,
                    )
                }
            };

            let p_normal = pileup.control_likelihood(af_normal, None);
            let prob = p_tumor + p_normal;

            prob
        };

        let prob = if af_normal.is_singleton() {
            density(*af_normal.start)
        } else {
            LogProb::ln_simpsons_integrate_exp(
                density,
                *af_normal.observable_min(n_obs_normal),
                *af_normal.observable_max(n_obs_normal),
                grid_points_normal,
            )
        };

        prob
    }

    fn marginal_prob(
        &self,
        _pileup: &mut PairPileup<ContinuousAlleleFreqs, ContinuousAlleleFreqs, Self>,
    ) -> LogProb {
        panic!("deprecated");
    }

    fn map(
        &self,
        pileup: &mut PairPileup<ContinuousAlleleFreqs, ContinuousAlleleFreqs, Self>,
    ) -> (AlleleFreq, AlleleFreq) {
        let n_obs_tumor = pileup.case.len();
        let n_obs_normal = pileup.control.len();

        let af_tumor = linspace(
            *self.allele_freqs_tumor.observable_min(n_obs_tumor),
            *self.allele_freqs_tumor.observable_max(n_obs_tumor),
            Self::grid_points_tumor(n_obs_tumor),
        );
        // build the entire normal allele freq spectrum
        let af_normal = linspace(
            *self.allele_freqs_normal_somatic.observable_min(n_obs_normal),
            *self.allele_freqs_normal_somatic.observable_max(n_obs_normal),
            5,
        ).chain(self.allele_freqs_normal_germline.iter().map(|af| **af));

        let (map_normal, map_tumor) = af_normal
            .cartesian_product(af_tumor)
            .max_by_key(|&(af_normal, af_tumor)| {
                let af_tumor = AlleleFreq(af_tumor);
                let af_normal = AlleleFreq(af_normal);
                let p = pileup.case_likelihood(af_tumor, Some(af_normal))
                    + pileup.control_likelihood(af_normal, None);
                NotNan::new(*p).expect("posterior probability is NaN")
            }).expect("bug: prior has empty allele frequency spectrum");

        (AlleleFreq(map_tumor), AlleleFreq(map_normal))
    }

    fn allele_freqs(&self) -> (&ContinuousAlleleFreqs, &ContinuousAlleleFreqs) {
        panic!("deprecated");
    }
}
