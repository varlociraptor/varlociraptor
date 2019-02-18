use std::cmp;

use bio::stats::bayesian::model::Posterior;
use bio::stats::LogProb;

use crate::model::sample::Pileup;
use crate::model::{AlleleFreq, ContinuousAlleleFreqs};

/// Posterior model for Tumor-Normal sample pairs.
pub struct TumorNormalPosterior {}

impl TumorNormalPosterior {
    fn grid_points_tumor(n_obs: usize) -> usize {
        let mut n = cmp::max(n_obs + 1, 5);
        if n % 2 == 0 {
            n += 1;
        }
        n
    }
}

impl Posterior for TumorNormalPosterior {
    type BaseEvent = (AlleleFreq, AlleleFreq);
    type Event = (ContinuousAlleleFreqs, ContinuousAlleleFreqs);
    type Data = Vec<Pileup>;

    fn compute<F: FnMut(&(AlleleFreq, AlleleFreq), &Self::Data) -> LogProb>(
        &self,
        allele_freqs: &(ContinuousAlleleFreqs, ContinuousAlleleFreqs),
        pileups: &Self::Data,
        joint_prob: &F,
    ) -> LogProb {
        let (af_tumor, af_normal) = allele_freqs;
        assert_eq!(pileups.len(), 2, "invalid number of pileups");
        let (ref pileup_tumor, ref pileup_normal) = (pileups[0], pileups[1]);

        let n_obs_tumor = pileup_tumor.len();
        let n_obs_normal = pileup_normal.len();
        let grid_points_normal = 5;
        let grid_points_tumor = Self::grid_points_tumor(n_obs_tumor);

        let mut density = |af_normal| {
            let af_normal = AlleleFreq(af_normal);

            let p = {
                let mut tumor_density = |af_tumor| {
                    let af_tumor = AlleleFreq(af_tumor);
                    let p = joint_prob(&(af_tumor, af_normal), pileups);
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

            p
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
}
