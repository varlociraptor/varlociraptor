use std::cmp;

use bio::stats::bayesian::model::{Likelihood, Posterior, Prior};
use bio::stats::LogProb;

use crate::model::likelihood;
use crate::model::sample::Pileup;
use crate::model::{AlleleFreq, ContinuousAlleleFreqs};

#[derive(Debug, Clone)]
pub struct TumorNormalPair<T> {
    pub tumor: T,
    pub normal: T,
}

impl<T> Into<Vec<T>> for TumorNormalPair<T> {
    fn into(self) -> Vec<T> {
        vec![self.tumor, self.normal]
    }
}

pub trait TumorNormalPairView<T> {
    fn tumor(&self) -> &T;

    fn normal(&self) -> &T;
}

impl<T> TumorNormalPairView<T> for Vec<T> {
    fn tumor(&self) -> &T {
        &self[0]
    }

    fn normal(&self) -> &T {
        &self[1]
    }
}

/// Posterior model for Tumor-Normal sample pairs.
#[derive(Debug, Clone, Default)]
pub struct TumorNormalPosterior {}

impl TumorNormalPosterior {
    pub fn new() -> Self {
        TumorNormalPosterior {}
    }

    fn grid_points_tumor(n_obs: usize) -> usize {
        let mut n = cmp::max(n_obs + 1, 5);
        if n % 2 == 0 {
            n += 1;
        }
        n
    }
}

impl Posterior for TumorNormalPosterior {
    type BaseEvent = Vec<AlleleFreq>;
    type Event = Vec<ContinuousAlleleFreqs>;
    type Data = Vec<Pileup>;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        allele_freqs: &Self::Event,
        pileups: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        let n_obs_tumor = pileups.tumor().len();
        let n_obs_normal = pileups.normal().len();
        let grid_points_normal = 5;
        let grid_points_tumor = Self::grid_points_tumor(n_obs_tumor);

        let mut density = |af_normal| {
            let af_normal = AlleleFreq(af_normal);

            let p = {
                let mut tumor_density = |af_tumor| {
                    let af_tumor = AlleleFreq(af_tumor);
                    let p = joint_prob(
                        &TumorNormalPair {
                            tumor: af_tumor,
                            normal: af_normal,
                        }
                        .into(),
                        pileups,
                    );
                    p
                };

                if allele_freqs.tumor().is_singleton() {
                    tumor_density(*allele_freqs.tumor().start)
                } else {
                    LogProb::ln_simpsons_integrate_exp(
                        tumor_density,
                        *allele_freqs.tumor().observable_min(n_obs_tumor),
                        *allele_freqs.tumor().observable_max(n_obs_tumor),
                        grid_points_tumor,
                    )
                }
            };

            p
        };

        let prob = if allele_freqs.normal().is_singleton() {
            density(*allele_freqs.normal().start)
        } else {
            LogProb::ln_simpsons_integrate_exp(
                density,
                *allele_freqs.normal().observable_min(n_obs_normal),
                *allele_freqs.normal().observable_max(n_obs_normal),
                grid_points_normal,
            )
        };

        prob
    }
}

#[derive(Debug, Clone, Default)]
pub struct TumorNormalLikelihood {
    tumor_likelihood: likelihood::ContaminatedSampleLikelihoodModel,
    normal_likelihood: likelihood::SampleLikelihoodModel,
}

impl TumorNormalLikelihood {
    pub fn new(purity: f64) -> Self {
        TumorNormalLikelihood {
            tumor_likelihood: likelihood::ContaminatedSampleLikelihoodModel::new(purity),
            normal_likelihood: likelihood::SampleLikelihoodModel::new(),
        }
    }
}

impl Likelihood for TumorNormalLikelihood {
    type Event = Vec<AlleleFreq>;
    type Data = Vec<Pileup>;

    fn compute(&self, allele_freq: &Self::Event, pileups: &Self::Data) -> LogProb {
        self.tumor_likelihood.compute(allele_freq, &pileups.tumor())
            + self
                .normal_likelihood
                .compute(&allele_freq.normal(), &pileups.normal())
    }
}
