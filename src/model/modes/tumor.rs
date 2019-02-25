// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;

use bio::stats::bayesian::model::{Likelihood, Posterior};
use bio::stats::LogProb;

use crate::model::likelihood;
use crate::model::sample::Pileup;
use crate::model::{AlleleFreq, ContinuousAlleleFreqs};
use crate::model;

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
    type BaseEvent = Vec<likelihood::Event>;
    type Event = Vec<model::Event<ContinuousAlleleFreqs>>;
    type Data = Vec<Pileup>;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        event: &Self::Event,
        pileups: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        let n_obs_tumor = pileups.tumor().len();
        let n_obs_normal = pileups.normal().len();
        let grid_points_normal = 5;
        let grid_points_tumor = Self::grid_points_tumor(n_obs_tumor);
        let strand_bias_tumor = event.tumor().strand_bias;
        let strand_bias_normal = event.normal().strand_bias;

        let mut density = |af_normal| {
            let af_normal = AlleleFreq(af_normal);

            let p = {
                let mut joint_density = |af_tumor| {
                    let af_tumor = AlleleFreq(af_tumor);
                    let p = joint_prob(
                        &TumorNormalPair {
                            tumor: likelihood::Event {
                                allele_freq: af_tumor,
                                strand_bias: strand_bias_normal,
                            },
                            normal: likelihood::Event {
                                allele_freq: af_normal,
                                strand_bias: strand_bias_tumor,
                            }
                        }
                        .into(),
                        pileups,
                    );
                    p
                };

                if event.tumor().allele_freqs.is_singleton() {
                    joint_density(*event.tumor().allele_freqs.start)
                } else {
                    LogProb::ln_simpsons_integrate_exp(
                        joint_density,
                        *event.tumor().allele_freqs.observable_min(n_obs_tumor),
                        *event.tumor().allele_freqs.observable_max(n_obs_tumor),
                        grid_points_tumor,
                    )
                }
            };

            p
        };

        let prob = if event.normal().allele_freqs.is_singleton() {
            density(*event.normal().allele_freqs.start)
        } else {
            LogProb::ln_simpsons_integrate_exp(
                density,
                *event.normal().allele_freqs.observable_min(n_obs_normal),
                *event.normal().allele_freqs.observable_max(n_obs_normal),
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
    type Event = Vec<likelihood::Event>;
    type Data = Vec<Pileup>;

    fn compute(&self, event: &Self::Event, pileups: &Self::Data) -> LogProb {
        let p = self.tumor_likelihood.compute(event, &pileups.tumor())
            + self
                .normal_likelihood
                .compute(&event.normal(), &pileups.normal());

        p
    }
}
