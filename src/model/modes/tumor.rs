use std::cmp;

use bio::stats::bayesian::model::{Likelihood, Posterior, Prior};
use bio::stats::LogProb;

use crate::model::likelihood;
use crate::model::sample::Pileup;
use crate::model::{AlleleFreq, ContinuousAlleleFreqs};
use crate::Event;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default)]
pub struct TumorNormalEvent {
    pub name: String,
    pub tumor: ContinuousAlleleFreqs,
    pub normal: ContinuousAlleleFreqs,
}

impl Event for TumorNormalEvent {
    fn name(&self) -> &str {
        &self.name
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct TumorNormalAlleleFreq {
    tumor: AlleleFreq,
    normal: AlleleFreq,
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
    type BaseEvent = TumorNormalAlleleFreq;
    type Event = TumorNormalEvent;
    type Data = Vec<Pileup>;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        allele_freqs: &TumorNormalEvent,
        pileups: &Vec<Pileup>,
        joint_prob: &mut F,
    ) -> LogProb {
        assert_eq!(pileups.len(), 2, "invalid number of pileups");
        let (pileup_tumor, pileup_normal) = (&pileups[0], &pileups[1]);

        let n_obs_tumor = pileup_tumor.len();
        let n_obs_normal = pileup_normal.len();
        let grid_points_normal = 5;
        let grid_points_tumor = Self::grid_points_tumor(n_obs_tumor);

        let mut density = |af_normal| {
            let af_normal = AlleleFreq(af_normal);

            let p = {
                let mut tumor_density = |af_tumor| {
                    let af_tumor = AlleleFreq(af_tumor);
                    let p = joint_prob(
                        &TumorNormalAlleleFreq {
                            tumor: af_tumor,
                            normal: af_normal,
                        },
                        pileups,
                    );
                    p
                };

                if allele_freqs.tumor.is_singleton() {
                    tumor_density(*allele_freqs.tumor.start)
                } else {
                    LogProb::ln_simpsons_integrate_exp(
                        tumor_density,
                        *allele_freqs.tumor.observable_min(n_obs_tumor),
                        *allele_freqs.tumor.observable_max(n_obs_tumor),
                        grid_points_tumor,
                    )
                }
            };

            p
        };

        let prob = if allele_freqs.normal.is_singleton() {
            density(*allele_freqs.normal.start)
        } else {
            LogProb::ln_simpsons_integrate_exp(
                density,
                *allele_freqs.normal.observable_min(n_obs_normal),
                *allele_freqs.normal.observable_max(n_obs_normal),
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
    type Event = TumorNormalAlleleFreq;
    type Data = Vec<Pileup>;

    fn compute(&self, allele_freq: &TumorNormalAlleleFreq, pileups: &Vec<Pileup>) -> LogProb {
        assert_eq!(pileups.len(), 2, "invalid number of pileups");
        let (pileup_tumor, pileup_normal) = (&pileups[0], &pileups[1]);

        self.tumor_likelihood
            .compute(&(allele_freq.tumor, allele_freq.normal), pileup_tumor)
            + self
                .normal_likelihood
                .compute(&allele_freq.normal, pileup_normal)
    }
}

#[derive(Debug, Clone, Default)]
pub struct TumorNormalFlatPrior {}

impl TumorNormalFlatPrior {
    pub fn new() -> Self {
        TumorNormalFlatPrior {}
    }
}

impl Prior for TumorNormalFlatPrior {
    type Event = TumorNormalAlleleFreq;

    fn compute(&self, _event: &TumorNormalAlleleFreq) -> LogProb {
        LogProb::ln_one()
    }
}
