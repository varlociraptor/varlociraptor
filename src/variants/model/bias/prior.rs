use bio::stats::{
    bayesian::model::{Likelihood, Posterior, Prior},
    LogProb,
};
use itertools::Itertools;

use crate::{
    grammar::VAFSpectrum,
    variants::{
        evidence::observation::{Observation, ReadPosition},
        model::AlleleFreq,
        sample::Pileup,
    },
};

pub(crate) fn model<L>(
    likelihood_func: L,
) -> bio::stats::bayesian::Model<HyperLikelihood<L>, HyperPrior, HyperPosterior>
where
    L: Fn(AlleleFreq, &Observation<ReadPosition>) -> LogProb,
{
    bio::stats::bayesian::Model::new(
        HyperLikelihood::new(likelihood_func),
        HyperPrior::new(),
        HyperPosterior::new(),
    )
}

#[derive(new, Debug, Clone)]
pub(crate) struct HyperLikelihood<L>
where
    L: Fn(AlleleFreq, &Observation<ReadPosition>) -> LogProb,
{
    likelihood: L,
}

impl<L> Likelihood for HyperLikelihood<L>
where
    L: Fn(AlleleFreq, &Observation<ReadPosition>) -> LogProb,
{
    type Event = AlleleFreq;
    type Data = Pileup;

    fn compute(&self, allele_freq: &AlleleFreq, data: &Pileup, _payload: &mut ()) -> LogProb {
        data.iter()
            .map(|obs| (&self.likelihood)(*allele_freq, obs))
            .sum()
    }
}

#[derive(new, Debug, Clone)]
pub(crate) struct HyperPrior;

impl Prior for HyperPrior {
    type Event = AlleleFreq;

    fn compute(&self, _: &Self::Event) -> LogProb {
        LogProb::ln_one()
    }
}

#[derive(new, Debug, Clone)]
pub(crate) struct HyperPosterior;

impl Posterior for HyperPosterior {
    type BaseEvent = AlleleFreq;
    type Event = VAFSpectrum;
    type Data = Pileup;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        event: &Self::Event,
        data: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        match event {
            VAFSpectrum::Range(vafs) => LogProb::ln_simpsons_integrate_exp(
                |_, vaf| joint_prob(&AlleleFreq(vaf), data),
                *vafs.start,
                *vafs.end,
                11,
            ),
            VAFSpectrum::Set(vafs) => {
                LogProb::ln_sum_exp(&vafs.iter().map(|vaf| joint_prob(vaf, data)).collect_vec())
            }
        }
    }
}
