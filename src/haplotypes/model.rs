use bio::stats::{bayesian::model, LogProb};

use crate::{
    calling::haplotypes::{KallistoEstimate, KallistoEstimates},
    variants::model::AlleleFreq,
};

pub(crate) struct HaplotypeFractions(Vec<AlleleFreq>);

impl HaplotypeFractions {
    pub(crate) fn likely(kallisto_estimates: &KallistoEstimates) -> Vec<Self> {
        // TODO: return all combinations of haplotype fractions that are somehow likely
        // given the callisto estimates. E.g., omit fractions > 0 for haplotypes that
        // do not occur at all.
        todo!();
    }
}

#[derive(Debug, new)]
pub(crate) struct Data {
    kallisto_estimates: Vec<KallistoEstimate>,
}

#[derive(Debug, new)]
pub(crate) struct Likelihood;

impl model::Likelihood<Cache> for Likelihood {
    type Event = HaplotypeFractions;
    type Data = Data;

    fn compute(&self, event: &Self::Event, data: &Self::Data, payload: &mut Cache) -> LogProb {
        self.compute_kallisto(event, data, cache) + self.compute_varlociraptor(event, data, cache)
    }
}

impl Likelihood {
    fn compute_kallisto(
        &self,
        event: &Self::Event,
        data: &Self::Data,
        cache: &mut Cache,
    ) -> LogProb {
        // TODO compute likelihood using neg_binom on the counts and dispersion
        // in the data and the fractions in the events.
        // Later: use the cache to avoid redundant computations.
        todo!()
    }

    fn compute_varlociraptor(
        &self,
        event: &Self::Event,
        data: &Self::Data,
        cache: &mut Cache,
    ) -> LogProb {
        // TODO compute likelihood based on Varlociraptor VAFs.
        // Let us postpone this until we have a working version with kallisto only.
        LogProb::ln_one()
    }
}

#[derive(Debug, new)]
pub(crate) struct Prior;

impl model::Prior for Prior {
    type Event;

    fn compute(&self, event: &Self::Event) -> LogProb {
        // flat prior for now
        LogProb::ln_one()
    }
}

#[derive(Debug, new)]
pub(crate) struct Posterior;

impl model::Posterior for Posterior {
    type Event = HaplotypeFractions;

    type BaseEvent = HaplotypeFractions;

    type Data = Data;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        event: &Self::Event,
        data: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        // joint_prob calculates the joint probability from likelihood and prior
        joint_prob(event, data)
    }
}

#[derive(Debug, Derefable, Default)]
pub(crate) struct Cache(#[deref] HashMap<usize, HashMap<AlleleFreq, LogProb>>);
