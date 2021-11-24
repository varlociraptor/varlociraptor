use crate::{
    calling::haplotypes::{KallistoEstimate, KallistoEstimates},
    variants::model::AlleleFreq,
};
use bio::stats::bayesian;
use bio::stats::{bayesian::model, LogProb};

use crate::utils::adaptive_integration;
use ordered_float::NotNaN;
use statrs::function::beta::ln_beta;
use std::collections::HashMap;
use std::mem;

#[derive(Hash, PartialEq, Eq, Clone, Debug, Derefable)]
pub(crate) struct HaplotypeFractions(#[deref] Vec<AlleleFreq>);

impl HaplotypeFractions {
    // pub(crate) fn likely(kallisto_estimates: &KallistoEstimates) -> Vec<Self> {
    //     // TODO: return all combinations of haplotype fractions that are somehow likely
    //     // given the callisto estimates. E.g., omit fractions > 0 for haplotypes that
    //     // do not occur at all.
    //     todo!();
    // }
}

#[derive(Debug, new)]
pub(crate) struct Marginal {
    n_haplotypes: usize,
}

impl Marginal {
    pub(crate) fn calc_marginal<
        F: FnMut(&<Self as model::Marginal>::Event, &<Self as model::Marginal>::Data) -> LogProb,
    >(
        &self,
        data: &Data,
        haplotype_index: usize,
        fractions: &mut Vec<AlleleFreq>,
        joint_prob: &mut F,
    ) -> LogProb {
        if haplotype_index == self.n_haplotypes {
            let event = HaplotypeFractions(fractions.to_vec());
            joint_prob(&event, data)
        } else {
            let density = |fraction| {
                let mut fractions = fractions.clone();
                fractions.push(fraction);
                self.calc_marginal(data, haplotype_index + 1, &mut fractions, joint_prob)
            };
            adaptive_integration::ln_integrate_exp(
                density,
                NotNaN::new(0.0).unwrap(),
                NotNaN::new(1.0).unwrap(),
                NotNaN::new(0.1).unwrap(),
            )
        }
    }
}

impl model::Marginal for Marginal {
    type Event = HaplotypeFractions;
    type Data = Data;
    type BaseEvent = HaplotypeFractions;

    fn compute<F: FnMut(&Self::Event, &Self::Data) -> LogProb>(
        &self,
        data: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        let mut fractions: Vec<AlleleFreq> = Vec::new();
        self.calc_marginal(data, 0, &mut fractions, joint_prob)
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
        self.compute_kallisto(event, data, payload)
            + self.compute_varlociraptor(event, data, payload)
    }
}

impl Likelihood {
    fn compute_kallisto(
        &self,
        event: &HaplotypeFractions,
        data: &Data,
        cache: &mut Cache,
    ) -> LogProb {
        // TODO compute likelihood using neg_binom on the counts and dispersion
        // in the data and the fractions in the events.
        //Later: use the cache to avoid redundant computations.
        event
            .iter()
            .zip(data.kallisto_estimates.iter())
            .map(|(fraction, estimate)| {
                neg_binom(
                    estimate.count,
                    NotNaN::into_inner(*fraction),
                    estimate.dispersion,
                )
            })
            .sum()
    }

    fn compute_varlociraptor(
        &self,
        event: &HaplotypeFractions,
        data: &Data,
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
    type Event = HaplotypeFractions;

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

// TODO move into model
pub(crate) fn neg_binom(x: f64, mu: f64, theta: f64) -> LogProb {
    let n = 1.0 / theta;
    let p = n / (n + mu);
    let mut p1 = if n > 0.0 { n * p.ln() } else { 0.0 };
    let mut p2 = if x > 0.0 { x * (1.0 - p).ln() } else { 0.0 };
    let b = ln_beta(x + 1.0, n);
    if p1 < p2 {
        mem::swap(&mut p1, &mut p2);
    }
    LogProb((p1 - b + p2) - (x + n).ln())
}
