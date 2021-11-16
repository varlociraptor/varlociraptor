use bio::stats::{LogProb, PHREDProb, Prob, bayesian};
use counter::Counter;
use itertools::Itertools;
use rgsl::randist::{binomial::binomial_pdf, multinomial::multinomial_lnpdf};

use crate::variants::{evidence::observation::{Observation, ReadPosition}, model::AlleleFreq};

pub(crate) type Mapq = u32;

#[derive(Debug, Clone)]
pub(crate) struct ObservedFragments {
    unique: u32,
    total: u32,
}

#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub(crate) enum LocusType {
    Unique,
    Homologous,
}

#[derive(Clone, Debug)]
pub(crate) struct Prior;

impl bayesian::model::Prior for Prior {
    type Event = AlleleFreq;

    fn compute(&self, _event: &Self::Event) -> bio::stats::LogProb {
        LogProb::ln_one() // flat prior
    }
}

#[derive(Clone, Debug)]
pub(crate) struct Posterior;

impl bayesian::model::Posterior for Posterior {
    type Event = LocusType;

    type BaseEvent = AlleleFreq;

    type Data = ObservedFragments;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        event: &Self::Event,
        data: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        match event {
            LocusType::Unique => joint_prob(&AlleleFreq(1.0), data),
            LocusType::Homologous => LogProb::ln_simpsons_integrate_exp(|i, freq| {
                dbg!(freq);
                joint_prob(&AlleleFreq(freq), data)
            }, 0.0, 1.0, 11)
        }
        
    }
}

#[derive(Clone, Debug)]
pub(crate) struct Likelihood;

impl bayesian::model::Likelihood for Likelihood {
    type Event = AlleleFreq;

    type Data = ObservedFragments;

    fn compute(&self, event: &Self::Event, data: &Self::Data, _payload: &mut ()) -> LogProb {
        let p = LogProb::from(Prob(binomial_pdf(data.unique, **event, data.total)));
        dbg!(p);

        p
    }
}

pub(crate) fn adjust_prob_mapping(pileup: &mut [Observation<ReadPosition>]) {
    if !pileup.is_empty() {
        // METHOD: we use a simple multinomial model to represent two cases:
        // (a) the locus is a homolog of something else: this means that MAPQs can be expected to be somehow uniformly distributed.
        // (b) the locus is unique: this means that MAPQs should be all at the maximum.
        // The rationale for the former case is that we always expect some noise of artificially high MAPQs because alleles can be unknown by the mapper
        // and sometimes base errors make a read randomly fit better than it should.

        let mapqs: Vec<_> = pileup
            .iter()
            .map(|obs| {PHREDProb::from(obs.prob_mapping_orig().ln_one_minus_exp()).round() as u32})
            .collect();
        let max_mapq = *mapqs.iter().max().unwrap();
        dbg!(max_mapq);
        assert!(mapqs.len() <= std::u32::MAX as usize);

        let data = ObservedFragments {
            unique: mapqs.iter().filter(|mapq| **mapq == max_mapq).count() as u32,
            total: mapqs.len() as u32,
        };
        dbg!(&data);

        let model = bayesian::model::Model::new(Likelihood, Prior, Posterior);

        let m = model.compute(vec![LocusType::Unique, LocusType::Homologous], &data);

        // Calculate expected value of MAPQ given the two posteriors and their repective MAPQs (0 and max_mapq).
        let adjusted =
            m.posterior(&LocusType::Unique).unwrap() + LogProb::from(PHREDProb(max_mapq as f64));
        dbg!(adjusted);
        for obs in pileup {
            if adjusted < obs.prob_mapping_orig() {
                obs.set_prob_mapping_adj(adjusted);
            }
        }
    }
}
