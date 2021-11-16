use bio::stats::{bayesian, LogProb, PHREDProb};
use counter::Counter;
use itertools::Itertools;
use rgsl::randist::multinomial::multinomial_lnpdf;

use crate::variants::evidence::observation::{Observation, ReadPosition};

pub(crate) type Mapq = u32;
pub(crate) type ObservedMapqs = Counter<Mapq>;

#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub(crate) enum LocusType {
    Unique,
    Homologous,
}

#[derive(Clone, Debug)]
pub(crate) struct Prior;

impl bayesian::model::Prior for Prior {
    type Event = LocusType;

    fn compute(&self, _event: &Self::Event) -> bio::stats::LogProb {
        LogProb::ln_one() // flat prior
    }
}

#[derive(Clone, Debug)]
pub(crate) struct Posterior;

impl bayesian::model::Posterior for Posterior {
    type Event = LocusType;

    type BaseEvent = LocusType;

    type Data = ObservedMapqs;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        event: &Self::Event,
        data: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        joint_prob(event, data)
    }
}

#[derive(Clone, Debug, new)]
pub(crate) struct Likelihood {
    max_mapq: Mapq,
}

impl bayesian::model::Likelihood for Likelihood {
    type Event = LocusType;

    type Data = ObservedMapqs;

    fn compute(&self, event: &Self::Event, data: &Self::Data, _payload: &mut ()) -> LogProb {
        let obs = (0..self.max_mapq)
            .into_iter()
            .map(|mapq| data.get(&mapq).cloned().unwrap_or(0) as u32)
            .collect_vec();
        match event {
            LocusType::Unique => LogProb(multinomial_lnpdf(
                &(0..self.max_mapq).into_iter().map(|mapq| 0.0).collect_vec(),
                &obs,
            )),
            LocusType::Homologous => LogProb(multinomial_lnpdf(
                &(0..self.max_mapq + 1)
                    .into_iter()
                    .map(|mapq| 1.0 / self.max_mapq as f64)
                    .collect_vec(),
                &obs,
            )),
        }
    }
}

pub(crate) fn adjust_prob_mapping(pileup: &mut [Observation<ReadPosition>]) {
    if !pileup.is_empty() {
        // METHOD: we use a simple multinomial model to represent two cases:
        // (a) the locus is a homolog of something else: this means that MAPQs can be expected to be somehow uniformly distributed.
        // (b) the locus is unique: this means that MAPQs should be all at the maximum.
        // The rationale for the former case is that we always expect some noise of artificially high MAPQs because alleles can be unknown by the mapper
        // and sometimes base errors make a read randomly fit better than it should.

        let data: Counter<_> = pileup
            .iter()
            .map(|obs| PHREDProb::from(obs.prob_mapping_orig()).round() as u32)
            .collect();
        let max_mapq = *data.keys().max().unwrap();

        let model = bayesian::model::Model::new(Likelihood::new(max_mapq), Prior, Posterior);

        let m = model.compute(vec![LocusType::Unique, LocusType::Unique], &data);

        // Calculate expected value of MAPQ given the two posteriors and their repective MAPQs (0 and max_mapq).
        let adjusted =
            m.posterior(&LocusType::Unique).unwrap() + LogProb::from(PHREDProb(max_mapq as f64));
        for obs in pileup {
            if adjusted < obs.prob_mapping_orig() {
                obs.set_prob_mapping_adj(adjusted);
            }
        }
    }
}
