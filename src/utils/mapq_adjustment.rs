use bio::stats::{bayesian, LogProb, PHREDProb};
use counter::Counter;
use itertools::Itertools;
use rgsl::randist::multinomial::multinomial_lnpdf;

use crate::variants::evidence::observation::{Observation, ReadPosition};

pub(crate) type Mapq = u32;
pub(crate) type ObservedMapqs = Vec<Mapq>;

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

#[derive(Clone, Debug)]
pub(crate) struct Likelihood {
    p_unique: Vec<f64>,
    p_homologous: Vec<f64>,
}

impl Likelihood {
    pub(crate) fn new(max_mapq: Mapq) -> Self {
        let mut p_unique = (0..max_mapq).into_iter().map(|mapq| 0.0).collect_vec();
        p_unique.push(1.0);
        Likelihood {
            p_unique,
            p_homologous: (0..=max_mapq)
                .into_iter()
                .map(|mapq| 1.0 / max_mapq as f64)
                .collect_vec(),
        }
    }
}

impl bayesian::model::Likelihood for Likelihood {
    type Event = LocusType;

    type Data = ObservedMapqs;

    fn compute(&self, event: &Self::Event, data: &Self::Data, _payload: &mut ()) -> LogProb {
        let p = match event {
            LocusType::Unique => LogProb(multinomial_lnpdf(&self.p_unique, data)),
            LocusType::Homologous => LogProb(multinomial_lnpdf(&self.p_homologous, data)),
        };
        dbg!((data, p, &self.p_unique, &self.p_homologous));
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

        let mapqs: Counter<_> = pileup
            .iter()
            .map(|obs| {PHREDProb::from(obs.prob_mapping_orig().ln_one_minus_exp()).round() as u32})
            .collect();
        let max_mapq = *mapqs.keys().max().unwrap();
        dbg!(max_mapq);
        let data = (0..=max_mapq)
            .into_iter()
            .map(|mapq| mapqs.get(&mapq).cloned().unwrap_or(0) as u32)
            .collect_vec();

        let model = bayesian::model::Model::new(Likelihood::new(max_mapq), Prior, Posterior);

        let m = model.compute(vec![LocusType::Unique, LocusType::Unique], &data);

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
