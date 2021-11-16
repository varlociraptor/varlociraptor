use bio::stats::{bayesian, LogProb, PHREDProb, Prob};
use counter::Counter;
use itertools::Itertools;
use ordered_float::NotNan;
use rgsl::randist::{binomial::binomial_pdf, multinomial::multinomial_lnpdf};
use statrs::function::{beta::ln_beta, factorial::ln_binomial};

use crate::variants::{
    evidence::observation::{Observation, ReadPosition},
    model::AlleleFreq,
};

use super::adaptive_integration;

pub(crate) type Mapq = u64;
pub(crate) type ObservedMapqs = Counter<Mapq>;

#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub(crate) struct Shape {
    alpha: NotNan<f64>,
    beta: NotNan<f64>,
}

pub(crate) fn beta_binomial_pmf(k: u64, n: u64, alpha: f64, beta: f64) -> LogProb {
    LogProb(
        ln_binomial(n, k) + ln_beta(k as f64 + alpha, (n - k) as f64 + beta) - ln_beta(alpha, beta),
    )
}

#[derive(Clone, Debug)]
pub(crate) struct Prior;

impl bayesian::model::Prior for Prior {
    type Event = Shape;

    fn compute(&self, _event: &Self::Event) -> bio::stats::LogProb {
        LogProb::ln_one() // flat prior
    }
}

#[derive(Clone, Debug)]
pub(crate) struct Posterior {
    max_mapq: Mapq,
}

impl bayesian::model::Posterior for Posterior {
    type Event = Mapq;

    type BaseEvent = Shape;

    type Data = ObservedMapqs;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        event: &Self::Event,
        data: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        adaptive_integration::ln_integrate_exp(
            |alpha| {
                // TODO what about event == 0?? this leads to a division by zero here. I guess we need an offset of 1 in the distribution
                let beta = NotNan::new(
                    (self.max_mapq as f64 * *alpha - *event as f64 * *alpha) / *event as f64,
                )
                .unwrap();
                joint_prob(&Shape { alpha, beta }, data)
            },
            NotNan::new(0.01).unwrap(),
            NotNan::new(200.0).unwrap(),
            NotNan::new(0.01).unwrap(),
        )
    }
}

#[derive(Clone, Debug)]
pub(crate) struct Likelihood {
    max_mapq: u64,
}

impl bayesian::model::Likelihood for Likelihood {
    type Event = Shape;

    type Data = ObservedMapqs;

    fn compute(&self, event: &Self::Event, data: &Self::Data, _payload: &mut ()) -> LogProb {
        LogProb(
            data.iter()
                .map(|(mapq, count)| {
                    *beta_binomial_pmf(*mapq, self.max_mapq, *event.alpha, *event.beta)
                        * (*count as f64).ln()
                })
                .sum(),
        )
    }
}

pub(crate) fn adjust_prob_mapping(pileup: &mut [Observation<ReadPosition>], max_mapq: u64) {
    if !pileup.is_empty() {
        // METHOD: we use a simple multinomial model to represent two cases:
        // (a) the locus is a homolog of something else: this means that MAPQs can be expected to be somehow uniformly distributed.
        // (b) the locus is unique: this means that MAPQs should be all at the maximum.
        // The rationale for the former case is that we always expect some noise of artificially high MAPQs because alleles can be unknown by the mapper
        // and sometimes base errors make a read randomly fit better than it should.

        let mapqs: Counter<_> = pileup
            .iter()
            .map(|obs| PHREDProb::from(obs.prob_mapping_orig().ln_one_minus_exp()).round() as u64)
            .collect();
        let adjusted = calc_adjusted(&mapqs, max_mapq);
        for obs in pileup {
            if adjusted < obs.prob_mapping_orig() {
                obs.set_prob_mapping_adj(adjusted);
            }
        }
    }
}

fn calc_adjusted(mapqs: &Counter<Mapq>, max_mapq: u64) -> LogProb {
    let model = bayesian::model::Model::new(
        Likelihood {
            max_mapq: max_mapq + 1,
        },
        Prior,
        Posterior {
            max_mapq: max_mapq + 1,
        },
    );

    let m = model.compute(1..=max_mapq, mapqs);

    dbg!((1..=max_mapq)
        .into_iter()
        .map(|mapq| { m.posterior(&mapq).unwrap().exp() })
        .collect_vec());

    let exp_values = (1..=max_mapq)
        .into_iter()
        .map(|mapq| m.posterior(&mapq).unwrap().exp() * (mapq as f64))
        .collect_vec();
    dbg!(&exp_values);
    dbg!(exp_values.iter().cloned().sum::<f64>());

    let adjusted = (1..=max_mapq)
        .into_iter()
        .map(|mapq| m.posterior(&mapq).unwrap().exp() * (mapq as f64))
        .sum();
    let adjusted = LogProb::from(PHREDProb(adjusted)).ln_one_minus_exp();

    dbg!(adjusted);
    adjusted
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_calc_adjusted() {
        let mut mapqs = Counter::default();
        //mapqs.insert(60, 4);
        mapqs.insert(0, 4);
        //mapqs.insert(28, 0);

        calc_adjusted(&mapqs, 60);
    }
}
