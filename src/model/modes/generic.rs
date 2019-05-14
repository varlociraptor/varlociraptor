use std::cmp;
use std::rc::Rc;

use bio::stats::bayesian::model::{Likelihood, Model, Posterior, Prior};
use bio::stats::LogProb;
use vec_map::VecMap;
use itertools::Itertools;

use crate::model;
use crate::model::likelihood;
use crate::model::sample::Pileup;
use crate::model::{AlleleFreq, Contamination, StrandBias};
use crate::grammar;

#[derive(Debug)]
pub enum CacheEntry {
    ContaminatedSample(likelihood::ContaminatedSampleCache),
    SingleSample(likelihood::SingleSampleCache),
}

impl CacheEntry {
    fn new(contaminated: bool) -> Self {
        if contaminated {
            CacheEntry::ContaminatedSample(likelihood::ContaminatedSampleCache::default())
        } else {
            CacheEntry::SingleSample(likelihood::SingleSampleCache::default())
        }
    }
}

pub type Cache = VecMap<CacheEntry>;

#[derive(Default, Debug, Clone)]
pub struct GenericModelBuilder<P> {
    resolutions: Vec<usize>,
    contaminations: Vec<Option<Contamination>>,
    prior: P,
}

impl<P: Prior> GenericModelBuilder<P>
where
    P: Prior<Event = VecMap<likelihood::Event>>,
{
    pub fn push_sample(mut self, resolution: usize, contamination: Option<Contamination>) -> Self {
        self.contaminations.push(contamination);
        self.resolutions.push(resolution);

        self
    }

    pub fn prior(mut self, prior: P) -> Self {
        self.prior = prior;

        self
    }

    pub fn build(self) -> Result<Model<GenericLikelihood, P, GenericPosterior, Cache>, String> {
        let posterior = GenericPosterior::new(self.resolutions);
        let likelihood = GenericLikelihood::new(self.contaminations);
        Ok(Model::new(likelihood, self.prior, posterior))
    }
}

#[derive(new, Default, Clone, Debug)]
pub struct GenericPosterior {
    resolutions: Vec<usize>,
}

impl GenericPosterior {
    fn grid_points(&self, pileups: &[Pileup]) -> Vec<usize> {
        pileups
            .iter()
            .zip(self.resolutions.iter())
            .map(|(pileup, res)| {
                let n_obs = pileup.len();
                let mut n = cmp::min(cmp::max(n_obs + 1, 5), *res);
                if n % 2 == 0 {
                    n += 1;
                }
                n
            })
            .collect()
    }

    fn density<F: FnMut(&<Self as Posterior>::BaseEvent, &<Self as Posterior>::Data) -> LogProb>(
        &self,
        formula: &grammar::Formula<usize>,
        operand: Option<usize>,
        base_events: Rc<VecMap<likelihood::Event>>,
        sample_grid_points: &[usize],
        pileups: &<Self as Posterior>::Data,
        strand_bias: StrandBias,
        joint_prob: &mut F,
    ) -> LogProb {
        let (subformula, is_leaf) = match (formula, operand) {
            (grammar::Formula::Conjunction { operands }, Some(operand)) => {
                let is_last = operand == operands.len() - 1;
                (operands[operand].as_ref(), is_last)
            },
            (grammar::Formula::Atom { .. }, None) => (formula, true),
            _ => (formula, false),
        };

        match formula {
            grammar::Formula::Atom { sample, ref vafs } => {
                let mut push_base_event = |allele_freq, base_events: Rc<VecMap<likelihood::Event>>| {
                    base_events.insert(
                        *sample,
                        likelihood::Event {
                            allele_freq: allele_freq,
                            strand_bias: strand_bias,
                        }
                    );
                };

                let subdensity = |base_events: Rc<VecMap<likelihood::Event>>| {
                    if is_leaf {
                        joint_prob(base_events.as_ref(), pileups)
                    } else {
                        self.density(
                            formula,
                            operand.map(|o| o + 1),
                            Rc::clone(&base_events),
                            sample_grid_points,
                            pileups,
                            strand_bias,
                            joint_prob,
                        )
                    }
                };
                match vafs {
                    grammar::VAFSpectrum::Set(vafs) => {
                        if vafs.len() == 1 {
                            push_base_event(*vafs.iter().next().unwrap(), base_events);
                            subdensity(base_events)
                        } else {
                            LogProb::ln_sum_exp(&vafs.iter().cloned().map(|vaf| {
                                let base_events = base_events.clone();
                                push_base_event(vaf, base_events);
                                subdensity(base_events)
                            }).collect_vec())
                        }
                    }
                    grammar::VAFSpectrum::Range(range) => {
                        let n_obs = pileups[*sample].len();
                        LogProb::ln_simpsons_integrate_exp(
                            |_, vaf| {
                                let base_events = base_events.clone();
                                push_base_event(AlleleFreq(vaf), base_events);
                                subdensity(base_events)
                            },
                            *range.observable_min(n_obs),
                            *range.observable_max(n_obs),
                            sample_grid_points[*sample],
                        )
                    }
                }
            }
            grammar::Formula::Conjunction { operands } => {
                self.density(
                    formula,
                    Some(operand.map_or(0, |o| o + 1)),
                    base_events,
                    sample_grid_points,
                    pileups,
                    strand_bias,
                    joint_prob,
                )
            }
            grammar::Formula::Disjunction { operands } => {
                LogProb::ln_sum_exp(&operands.iter().map(|o| {
                    self.density(
                        o.as_ref(),
                        None,
                        base_events,
                        sample_grid_points,
                        pileups,
                        strand_bias,
                        joint_prob,
                    )
                }).collect_vec())
            }
        }
    }
}

impl Posterior for GenericPosterior {
    type BaseEvent = VecMap<likelihood::Event>;
    type Event = model::Event;
    type Data = Vec<Pileup>;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        event: &Self::Event,
        pileups: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        let grid_points = self.grid_points(pileups);

        self.density(&event.formula, None, Rc::new(VecMap::new()), &grid_points, pileups, event.strand_bias, joint_prob)
    }
}

#[derive(Clone, Debug)]
enum SampleModel {
    Contaminated {
        likelihood_model: likelihood::ContaminatedSampleLikelihoodModel,
        by: usize,
    },
    Normal(likelihood::SampleLikelihoodModel),
}

#[derive(Default, Clone, Debug)]
pub struct GenericLikelihood {
    inner: Vec<SampleModel>,
}

impl GenericLikelihood {
    pub fn new(contaminations: Vec<Option<Contamination>>) -> Self {
        let mut inner = Vec::new();
        for contamination in contaminations.iter() {
            if let Some(contamination) = contamination {
                inner.push(SampleModel::Contaminated {
                    likelihood_model: likelihood::ContaminatedSampleLikelihoodModel::new(
                        1.0 - contamination.fraction,
                    ),
                    by: contamination.by,
                });
            } else {
                inner.push(SampleModel::Normal(likelihood::SampleLikelihoodModel::new()));
            }
        }
        GenericLikelihood { inner }
    }
}

impl Likelihood<Cache> for GenericLikelihood {
    type Event = VecMap<likelihood::Event>;
    type Data = Vec<Pileup>;

    fn compute(&self, events: &Self::Event, pileups: &Self::Data, cache: &mut Cache) -> LogProb {
        let mut p = LogProb::ln_one();

        for (sample, event) in events {
            let pileup = &pileups[sample];
            let inner = &self.inner[sample];
            p += match inner {
                &SampleModel::Contaminated {
                    ref likelihood_model,
                    by,
                } => {
                    if let CacheEntry::ContaminatedSample(ref mut cache) =
                        cache.entry(sample).or_insert_with(|| CacheEntry::new(true))
                    {
                        likelihood_model.compute(
                            &likelihood::ContaminatedSampleEvent {
                                primary: event.clone(),
                                // TODO ensure that VAF is always defined for contaminating sample
                                secondary: events.get(by).expect("bug: no allele frequency given for contaminating sample").clone(),
                            },
                            pileup,
                            cache,
                        )
                    } else {
                        unreachable!();
                    }
                }
                &SampleModel::Normal(ref likelihood_model) => {
                    if let CacheEntry::SingleSample(ref mut cache) = cache
                        .entry(sample)
                        .or_insert_with(|| CacheEntry::new(false))
                    {
                        likelihood_model.compute(event, pileup, cache)
                    } else {
                        unreachable!();
                    }
                }
            }
        }

        p
    }
}

#[derive(Default, Clone, Debug)]
pub struct FlatPrior {}

impl FlatPrior {
    pub fn new() -> Self {
        FlatPrior {}
    }
}

impl Prior for FlatPrior {
    type Event = VecMap<likelihood::Event>;

    fn compute(&self, _event: &Self::Event) -> LogProb {
        LogProb::ln_one()
    }
}
