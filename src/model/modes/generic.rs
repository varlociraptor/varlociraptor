use std::cmp;

use bio::stats::bayesian::model::{Likelihood, Model, Posterior, Prior};
use bio::stats::LogProb;
use vec_map::VecMap;

use crate::model;
use crate::model::likelihood;
use crate::model::sample::Pileup;
use crate::model::{AlleleFreq, Contamination, ContinuousAlleleFreqs};
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
    P: Prior<Event = Vec<likelihood::Event>>,
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
        dbg!(&posterior);
        dbg!(&likelihood);
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
        i: usize,
        clause: &grammar::Clause,
        base_events: VecMap<likelihood::Event>,
        sample_grid_points: &[usize],
        pileups: &<Self as Posterior>::Data,
        joint_prob: &mut F,
        strand_bias: model::StrandBias,
    ) -> LogProb {
        let atom = &clause[i];
        let mut subdensity = |allele_freq| {
            let mut base_events = base_events.clone();
            base_events.insert(
                atom.sample,
                likelihood::Event {
                    allele_freq: allele_freq,
                    strand_bias: strand_bias,
                }
            );

            if i == clause.len() - 1 {
                joint_prob(&base_events, pileups)
            } else {
                self.density(
                    i + 1,
                    base_events,
                    sample_grid_points,
                    pileups,
                    joint_prob,
                )
            }
        };

        match atom.vafs {
            grammar::VAFSpectrum::Set(vafs) => {
                if vafs.len() == 1 {
                    subdensity(vafs.iter().next().unwrap())
                } else {
                    LogProb::ln_sum_exp(&vafs.iter().cloned().map(|vaf| subdensity(vaf)).collect_vec())
                }
            }
            grammar::VAFSpectrum::Range(range) => {
                let n_obs = pileups[atom.sample].len();
                LogProb::ln_simpsons_integrate_exp(
                    |_, allele_freq| subdensity(AlleleFreq(allele_freq)),
                    range.observable_min(n_obs),
                    range.observable_max(n_obs),
                    sample_grid_points[sample],
                )
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

        self.density(0, Vec::new(), &grid_points, event, pileups, joint_prob)
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
                                // TODO check this before!!
                                secondary: events.get(by).expect("bug: contamination sample VAF not defined in event clause").clone(),
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
    type Event = Vec<likelihood::Event>;

    fn compute(&self, _event: &Self::Event) -> LogProb {
        LogProb::ln_one()
    }
}
