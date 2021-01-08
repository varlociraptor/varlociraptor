use std::cmp;

use bio::stats::bayesian::model::{Likelihood, Model, Posterior, Prior};
use bio::stats::LogProb;
use derive_builder::Builder;
use itertools::Itertools;
use vec_map::VecMap;

use crate::grammar;
use crate::utils::PROB_HALF;
use crate::variants::model;
use crate::variants::model::likelihood;
use crate::variants::model::{bias::Biases, AlleleFreq, Contamination};
use crate::variants::sample::Pileup;

#[derive(new, Clone, Debug)]
pub(crate) struct SNV {
    refbase: u8,
    altbase: u8,
}

#[derive(new, Clone, Debug, Getters)]
#[get = "pub"]
pub(crate) struct Data {
    pileups: Vec<Pileup>,
    snv: Option<SNV>,
}

impl Data {
    pub(crate) fn into_pileups(self) -> Vec<Pileup> {
        self.pileups
    }
}

#[derive(Debug)]
pub(crate) enum CacheEntry {
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

pub(crate) type Cache = VecMap<CacheEntry>;

#[derive(Default, Debug, Clone, Builder)]
pub(crate) struct GenericModelBuilder<P>
where
    P: Prior<Event = Vec<likelihood::Event>>,
{
    resolutions: Option<grammar::SampleInfo<usize>>,
    contaminations: Option<grammar::SampleInfo<Option<Contamination>>>,
    prior: P,
}

impl<P> GenericModelBuilder<P>
where
    P: Prior<Event = Vec<likelihood::Event>>,
{
    pub(crate) fn resolutions(mut self, resolutions: grammar::SampleInfo<usize>) -> Self {
        self.resolutions = Some(resolutions);

        self
    }

    pub(crate) fn contaminations(
        mut self,
        contaminations: grammar::SampleInfo<Option<Contamination>>,
    ) -> Self {
        self.contaminations = Some(contaminations);

        self
    }

    pub(crate) fn prior(mut self, prior: P) -> Self {
        self.prior = prior;

        self
    }

    pub(crate) fn build(
        self,
    ) -> Result<Model<GenericLikelihood, P, GenericPosterior, Cache>, String> {
        let posterior = GenericPosterior::new(
            self.resolutions
                .expect("GenericModelBuilder: need to call resolutions() before build()"),
        );
        let likelihood = GenericLikelihood::new(
            self.contaminations
                .expect("GenericModelBuilder: need to call contaminations() before build()"),
        );
        Ok(Model::new(likelihood, self.prior, posterior))
    }
}

#[derive(new, Clone, Debug, Default)]
pub(crate) struct GenericPosterior {
    resolutions: grammar::SampleInfo<usize>,
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
        vaf_tree_node: &grammar::vaftree::Node,
        base_events: &mut VecMap<likelihood::Event>,
        sample_grid_points: &[usize],
        data: &<Self as Posterior>::Data,
        biases: &Biases,
        joint_prob: &mut F,
    ) -> LogProb {
        let mut subdensity = |base_events: &mut VecMap<likelihood::Event>| {
            if vaf_tree_node.is_leaf() {
                joint_prob(&base_events.values().cloned().collect(), data)
            } else if vaf_tree_node.is_branching() {
                LogProb::ln_sum_exp(
                    &vaf_tree_node
                        .children()
                        .iter()
                        .map(|child| {
                            self.density(
                                child,
                                &mut base_events.clone(),
                                sample_grid_points,
                                data,
                                biases,
                                joint_prob,
                            )
                        })
                        .collect_vec(),
                )
            } else {
                self.density(
                    &vaf_tree_node.children()[0],
                    base_events,
                    sample_grid_points,
                    data,
                    biases,
                    joint_prob,
                )
            }
        };

        match vaf_tree_node.kind() {
            grammar::vaftree::NodeKind::Sample { sample, vafs } => {
                let push_base_event = |allele_freq, base_events: &mut VecMap<likelihood::Event>| {
                    base_events.insert(
                        *sample,
                        likelihood::Event {
                            allele_freq,
                            biases: biases.clone(),
                        },
                    );
                };

                match vafs {
                    grammar::VAFSpectrum::Set(vafs) => {
                        if vafs.len() == 1 {
                            push_base_event(*vafs.iter().next().unwrap(), base_events);
                            subdensity(base_events)
                        } else {
                            LogProb::ln_sum_exp(
                                &vafs
                                    .iter()
                                    .map(|vaf| {
                                        let mut base_events = base_events.clone();
                                        push_base_event(*vaf, &mut base_events);
                                        subdensity(&mut base_events)
                                    })
                                    .collect_vec(),
                            )
                        }
                    }
                    grammar::VAFSpectrum::Range(vafs) => {
                        let n_obs = data.pileups[*sample].len();
                        LogProb::ln_simpsons_integrate_exp(
                            |_, vaf| {
                                let mut base_events = base_events.clone();
                                push_base_event(AlleleFreq(vaf), &mut base_events);
                                subdensity(&mut base_events)
                            },
                            *vafs.observable_min(n_obs),
                            *vafs.observable_max(n_obs),
                            sample_grid_points[*sample],
                        )
                    }
                }
            }
            grammar::vaftree::NodeKind::Variant {
                positive,
                refbase: given_refbase,
                altbase: given_altbase,
            } => {
                if let Some(SNV { refbase, altbase }) = data.snv {
                    let contains =
                        given_refbase.contains(refbase) && given_altbase.contains(altbase);
                    if (*positive && !contains) || (!*positive && contains) {
                        // abort computation, branch does not allow this variant
                        LogProb::ln_zero()
                    } else {
                        // skip this node
                        subdensity(base_events)
                    }
                } else {
                    // skip this node
                    subdensity(base_events)
                }
            }
        }
    }
}

impl Posterior for GenericPosterior {
    type BaseEvent = Vec<likelihood::Event>;
    type Event = model::Event;
    type Data = Data;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        event: &Self::Event,
        data: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        let grid_points = self.grid_points(&data.pileups);
        let vaf_tree = &event.vafs;
        let bias_prior = if event.is_artifact() {
            *PROB_HALF + LogProb((1.0 / event.biases.len() as f64).ln())
        } else {
            *PROB_HALF
        };

        // METHOD: filter out biases that are impossible to observe, (e.g. + without any + observation).
        let possible_biases = event
            .biases
            .iter()
            .filter(|bias| bias.is_possible(&data.pileups));
        LogProb::ln_sum_exp(
            &possible_biases
                .cartesian_product(vaf_tree)
                .map(|(biases, node)| {
                    let mut base_events = VecMap::with_capacity(data.pileups.len());
                    bias_prior
                        + self.density(
                            node,
                            &mut base_events,
                            &grid_points,
                            data,
                            biases,
                            joint_prob,
                        )
                })
                .collect_vec(),
        )
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

#[derive(Clone, Debug, Default)]
pub(crate) struct GenericLikelihood {
    inner: grammar::SampleInfo<SampleModel>,
}

impl GenericLikelihood {
    pub(crate) fn new(contaminations: grammar::SampleInfo<Option<Contamination>>) -> Self {
        let inner = contaminations.map(|contamination| {
            if let Some(contamination) = contamination {
                SampleModel::Contaminated {
                    likelihood_model: likelihood::ContaminatedSampleLikelihoodModel::new(
                        1.0 - contamination.fraction,
                    ),
                    by: contamination.by,
                }
            } else {
                SampleModel::Normal(likelihood::SampleLikelihoodModel::new())
            }
        });

        GenericLikelihood { inner }
    }
}

impl Likelihood<Cache> for GenericLikelihood {
    type Event = Vec<likelihood::Event>;
    type Data = Data;

    fn compute(&self, events: &Self::Event, data: &Self::Data, cache: &mut Cache) -> LogProb {
        let mut p = LogProb::ln_one();

        for (((sample, event), pileup), inner) in events
            .iter()
            .enumerate()
            .zip(data.pileups.iter())
            .zip(self.inner.iter())
        {
            p += match *inner {
                SampleModel::Contaminated {
                    ref likelihood_model,
                    by,
                } => {
                    if let CacheEntry::ContaminatedSample(ref mut cache) =
                        cache.entry(sample).or_insert_with(|| CacheEntry::new(true))
                    {
                        likelihood_model.compute(
                            &likelihood::ContaminatedSampleEvent {
                                primary: event.clone(),
                                secondary: events[by].clone(),
                            },
                            pileup,
                            cache,
                        )
                    } else {
                        unreachable!();
                    }
                }
                SampleModel::Normal(ref likelihood_model) => {
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
pub(crate) struct FlatPrior {
    universe: Option<grammar::SampleInfo<grammar::VAFUniverse>>,
}

impl FlatPrior {
    pub(crate) fn new() -> Self {
        FlatPrior { universe: None }
    }
}

impl Prior for FlatPrior {
    type Event = Vec<likelihood::Event>;

    fn compute(&self, event: &Self::Event) -> LogProb {
        if event
            .iter()
            .zip(self.universe.as_ref().unwrap().iter())
            .any(|(e, u)| !u.contains(e.allele_freq))
        {
            // if any of the events is not allowed in the universe of the corresponding sample, return probability zero.
            return LogProb::ln_zero();
        }
        LogProb::ln_one()
    }
}

impl model::modes::UniverseDrivenPrior for FlatPrior {
    fn set_universe(&mut self, universe: grammar::SampleInfo<grammar::VAFUniverse>) {
        self.universe = Some(universe);
    }
}
