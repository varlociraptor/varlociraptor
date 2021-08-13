use std::cmp;

use bio::stats::bayesian::model::{Likelihood, Model, Posterior, Prior};
use bio::stats::LogProb;
use derive_builder::Builder;
use itertools::{Itertools, MinMaxResult};
use vec_map::VecMap;
use itertools_num::linspace;

use crate::grammar;
use crate::utils::PROB_05;
use crate::utils::adaptive_integration;
use crate::variants::model;
use crate::variants::model::likelihood;
use crate::variants::model::{bias::Biases, AlleleFreq, Contamination, VariantType};
use crate::variants::sample::Pileup;

#[derive(new, Clone, Debug)]
pub(crate) struct Snv {
    refbase: u8,
    altbase: u8,
}

#[derive(new, Clone, Debug, Getters)]
#[get = "pub"]
pub(crate) struct Data {
    pileups: Vec<Pileup>,
    snv: Option<Snv>,
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
    resolutions: Option<grammar::SampleInfo<grammar::Resolution>>,
    contaminations: Option<grammar::SampleInfo<Option<Contamination>>>,
    prior: P,
}

impl<P> GenericModelBuilder<P>
where
    P: Prior<Event = Vec<likelihood::Event>>,
{
    pub(crate) fn resolutions(mut self, resolutions: grammar::SampleInfo<grammar::Resolution>) -> Self {
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
    resolutions: grammar::SampleInfo<grammar::Resolution>,
}

impl GenericPosterior {
    fn grid_points(&self, pileups: &[Pileup]) -> Vec<grammar::Resolution> {
        pileups
            .iter()
            .zip(self.resolutions.iter())
            .map(|(pileup, res)| {
                if let grammar::Resolution::Uniform(res) = res {
                    let n_obs = pileup.len();
                    let mut n = cmp::min(cmp::max(n_obs + 1, 5), *res);
                    if n % 2 == 0 {
                        n += 1;
                    }
                    grammar::Resolution::Uniform(n)
                } else {
                    grammar::Resolution::Adaptive
                }
            })
            .collect()
    }

    fn density<F: FnMut(&<Self as Posterior>::BaseEvent, &<Self as Posterior>::Data) -> LogProb>(
        &self,
        vaf_tree_node: &grammar::vaftree::Node,
        base_events: &mut VecMap<likelihood::Event>,
        sample_grid_points: &[grammar::Resolution],
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
            grammar::vaftree::NodeKind::False => LogProb::ln_zero(),
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
                        let resolution = sample_grid_points[*sample];
                        let min_vaf = vafs.observable_min(n_obs);
                        let max_vaf = vafs.observable_max(n_obs);
                        let density = |vaf| {
                            let mut base_events = base_events.clone();
                            push_base_event(AlleleFreq(vaf), &mut base_events);
                            subdensity(&mut base_events)
                        };

                        if (max_vaf - min_vaf) < *resolution {
                            // METHOD: Interval too small for desired resolution.
                            // Just use 3 grid points.
                            LogProb::ln_simpsons_integrate_exp(
                                |_, vaf| density,
                                *min_vaf,
                                *max_vaf,
                                3,
                            )
                        } else if n_obs < 5 {
                            // METHOD: Not enough observations to expect a unimodal density.
                            // Use 11 grid points.
                            LogProb::ln_simpsons_integrate_exp(
                                |_, vaf| density,
                                *min_vaf,
                                *max_vaf,
                                11,
                            )
                        } else {
                            // METHOD: enough data and large enough interval, use adaptive integration
                            // at the desired resolution.
                            adaptive_integration::ln_integrate_exp(density, min_vaf, max_vaf, *resolution)
                        }
                        
                        match sample_grid_points[*sample] {
                            grammar::Resolution::Uniform(n) => {
                                LogProb::ln_simpsons_integrate_exp(
                                    |_, vaf| density,
                                    *min_vaf,
                                    *max_vaf,
                                    n,
                                )
                            }
                            grammar::Resolution::Adaptive if (max_vaf - min_vaf) > AlleleFreq(0.01)  => {
                                
                            }
                            grammar::Resolution::
                        }


                            // select number of points such that step size is 0.01, at least 5
                            let get_n = |min_vaf, max_vaf, step_size| {
                                let n = cmp::max(((max_vaf - min_vaf) / step_size).round() as usize, 5);
                                if n % 2 == 0 {
                                    n + 1
                                } else {
                                    n
                                }
                            };
                            let n = get_n(min_vaf, max_vaf, 0.01);

                            // generate grid
                            let grid_points: BTreeSet<_> = linspace(min_vaf, max_vaf, n).map(|vaf| AlleleFreq(vaf)).collect();
                            
                            // take 5 points for probing
                            let probes: BTreeMap<_, _> = grid_points.iter().map(|vaf| (AlleleFreq(vaf), density(vaf))).collect();

                            // get maximum probability
                            if let Some((vaf, max_prob)) = probes.iter().max_by_key(|(_, prob)| prob) {
                                // get standard deviation around max (binomial with n_obs, see paper)
                                let stddev = 1.0 / n_obs * (n_obs * *vaf * (1.0 - *vaf)).sqrt();
                                // get grid points enclosed by the stddev
                                let std_range: Vec<_> = grid_points.range(AlleleFreq(cmp::max(0.0, *vaf - stddev))..AlleleFreq(cmp::min(*vaf + stddev, 1.0))).collect();

                                let mut n = range.len();

                                let probed_density = |_, vaf| {
                                    // if already probed, reuse, otherwise calculate
                                    probes.get(NotNan(vaf)).unwrap_or_else(|| {
                                        density(vaf)
                                    })
                                };

                                // integrate over standard deviation interval
                                let std_prob = LogProb::ln_trapezoid_integrate_grid_exp(
                                    probed_density,
                                    &std_range
                                );

                                let range = grid_points.range(..std_range.first());
                                let left_tail_prob = LogProb::ln_trapezoid_integrate_grid_exp(
                                    probed_density,
                                    range.step_by(range.len() / 3)
                                );

                                let right_tail_prob = LogProb::ln_trapezoid_integrate_grid_exp(
                                    probed_density,
                                    grid_points.range(std_range.last()..)
                                );

                                left_tail_prob + std_prob + right_tail_prob
                            } else {
                                panic!("bug: probes map is supposed to be non empty")
                            }

                            


                        }
                    }
                }
            }
            grammar::vaftree::NodeKind::Variant {
                positive,
                refbase: given_refbase,
                altbase: given_altbase,
            } => {
                if let Some(Snv { refbase, altbase }) = data.snv {
                    let contains =
                        given_refbase.contains(refbase) && given_altbase.contains(altbase);
                    if (*positive && !contains) || (!*positive && contains) {
                        // abort computation, branch does not allow this variant
                        LogProb::ln_zero()
                    } else {
                        // skip this node
                        subdensity(base_events)
                    }
                } else if *positive {
                    // no SNV but branch requires the defined SNV, hence abort with prob 0
                    LogProb::ln_zero()
                } else {
                    // skip this node, as we don't have the defined SNV but it is negated
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
            *PROB_05 + LogProb((1.0 / event.biases.len() as f64).ln())
        } else {
            *PROB_05
        };

        // METHOD: filter out biases that are impossible to observe, (e.g. + without any + observation).
        let possible_biases = event.biases.iter().filter(|bias| {
            bias.is_possible(&data.pileups)
                && bias.is_informative(&data.pileups)
                && bias.is_likely(&data.pileups)
        });

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

// TODO: remove the following in favor of the new universal prior.

#[derive(Default, Clone, Debug)]
pub(crate) struct FlatPrior {
    universe: Option<grammar::SampleInfo<grammar::VAFUniverse>>,
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

impl model::prior::UpdatablePrior for FlatPrior {
    fn set_universe_and_ploidies(
        &mut self,
        universe: grammar::SampleInfo<grammar::VAFUniverse>,
        _ploidies: grammar::SampleInfo<Option<u32>>,
    ) {
        self.universe = Some(universe);
    }

    fn set_variant_type(&mut self, _: VariantType) {}
}
