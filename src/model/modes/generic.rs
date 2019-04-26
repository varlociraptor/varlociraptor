use std::cmp;

use bio::stats::bayesian::model::{Model, Likelihood, Posterior, Prior};
use bio::stats::LogProb;
use petgraph::{Graph, graph::NodeIndex, algo::kosaraju_scc};
use itertools::Itertools;
use derive_builder::Builder;
use vec_map::VecMap;

use crate::model;
use crate::model::likelihood;
use crate::model::sample::Pileup;
use crate::model::{AlleleFreq, ContinuousAlleleFreqs, Contamination};
use crate::utils;

pub struct GenericModelBuilder<P: Prior> {
    resolutions: Vec<usize>,
    contaminations: Vec<Option<Contamination>>,
    prior: P
}

impl<P: Prior> GenericModelBuilder<P> {
    pub fn push_sample(&mut self, resolution: usize, contamination: Option<Contamination>) -> &mut Self {
        self.contaminations.push(contamination);
        self.resolutions.push(resolution);

        self
    }

    pub fn prior(&mut self, prior: P) -> &mut Self {
        self.prior = prior;

        self
    }

    pub fn build(self) -> Model<GenericLikelihood, P, GenericPosterior> {
        let posterior = GenericPosteriorBuilder::default().contaminations(&self.contaminations).resolutions(self.resolutions).build();
        let likelihood = GenericLikelihood::new(self.contaminations);

        Model::new(likelihood, prior, posterior)
    }
}

#[derive(Builder)]
pub struct GenericPosterior {
    resolutions: Vec<usize>,
    /// Connected components of the contamination graph.
    /// These can be called independently.
    #[builder(private)]
    groups: Vec<Vec<usize>>,
}

impl GenericPosteriorBuilder {
    pub fn contaminations(&mut self, contaminations: &[Option<Contamination>]) -> &mut Self {
        let mut contamination_graph = Graph::new_undirected();
        for sample in 0..contaminations.len() {
            contamination_graph.add_node(());
        }
        for (sample, contamination) in contaminations.iter().enumerate() {
            if let Some(contamination) = contamination {
                contamination_graph.add_edge(NodeIndex::new(sample), NodeIndex::new(contamination.by), ());
            }
        }

        self.groups(kosaraju_scc(&contamination_graph).into_iter().map(|group| {
            group.into_iter().map(|i| i.index()).collect_vec()
        }).collect_vec())
    }
}

impl GenericPosterior {
    fn grid_points(&self, pileups: &[Pileup]) -> Vec<usize> {
        pileups.iter().zip(self.resolutions.iter()).map(|(pileup, res)| {
            let n_obs = pileup.len();
            let mut n = cmp::min(cmp::max(n_obs + 1, 5), *res);
            if n % 2 == 0 {
                n += 1;
            }
            n
        }).collect()
    }
}

impl Posterior for GenericPosterior {
    type BaseEvent = GroupEvent;
    type Event = Vec<model::Event<ContinuousAlleleFreqs>>;
    type Data = Vec<Pileup>;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        event: &Self::Event,
        pileups: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        let grid_points = self.grid_points(pileups);

        let mut p = LogProb::ln_one();
        for group in &self.groups {
            if group.len() == 1 {
                // Simple case, just one sample that is independent of the rest.

                let s = group[0];
                let sample_event = &event[s];
                if sample_event.allele_freqs.is_singleton() {
                    let mut group_event = GroupEvent::default();
                    group_event.push(s, likelihood::Event {
                        allele_freq: sample_event.allele_freqs.start,
                        strand_bias: sample_event.strand_bias,
                    });
                    p += joint_prob(&group_event, pileups);
                } else {
                    let n_obs = pileups[s].len();
                    p += LogProb::ln_simpsons_integrate_exp(
                        |_, allele_freq| {
                            let mut group_event = GroupEvent::default();
                            group_event.push(s, likelihood::Event {
                                allele_freq: AlleleFreq(allele_freq),
                                strand_bias: sample_event.strand_bias,
                            });
                            joint_prob(
                                &group_event,
                                pileups
                            )
                        },
                        *sample_event.allele_freqs.observable_min(n_obs),
                        *sample_event.allele_freqs.observable_max(n_obs),
                        grid_points[s],
                    );
                }
            } else {
                // Multiple samples contaminating each other.

                fn density<F: FnMut(&GroupEvent, &Vec<Pileup>) -> LogProb>(s: usize, group: &[usize], sample_events: &[model::Event<ContinuousAlleleFreqs>], sample_grid_points: &[usize], pileups: &Vec<Pileup>, base_events: Vec<likelihood::Event>, joint_prob: &mut F) -> LogProb {
                    let e = &sample_events[s];
                    let mut subdensity = |allele_freq| {
                        let mut base_events = base_events.clone();
                        base_events.push(likelihood::Event {
                            allele_freq: allele_freq,
                            strand_bias: e.strand_bias,
                        });

                        if s == sample_events.len() - 1 {
                            let mut group_event = GroupEvent::default();
                            for (s, base_event) in group.iter().zip(base_events.into_iter()) {
                                group_event.push(*s, base_event);
                            }

                            joint_prob(&group_event, pileups)
                        } else {
                            density(s + 1, group, sample_events, sample_grid_points, pileups, base_events, joint_prob)
                        }
                    };
                    if e.allele_freqs.is_singleton() {
                        subdensity(e.allele_freqs.start)
                    } else {
                        let n_obs = pileups[s].len();
                        LogProb::ln_simpsons_integrate_exp(
                            |_, allele_freq| {
                                subdensity(AlleleFreq(allele_freq))
                            },
                            *e.allele_freqs.observable_min(n_obs),
                            *e.allele_freqs.observable_max(n_obs),
                            sample_grid_points[s],
                        )
                    }
                }

                if group.len() == pileups.len() {
                    p += density(0, &group, event, &grid_points, pileups, Vec::new(), joint_prob);
                } else {
                    let group_events = utils::select(group, event).collect_vec();
                    let group_grid_points = utils::select(group, &grid_points).collect_vec();

                    p += density(0, &group, &group_events, &group_grid_points, pileups, Vec::new(), joint_prob);
                }
            }
        }
        p
    }
}

#[derive(Default)]
pub struct GroupEvent {
    events: VecMap<likelihood::Event>,
    group: Vec<usize>
}

impl GroupEvent {
    pub fn push(&mut self, sample: usize, event: likelihood::Event) {
        self.events.insert(sample, event);
        self.group.push(sample);
    }

    pub fn event(&self, sample: usize) -> &likelihood::Event {
        self.events.get(sample).unwrap()
    }
}


pub struct GenericLikelihood {
    contaminations: Vec<Option<usize>>,
    contaminated: VecMap<likelihood::ContaminatedSampleLikelihoodModel>,
    normal: VecMap<likelihood::SampleLikelihoodModel>,
}

impl GenericLikelihood {
    pub fn new(contaminations: Vec<Option<Contamination>>) -> Self {
        let mut contaminated = VecMap::new();
        let mut normal = VecMap::new();
        let mut _contaminations = Vec::new();
        for (sample, contamination) in contaminations.iter().enumerate() {
            _contaminations.push(
                if let Some(contamination) = contamination {
                    contaminated.insert(sample, likelihood::ContaminatedSampleLikelihoodModel::new(1.0 - contamination.fraction));
                    Some(contamination.by)
                } else {
                    normal.insert(sample, likelihood::SampleLikelihoodModel::new());
                    None
                }
            );
        }
        GenericLikelihood {
            contaminations: _contaminations,
            contaminated,
            normal,
        }
    }
}

impl Likelihood for GenericLikelihood {
    type Event = GroupEvent;
    type Data = Vec<Pileup>;

    fn compute(&self, event: &Self::Event, pileups: &Self::Data) -> LogProb {
        let mut p = LogProb::ln_one();
        for sample in event.group.iter() {
            let pileup = &pileups[*sample];
            p += if let Some(contaminant) = self.contaminations[*sample] {
                let event = vec![event.event(*sample).clone(), event.event(contaminant).clone()];
                self.contaminated.get(*sample).unwrap().compute(&event, pileup.as_ref())
            } else {
                self.normal.get(*sample).unwrap().compute(&event.events[*sample], pileup.as_ref())
            };
        }

        p
    }
}
