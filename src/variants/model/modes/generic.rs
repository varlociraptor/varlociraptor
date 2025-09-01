use bio::stats::bayesian::model::{Likelihood, Model, Posterior, Prior};
use bio::stats::LogProb;
use derive_builder::Builder;
use itertools::Itertools;
use vec_map::{Values, VecMap};

use crate::grammar::{self, VAFRange};
use crate::utils::adaptive_integration;
use crate::utils::log2_fold_change::{Log2FoldChange, Log2FoldChangePredicate};
use crate::utils::PROB_05;
use crate::variants::evidence::observations::pileup::Pileup;
use crate::variants::model;
use crate::variants::model::likelihood;
use crate::variants::model::likelihood::Event;
use crate::variants::model::{bias::Artifacts, AlleleFreq, Contamination};
use std::ops::Index;

#[derive(new, Clone, Debug)]
pub(crate) struct Snv {
    refbase: u8,
    altbase: u8,
}

#[derive(new, Debug, Getters)]
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
            CacheEntry::ContaminatedSample(likelihood::ContaminatedSampleCache::new(10000))
        } else {
            CacheEntry::SingleSample(likelihood::SingleSampleCache::new(10000))
        }
    }
}

pub(crate) type Cache = VecMap<CacheEntry>;

#[derive(Default, Debug, Clone, Builder)]
pub(crate) struct GenericModelBuilder<P>
where
    P: Prior<Event = LikelihoodOperands>,
{
    resolutions: Option<grammar::SampleInfo<grammar::Resolution>>,
    contaminations: Option<grammar::SampleInfo<Option<Contamination>>>,
    prior: P,
}

impl<P> GenericModelBuilder<P>
where
    P: Prior<Event = LikelihoodOperands>,
{
    pub(crate) fn resolutions(
        mut self,
        resolutions: grammar::SampleInfo<grammar::Resolution>,
    ) -> Self {
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

#[derive(Debug, Clone, Hash, PartialEq, Eq, Getters)]
#[getset(get = "pub(crate)")]
pub(crate) struct VafLfc {
    sample_a: usize,
    sample_b: usize,
    predicate: Log2FoldChangePredicate,
}

#[derive(Debug, Default, Clone, Hash, PartialEq, Eq, Getters)]
#[getset(get = "pub(crate)")]
pub(crate) struct LikelihoodOperands {
    events: VecMap<likelihood::Event>,
    lfcs: Vec<VafLfc>,
}

impl LikelihoodOperands {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    pub(crate) fn len(&self) -> usize {
        self.events.len()
    }

    pub(crate) fn push(&mut self, event: Event) {
        self.events.insert(self.events.len(), event);
    }

    pub(crate) fn iter(&self) -> Values<Event> {
        self.events.values()
    }

    pub(crate) fn is_discrete(&self) -> bool {
        self.events.values().all(|evt| evt.is_discrete)
    }

    pub(crate) fn is_absent(&self) -> bool {
        self.events.values().all(|evt| evt.is_absent())
    }

    pub(crate) fn lfc_bounds(&self, sample: usize) -> Option<VAFRange> {
        let mut acc: Option<VAFRange> = None;
        for lfc in &self.lfcs {
            let maybe_bounds = if lfc.sample_a == sample {
                // sample is on side A: invert the predicate against B’s allele_freq
                self.events
                    .get(lfc.sample_b)
                    .map(|evt_b| lfc.predicate.invert().infer_vaf_bounds(evt_b.allele_freq))
            } else if lfc.sample_b == sample {
                // sample is on side B: use the predicate directly against A’s allele_freq
                self.events
                    .get(lfc.sample_a)
                    .map(|evt_a| lfc.predicate.infer_vaf_bounds(evt_a.allele_freq))
            } else {
                // this LFC doesn’t involve the current sample
                None
            };

            if let Some(bounds) = maybe_bounds {
                acc = Some(match acc {
                    Some(prev) => prev.intersect(&bounds),
                    None => bounds,
                });
            }
        }
        acc
    }
}

impl Index<usize> for LikelihoodOperands {
    type Output = likelihood::Event;

    fn index(&self, index: usize) -> &Self::Output {
        &self.events[index]
    }
}

#[derive(new, Clone, Debug, Default)]
pub(crate) struct GenericPosterior {
    resolutions: grammar::SampleInfo<grammar::Resolution>,
}

impl GenericPosterior {
    fn density<F: FnMut(&<Self as Posterior>::BaseEvent, &<Self as Posterior>::Data) -> LogProb>(
        &self,
        vaf_tree_node: &grammar::vaftree::Node,
        likelihood_operands: &mut LikelihoodOperands,
        data: &<Self as Posterior>::Data,
        biases: &Artifacts,
        joint_prob: &mut F,
    ) -> LogProb {
        let mut subdensity = |likelihood_operands: &mut LikelihoodOperands| {
            let p = if vaf_tree_node.is_leaf() {
                joint_prob(likelihood_operands, data)
            } else if vaf_tree_node.is_branching() {
                LogProb::ln_sum_exp(
                    &vaf_tree_node
                        .children()
                        .iter()
                        .map(|child| {
                            self.density(
                                child,
                                &mut likelihood_operands.clone(),
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
                    likelihood_operands,
                    data,
                    biases,
                    joint_prob,
                )
            };

            assert!(!p.is_nan(), "bug: density has become NaN");
            p
        };

        match vaf_tree_node.kind() {
            grammar::vaftree::NodeKind::Log2FoldChange {
                sample_a,
                sample_b,
                predicate,
            } => {
                likelihood_operands.lfcs.push(VafLfc {
                    sample_a: *sample_a,
                    sample_b: *sample_b,
                    predicate: *predicate,
                });
                subdensity(likelihood_operands)
            }
            grammar::vaftree::NodeKind::False => LogProb::ln_zero(),
            grammar::vaftree::NodeKind::True => LogProb::ln_one(),
            grammar::vaftree::NodeKind::Sample { sample, vafs } => {
                let push_base_event =
                    |allele_freq, likelihood_operands: &mut LikelihoodOperands, is_discrete| {
                        likelihood_operands.events.insert(
                            *sample,
                            likelihood::Event {
                                allele_freq,
                                artifacts: biases.clone(),
                                is_discrete,
                            },
                        );
                    };

                let lfc_bounds = likelihood_operands.lfc_bounds(*sample);

                if let Some(bounds) = &lfc_bounds {
                    if bounds.is_empty() {
                        // METHOD: The current set of log fold changes is impossible to satisfy.
                        // Hence, we can immediately return a probability of zero.
                        return LogProb::ln_zero();
                    }
                }

                let (n_obs, is_clear_ref) = {
                    let pileup = &data.pileups[*sample];
                    if pileup.read_observations().is_empty()
                        && !pileup.depth_observations().is_empty()
                    {
                        // CNV case
                        // METHOD: We cannot determine the ref support in general without alt allele frequencies.
                        // Hence we always assume that there is no clear ref support and do the full
                        // evaluation below.
                        // This should be fine since there should be much less CNV calls than small variants.
                        (pileup.depth_observations().len(), false)
                    } else {
                        // normal variants, only consider read observations for these heuristic markers
                        let n_obs = pileup.read_observations().len();
                        let is_clear_ref = n_obs > 10
                            && pileup
                                .read_observations()
                                .iter()
                                .all(|obs| obs.is_positive_ref_support());
                        (n_obs, is_clear_ref)
                    }
                };

                match vafs {
                    grammar::VAFSpectrum::Set(vafs) => {
                        if is_clear_ref && vafs.iter().all(|vaf| **vaf > 0.0) {
                            // METHOD: shortcut for the case that all obs support the reference but the vafs
                            // in this event are > 0. Then, we don't need to recurse further and can
                            // immediately stop, returning a probability of zero.
                            return LogProb::ln_zero();
                        }
                        let vafs = if let Some(lfc_bounds) = &lfc_bounds {
                            vafs.iter()
                                .filter(|vaf| lfc_bounds.contains(**vaf))
                                .cloned()
                                .collect()
                        } else {
                            vafs.clone()
                        };

                        if vafs.len() == 1 {
                            push_base_event(
                                *vafs.iter().next().unwrap(),
                                likelihood_operands,
                                true,
                            );

                            subdensity(likelihood_operands)
                        } else {
                            LogProb::ln_sum_exp(
                                &vafs
                                    .iter()
                                    .map(|vaf| {
                                        let mut likelihood_operands = likelihood_operands.clone();
                                        push_base_event(*vaf, &mut likelihood_operands, true);
                                        subdensity(&mut likelihood_operands)
                                    })
                                    .collect_vec(),
                            )
                        }
                    }
                    grammar::VAFSpectrum::Range(vafs) => {
                        let vafs = if let Some(bounds) = &lfc_bounds {
                            vafs.intersect(bounds)
                        } else {
                            vafs.clone()
                        };
                        if vafs.is_empty() {
                            // METHOD: empty interval, integral must be zero.
                            return LogProb::ln_zero();
                        }
                        dbg!(&vafs);
                        if is_clear_ref && (*vafs.start > 0.0) {
                            // METHOD: shortcut for the case that all obs support the reference but the vaf
                            // range in this event is > 0. Then, we don't need to recurse further and can
                            // immediately stop, returning a probability of zero.
                            return LogProb::ln_zero();
                        }

                        if vafs.is_singleton() {
                            // METHOD: interval represents a single value, no need
                            // to integrate.
                            let vaf = vafs.start;
                            push_base_event(vaf, &mut likelihood_operands, true);
                            return subdensity(&mut likelihood_operands);
                        }

                        let resolution = &self.resolutions[*sample];
                        let min_vaf = vafs.observable_min(n_obs);
                        let max_vaf = vafs.observable_max(n_obs);
                        assert!(min_vaf <= max_vaf, "bug: min_vaf > max_vaf");
                        let mut density = |vaf| {
                            let mut likelihood_operands = likelihood_operands.clone();
                            push_base_event(vaf, &mut likelihood_operands, false);
                            subdensity(&mut likelihood_operands)
                        };

                        if (max_vaf - min_vaf) < **resolution {
                            // METHOD: Interval too small for desired resolution.
                            // Just use 3 grid points.
                            LogProb::ln_simpsons_integrate_exp(
                                |_, vaf| density(AlleleFreq(vaf)),
                                *min_vaf,
                                *max_vaf,
                                3,
                            )
                        } else if n_obs < 5 {
                            // METHOD: Not enough observations to expect a unimodal density.
                            // Use 11 grid points.
                            LogProb::ln_simpsons_integrate_exp(
                                |_, vaf| density(AlleleFreq(vaf)),
                                *min_vaf,
                                *max_vaf,
                                11,
                            )
                        } else {
                            // METHOD: enough data and large enough interval, use adaptive integration
                            // at the desired resolution.
                            adaptive_integration::ln_integrate_exp(
                                density,
                                min_vaf,
                                max_vaf,
                                **resolution,
                            )
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
                        subdensity(likelihood_operands)
                    }
                } else if *positive {
                    // no SNV but branch requires the defined SNV, hence abort with prob 0
                    LogProb::ln_zero()
                } else {
                    // skip this node, as we don't have the defined SNV but it is negated
                    subdensity(likelihood_operands)
                }
            }
        }
    }
}

impl Posterior for GenericPosterior {
    type BaseEvent = LikelihoodOperands;
    type Event = model::Event;
    type Data = Data;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        event: &Self::Event,
        data: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
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
                    let mut likelihood_operands = LikelihoodOperands::default();
                    bias_prior
                        + self.density(node, &mut likelihood_operands, data, biases, joint_prob)
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
    type Event = LikelihoodOperands;
    type Data = Data;

    fn compute(&self, operands: &Self::Event, data: &Self::Data, cache: &mut Cache) -> LogProb {
        // Step 1: Check if sample VAFs are compliant with any defined log fold changes.
        // If not, quickly return probability zero.
        for lfc in &operands.lfcs {
            let vaf_a = operands.events[lfc.sample_a].allele_freq;
            let vaf_b = operands.events[lfc.sample_b].allele_freq;
            if !lfc.predicate.is_true(&Log2FoldChange::new(vaf_a, vaf_b)) {
                return LogProb::ln_zero();
            }
        }

        let mut p = LogProb::ln_one();

        // Step 2: Calculate joint likelihood of sample VAFs.
        for (((sample, event), pileup), inner) in operands
            .events
            .iter()
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
                                secondary: operands.events[by].clone(),
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
