use std::cell::RefCell;
use std::cmp;
use std::collections::HashSet;
use std::str;

use anyhow::Result;
use bio::stats::bayesian;
use bio::stats::bayesian::model::Prior as PriorTrait;
use bio::stats::{LogProb, Prob};
use data_encoding::HEXUPPER;
use itertools::Itertools;
use itertools_num::linspace;
use lru::LruCache;
use ring::digest;
use serde_json::{self, json, Value};
use statrs::distribution::{self, Discrete};

use crate::errors;
use crate::grammar;
use crate::variants::model::modes::generic::LikelihoodOperands;
use crate::variants::model::{bias::Artifacts, likelihood, AlleleFreq, VariantType};

pub(crate) trait UpdatablePrior {
    fn set_universe_and_ploidies(
        &mut self,
        universe: grammar::SampleInfo<grammar::VAFUniverse>,
        ploidies: grammar::SampleInfo<Option<u32>>,
    );

    fn set_variant_type(&mut self, variant_type: VariantType);
}

pub(crate) trait CheckablePrior {
    fn check(&self) -> Result<()>;
}

#[derive(Debug, Clone)]
pub(crate) enum Inheritance {
    Mendelian { from: (usize, usize) },
    Clonal { from: usize, somatic: bool },
    Subclonal { from: usize },
}

#[derive(Derefable, Debug)]
pub struct Cache {
    #[deref(mutable)]
    inner: LruCache<Vec<AlleleFreq>, LogProb>,
}

impl Default for Cache {
    fn default() -> Self {
        Cache {
            inner: LruCache::new(1000),
        }
    }
}

#[derive(Debug, TypedBuilder, Default)]
pub(crate) struct Prior {
    uniform: grammar::SampleInfo<bool>,
    ploidies: Option<grammar::SampleInfo<Option<u32>>>,
    universe: Option<grammar::SampleInfo<grammar::VAFUniverse>>,
    germline_mutation_rate: grammar::SampleInfo<Option<f64>>,
    somatic_effective_mutation_rate: grammar::SampleInfo<Option<f64>>,
    heterozygosity: Option<LogProb>,
    inheritance: grammar::SampleInfo<Option<Inheritance>>,
    variant_type_fractions: grammar::VariantTypeFraction,
    #[builder(default)]
    variant_type: Option<VariantType>,
    #[builder(default)]
    cache: RefCell<Cache>,
}

impl Clone for Prior {
    fn clone(&self) -> Self {
        Prior {
            uniform: self.uniform.clone(),
            ploidies: self.ploidies.clone(),
            universe: self.universe.clone(),
            germline_mutation_rate: self.germline_mutation_rate.clone(),
            somatic_effective_mutation_rate: self.somatic_effective_mutation_rate.clone(),
            heterozygosity: self.heterozygosity,
            inheritance: self.inheritance.clone(),
            cache: RefCell::default(),
            variant_type_fractions: self.variant_type_fractions.clone(),
            variant_type: self.variant_type.clone(),
        }
    }
}

impl Prior {
    fn n_samples(&self) -> usize {
        self.germline_mutation_rate.len()
    }

    fn collect_events(&self, event: LikelihoodOperands, events: &mut Vec<LikelihoodOperands>) {
        if event.len() < self.universe.as_ref().unwrap().len() {
            let sample = event.len();
            let new_event = |vaf| likelihood::Event {
                allele_freq: vaf,
                artifacts: Artifacts::none(),
                is_discrete: false,
            };

            for vaf_spectrum in self.universe.as_ref().unwrap()[sample].iter() {
                match vaf_spectrum {
                    grammar::formula::VAFSpectrum::Set(vafs) => {
                        for vaf in vafs {
                            let mut next = event.clone();
                            next.push(new_event(*vaf));
                            self.collect_events(next, events);
                        }
                    }
                    grammar::formula::VAFSpectrum::Range(range) => {
                        for vaf in linspace(*range.start, *range.end, 5) {
                            let mut next = event.clone();
                            next.push(new_event(AlleleFreq(vaf)));
                            self.collect_events(next, events);
                        }
                    }
                }
            }
        } else {
            events.push(event);
        }
    }

    pub(crate) fn plot(
        &self,
        target_sample: &str,
        sample_names: &grammar::SampleInfo<String>,
    ) -> Result<()> {
        let mut blueprint =
            serde_json::from_str(include_str!("../../../templates/plots/prior.json"))?;
        let mut events = Vec::new();
        self.collect_events(LikelihoodOperands::new(), &mut events);

        let mut visited = HashSet::new();

        let data = events
            .iter()
            .filter_map(|event| {
                let prob = self.compute(event);

                let hash = HEXUPPER.encode(
                    digest::digest(
                        &digest::SHA256,
                        event
                            .iter()
                            .zip(sample_names.iter())
                            .filter_map(|(e, sample)| {
                                if sample != target_sample {
                                    Some(
                                        json!({
                                            "sample": sample,
                                            "vaf": *e.allele_freq,
                                        })
                                        .to_string(),
                                    )
                                } else {
                                    None
                                }
                            })
                            .join(",")
                            .as_bytes(),
                    )
                    .as_ref(),
                )[..8]
                    .to_owned();

                let mut any_prob = false;
                let records = event
                    .iter()
                    .zip(sample_names.iter())
                    .filter_map(|(e, sample)| {
                        if sample == target_sample {
                            let prob = *Prob::from(prob);
                            if prob == 0.0 {
                                return None;
                            }
                            any_prob = true;
                            Some(json!({
                                "sample": sample.to_owned(),
                                "prob": prob,
                                "vaf": *e.allele_freq,
                                "hash": hash.clone(),
                            }))
                        } else if !visited.contains(&hash) {
                            Some(json!({
                                "sample": sample.to_owned(),
                                "vaf": *e.allele_freq,
                                "hash": hash.clone(),
                            }))
                        } else {
                            None
                        }
                    })
                    .collect_vec();
                if any_prob {
                    visited.insert(hash);
                    Some(records)
                } else {
                    None
                }
            })
            .flatten()
            .collect_vec();

        if let Value::Object(ref mut blueprint) = blueprint {
            blueprint["data"]["values"] = json!(data);
            blueprint["spec"]["layer"][1]["transform"][0]["filter"]["equal"] = json!(target_sample);
            blueprint["spec"]["layer"][0]["transform"][0]["filter"] =
                json!(format!("datum.sample != '{}'", target_sample));
            // print to STDOUT
            println!("{}", serde_json::to_string_pretty(blueprint)?);
            Ok(())
        } else {
            unreachable!();
        }
    }

    fn is_valid_germline_vaf(&self, sample: usize, vaf: AlleleFreq) -> bool {
        let ploidy = self.ploidies.as_ref().unwrap()[sample].expect("bug: ploidy not set");
        let n_alt = ploidy as f64 * *vaf;
        relative_eq!(n_alt, n_alt.round())
    }

    fn variant_type_fraction(&self) -> f64 {
        self.variant_type_fractions.get(
            self.variant_type
                .as_ref()
                .expect("bug: variant type not set in prior"),
        )
    }

    fn vartype_somatic_effective_mutation_rate(&self, sample: usize) -> Option<f64> {
        self.somatic_effective_mutation_rate[sample].map(|rate| rate * self.variant_type_fraction())
    }

    fn vartype_germline_mutation_rate(&self, sample: usize) -> Option<f64> {
        self.germline_mutation_rate[sample].map(|rate| rate * self.variant_type_fraction())
    }

    fn vartype_heterozygosity(&self) -> Option<LogProb> {
        self.heterozygosity
            .map(|het| LogProb((het.exp() * self.variant_type_fraction()).ln()))
    }

    fn has_somatic_variation(&self, sample: usize) -> bool {
        self.somatic_effective_mutation_rate[sample].is_some()
    }

    fn has_germline_variation(&self, sample: usize) -> bool {
        self.germline_mutation_rate[sample].is_some()
    }

    fn has_ploidy(&self, sample: usize) -> bool {
        self.ploidies.as_ref().unwrap()[sample].is_some()
    }

    fn has_uniform_prior(&self, sample: usize) -> bool {
        self.uniform[sample]
    }

    fn effective_somatic_vaf(
        &self,
        sample: usize,
        event: &LikelihoodOperands,
        germline_vafs: &[AlleleFreq],
    ) -> AlleleFreq {
        event[sample].allele_freq - germline_vafs[sample]
    }

    fn calc_prob(&self, event: &LikelihoodOperands, germline_vafs: Vec<AlleleFreq>) -> LogProb {
        if germline_vafs.len() == event.len() {
            // recursion end

            // step 1: population
            let mut prob = if let Some(heterozygosity) = self.vartype_heterozygosity() {
                // calculate population prior
                let population_samples = self
                    .inheritance
                    .iter()
                    .zip(self.ploidies.as_ref().unwrap().iter())
                    .enumerate()
                    .filter_map(|(sample, (inheritance, ploidy))| {
                        if inheritance.is_none()
                            && ploidy.is_some()
                            && !self.has_uniform_prior(sample)
                        {
                            Some(sample)
                        } else {
                            None
                        }
                    })
                    .collect_vec();
                self.prob_population_germline(&population_samples, &germline_vafs, heterozygosity)
            } else {
                LogProb::ln_one()
            };

            // step 2: inheritance and somatic mutation
            prob += self
                .inheritance
                .iter()
                .enumerate()
                .filter_map(|(sample, inheritance)| {
                    if self.has_uniform_prior(sample) {
                        // if sample has a uniform prior, ignore any defined inheritance patterns
                        return None;
                    }
                    match inheritance {
                        Some(Inheritance::Mendelian { from: parents }) => {
                            Some(self.prob_mendelian_inheritance(
                                sample,
                                *parents,
                                event,
                                &germline_vafs,
                            ))
                        }
                        Some(Inheritance::Clonal {
                            from: parent,
                            somatic,
                        }) => Some(self.prob_clonal_inheritance(
                            sample,
                            *parent,
                            event,
                            &germline_vafs,
                            *somatic,
                        )),
                        Some(Inheritance::Subclonal { from: parent }) => {
                            // here, somatic variation is already included in the returned probability
                            Some(self.prob_subclonal_inheritance(
                                sample,
                                *parent,
                                event,
                                &germline_vafs,
                            ))
                        }
                        None => {
                            // no inheritance pattern defined
                            self.vartype_somatic_effective_mutation_rate(sample)
                                .map(|r| {
                                    self.prob_somatic_mutation(
                                        r,
                                        self.effective_somatic_vaf(sample, event, &germline_vafs),
                                    )
                                })
                        }
                    }
                })
                .sum(); // product in log space

            assert!(*prob <= 0.0);

            prob
        } else {
            // recursion

            let sample = germline_vafs.len();
            let sample_event = &event[sample];
            let sample_ploidy = self.ploidies.as_ref().unwrap()[sample];
            let push_vafs = |germline| {
                let mut germline_vafs = germline_vafs.clone();
                germline_vafs.push(germline);

                germline_vafs
            };

            if sample_ploidy == Some(0) && *sample_event.allele_freq != 0.0 {
                // Chromosome does not occur in sample (e.g. Y chromsome) but allele freq > 0,
                // that's impossible, hence stop and return 0.
                LogProb::ln_zero()
            } else if self.has_uniform_prior(sample) {
                // sample has a uniform prior
                if self.universe.as_ref().unwrap()[sample].contains(event[sample].allele_freq) {
                    // no explicit info about germline VAF, assume 0.0
                    let germline_vafs = push_vafs(AlleleFreq(0.0));
                    self.calc_prob(event, germline_vafs)
                } else {
                    LogProb::ln_zero()
                }
            } else if self.has_somatic_variation(sample) {
                if let Some(ploidy) = self.ploidies.as_ref().unwrap()[sample] {
                    let mut probs = Vec::with_capacity(ploidy as usize + 1);
                    for n_alt in 0..=ploidy {
                        // for each possible number of germline alt alleles, obtain necessary somatic VAF to get the event VAF.
                        let germline_vaf = if ploidy > 0 {
                            AlleleFreq(n_alt as f64 / ploidy as f64)
                        } else {
                            AlleleFreq(0.0)
                        };

                        let germline_vafs = push_vafs(germline_vaf);
                        probs.push(self.calc_prob(event, germline_vafs));
                    }
                    LogProb::ln_sum_exp(&probs)
                } else {
                    unreachable!("bug: sample with somatic mutation rate but no ploidy")
                }
            } else if sample_ploidy.is_some() && self.heterozygosity.is_some() {
                if self.is_valid_germline_vaf(sample, sample_event.allele_freq) {
                    let germline_vafs = push_vafs(sample_event.allele_freq);

                    self.calc_prob(event, germline_vafs)
                } else {
                    // this vaf is impossible without somatic mutation
                    LogProb::ln_zero()
                }
            } else {
                unreachable!("bug: not enough info for prior but no universe specified");
            }
        }
    }

    fn prob_somatic_mutation(
        &self,
        somatic_effective_mutation_rate: f64,
        somatic_vaf: AlleleFreq,
    ) -> LogProb {
        // METHOD: we do not apply the model of Williams et al. The reason is that
        // too much can happen in addition to the somatic mutation (e.g. an overlapping SV).
        // Instead, we simply assume a flat prior over VAFs > 0.0, and distinguish between
        // 0.0 and >0.0 via the given mutation rate. In other words, the mutation rate
        // models the rate of loci with somatic mutations, but for those, any VAFs >0.0
        // are a priori equally possible.
        if relative_eq!(*somatic_vaf, 0.0) {
            LogProb(somatic_effective_mutation_rate.ln()).ln_one_minus_exp()
        } else {
            LogProb(somatic_effective_mutation_rate.ln())
        }
    }

    fn prob_clonal_inheritance(
        &self,
        sample: usize,
        parent: usize,
        event: &LikelihoodOperands,
        germline_vafs: &[AlleleFreq],
        somatic: bool,
    ) -> LogProb {
        if !relative_eq!(*germline_vafs[sample], *germline_vafs[parent]) {
            LogProb::ln_zero()
        } else {
            match (
                somatic,
                self.vartype_somatic_effective_mutation_rate(sample),
            ) {
                (true, Some(somatic_mutation_rate)) => {
                    // METHOD: de novo somatic variation in the sample, anything is possible.
                    let denovo_vaf = event[sample].allele_freq
                        - germline_vafs[sample]
                        - self.effective_somatic_vaf(parent, event, germline_vafs);
                    self.prob_somatic_mutation(somatic_mutation_rate, denovo_vaf)
                }
                (true, None) => {
                    // METHOD: somatic variation has to stay the same since it is inherited and unmodified.
                    if relative_eq!(
                        *self.effective_somatic_vaf(sample, event, germline_vafs),
                        *self.effective_somatic_vaf(parent, event, germline_vafs)
                    ) {
                        LogProb::ln_one()
                    } else {
                        LogProb::ln_zero()
                    }
                }
                (false, Some(somatic_mutation_rate)) => {
                    // METHOD: no somatic inheritance, all effective somatic vaf must be de novo.
                    self.prob_somatic_mutation(
                        somatic_mutation_rate,
                        self.effective_somatic_vaf(sample, event, germline_vafs),
                    )
                }
                (false, None) => LogProb::ln_one(),
            }
        }
    }

    fn prob_subclonal_inheritance(
        &self,
        sample: usize,
        parent: usize,
        event: &LikelihoodOperands,
        germline_vafs: &[AlleleFreq],
    ) -> LogProb {
        let total_vaf = event[sample].allele_freq;
        let germline_vaf = germline_vafs[sample];
        if !relative_eq!(*germline_vaf, *germline_vafs[parent]) {
            LogProb::ln_zero()
        } else {
            let parent_somatic_vaf = self.effective_somatic_vaf(parent, event, germline_vafs);
            let parent_total_vaf = event[parent].allele_freq;

            match self.vartype_somatic_effective_mutation_rate(sample) {
                Some(somatic_mutation_rate) => {
                    if *parent_total_vaf == 0.0 && *germline_vaf == 0.0 {
                        // METHOD: variant is denovo in this sample, calculate somatic mutation probability
                        self.prob_somatic_mutation(somatic_mutation_rate, total_vaf)
                    } else {
                        // METHOD: inherited, consider uniform probability because anything can happen during subclonal inheritance
                        LogProb::ln_one()
                    }
                }
                None => {
                    // METHOD: somatic variation has to stay the same since it is inherited and unmodified (no own somatic mutation rate defined).
                    if relative_eq!(
                        *self.effective_somatic_vaf(sample, event, germline_vafs),
                        *self.effective_somatic_vaf(parent, event, germline_vafs)
                    ) {
                        LogProb::ln_one()
                    } else {
                        LogProb::ln_zero()
                    }
                }
            }
        }
    }

    fn prob_population_germline(
        &self,
        population_samples: &[usize],
        germline_vafs: &[AlleleFreq],
        heterozygosity: LogProb,
    ) -> LogProb {
        let m = population_samples
            .iter()
            .map(|sample| {
                // we control above that the vafs are valid for the ploidy, but the rounding ensures that there are no numeric glitches
                (self.ploidies.as_ref().unwrap()[*sample].unwrap() as f64 * *germline_vafs[*sample])
                    .round() as u32
            })
            .sum();

        let prob_m = |m| LogProb(*heterozygosity - (m as f64).ln());

        if m > 0 {
            // m alt alleles
            prob_m(m)
        } else {
            // no alt alleles
            let n: u32 = population_samples
                .iter()
                .map(|sample| self.ploidies.as_ref().unwrap()[*sample].unwrap())
                .sum();
            LogProb::ln_sum_exp(&(1..=n).into_iter().map(prob_m).collect_vec()).ln_one_minus_exp()
        }
    }

    fn prob_select_ref_alt_alleles(
        &self,
        ploidy: u32,
        source_alt: u32,
        target_alt: u32,
        target_ref: u32,
    ) -> LogProb {
        let urn = distribution::Hypergeometric::new(
            ploidy as u64,
            source_alt as u64,
            (target_alt + target_ref) as u64,
        )
        .unwrap();
        LogProb::from(Prob(urn.pmf(target_alt as u64)))
    }

    fn prob_mendelian_alt_counts(
        &self,
        source_ploidy: (u32, u32),
        target_ploidy: u32,
        source_alt: (u32, u32),
        target_alt: u32,
        germline_mutation_rate: f64,
    ) -> LogProb {
        let prob_after_meiotic_split = |first_split_ploidy, second_split_ploidy| {
            (0..=cmp::min(source_alt.0, first_split_ploidy))
                .into_iter()
                .cartesian_product(0..=cmp::min(source_alt.1, second_split_ploidy))
                .filter_map(|(alt_from_first, alt_from_second)| {
                    // There may not be more alts from first and second than in the target
                    // but there may be more alts in the target due to denovo mutations.
                    if alt_from_first + alt_from_second <= target_alt {
                        let ref_from_first = first_split_ploidy - alt_from_first;
                        let ref_from_second = second_split_ploidy - alt_from_second;
                        let prob = self.prob_select_ref_alt_alleles(
                            source_ploidy.0,
                            source_alt.0,
                            alt_from_first,
                            ref_from_first,
                        ) + self.prob_select_ref_alt_alleles(
                            source_ploidy.1,
                            source_alt.1,
                            alt_from_second,
                            ref_from_second,
                        );

                        let missing = target_alt as i32 - (alt_from_first + alt_from_second) as i32;
                        Some(prob + LogProb(germline_mutation_rate.ln() * missing as f64))
                    } else {
                        None
                    }
                })
                .collect_vec()
        };

        let parent_inheritance_cases = |parental_ploidy| {
            if parental_ploidy % 2 == 0 {
                vec![parental_ploidy / 2]
            } else {
                let half = parental_ploidy as f64 / 2.0;
                vec![half.floor() as u32, half.ceil() as u32]
            }
        };
        let inheritance_cases = |p1, p2| {
            parent_inheritance_cases(p1)
                .into_iter()
                .cartesian_product(parent_inheritance_cases(p2).into_iter())
        };
        let is_valid_inheritance =
            |p1, p2, c| inheritance_cases(p1, p2).any(|(p1, p2)| c == p1 + p2);

        if is_valid_inheritance(source_ploidy.0, source_ploidy.1, target_ploidy) {
            LogProb::ln_sum_exp(
                &inheritance_cases(source_ploidy.0, source_ploidy.1)
                    .filter_map(|(p1, p2)| {
                        if p1 + p2 == target_ploidy {
                            Some(prob_after_meiotic_split(p1, p2))
                        } else {
                            None
                        }
                    })
                    .flatten()
                    .collect::<Vec<_>>(),
            )
        } else {
            // something went wrong, chromosomes of child do not match the parents
            // TODO For the future:
            // case 1: no separation in the first meiotic split (choose from all chromosomes of that parent)
            // case 2: no separation in the second meiotic split (duplicate a parental chromosome)
            panic!(
                "ploidies of child and parents do not match ({}, {} => {}) chromosome duplication events \
                    (e.g. trisomy) are not yet supported by the mendelian inheritance model of varlociraptor", 
                source_ploidy.0, source_ploidy.1, target_ploidy
            );
        }
    }

    fn prob_mendelian_inheritance(
        &self,
        child: usize,
        parents: (usize, usize),
        event: &LikelihoodOperands,
        germline_vafs: &[AlleleFreq],
    ) -> LogProb {
        let ploidies = self.ploidies.as_ref().unwrap();

        let ploidy = |sample: usize| ploidies[sample].unwrap();
        // we control above that the vafs are valid for the ploidy, but the rounding ensures that there are no numeric glitches
        let n_alt = |sample: usize| (*germline_vafs[sample] * ploidy(sample) as f64).round() as u32;

        let mut prob = self.prob_mendelian_alt_counts(
            (ploidy(parents.0), ploidy(parents.1)),
            ploidy(child),
            (n_alt(parents.0), n_alt(parents.1)),
            n_alt(child),
            self.vartype_germline_mutation_rate(child)
                .expect("bug: no germline mutation rate for child"),
        );

        if let Some(somatic_mutation_rate) = self.vartype_somatic_effective_mutation_rate(child) {
            prob += self.prob_somatic_mutation(
                somatic_mutation_rate,
                self.effective_somatic_vaf(child, event, germline_vafs),
            );
        }

        prob
    }
}

impl bayesian::model::Prior for Prior {
    type Event = LikelihoodOperands;

    fn compute(&self, event: &Self::Event) -> LogProb {
        if event.is_discrete() {
            let key: Vec<_> = event
                .iter()
                .map(|sample_event| sample_event.allele_freq)
                .collect();

            if let Some(prob) = self.cache.borrow_mut().get(&key) {
                return *prob;
            }
            let prob = self.calc_prob(event, Vec::with_capacity(event.len()));
            self.cache.borrow_mut().put(key, prob);

            prob
        } else {
            // METHOD: No caching for events with continuous VAFs as they are unlikely to reoccur.
            self.calc_prob(event, Vec::with_capacity(event.len()))
        }
    }
}

impl UpdatablePrior for Prior {
    fn set_universe_and_ploidies(
        &mut self,
        universe: grammar::SampleInfo<grammar::VAFUniverse>,
        ploidies: grammar::SampleInfo<Option<u32>>,
    ) {
        self.cache.borrow_mut().clear();
        self.universe = Some(universe);
        self.ploidies = Some(ploidies);
    }

    fn set_variant_type(&mut self, variant_type: VariantType) {
        self.variant_type = Some(variant_type);
    }
}

impl CheckablePrior for Prior {
    fn check(&self) -> Result<()> {
        let err = |msg: &str| {
            Err(errors::Error::InvalidPriorConfiguration {
                msg: msg.to_owned(),
            }
            .into())
        };
        for sample in 0..self.n_samples() {
            if let Some(inheritance) = &self.inheritance[sample] {
                if match inheritance {
                    Inheritance::Mendelian { from: (p1, p2) }
                        if !self.has_ploidy(*p1) || !self.has_ploidy(*p2) =>
                    {
                        true
                    }
                    Inheritance::Clonal { from, .. } if !self.has_ploidy(*from) => true,
                    Inheritance::Subclonal { from, .. } if !self.has_ploidy(*from) => true,
                    _ => false,
                } {
                    return err("inheritance defined but parental samples do not have a ploidy: define ploidy for each sample or the species");
                }
                match inheritance {
                    Inheritance::Mendelian { .. } if !self.has_germline_variation(sample) => {
                        return err("mendelian inheritance but no germline mutation rate defined: define germline mutation rate for child samples or the species")
                    }
                    _ => ()
                }
            }
            if let Some(Inheritance::Subclonal { .. }) = &self.inheritance[sample] {
                if !self.has_somatic_variation(sample) {
                    return err("subclonal inheritance defined but no somatic mutation: define somatic effective mutation rate for sample that inherits");
                }
            }
        }
        Ok(())
    }
}
