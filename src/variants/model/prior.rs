use std::cell::RefCell;
use std::cmp;
use std::collections::BTreeMap;

use bio::stats::bayesian;
use bio::stats::bayesian::model::Prior as PriorTrait;
use bio::stats::{LogProb, Prob};
use itertools::Itertools;
use itertools_num::linspace;
use serde_json::json;
use statrs::function::factorial::ln_binomial;

use crate::grammar;
use crate::variants::model::{bias::Biases, likelihood, AlleleFreq};

pub(crate) trait UpdatablePrior {
    fn set_universe(&mut self, universe: grammar::SampleInfo<grammar::VAFUniverse>);
    fn set_ploidies(&mut self, ploidies: grammar::SampleInfo<Option<u32>>);
}

#[derive(Debug)]
pub(crate) enum Inheritance {
    Mendelian { from: (usize, usize) },
    Clonal { from: usize },
}

#[derive(Debug, TypedBuilder)]
pub(crate) struct Prior {
    ploidies: Option<grammar::SampleInfo<Option<u32>>>,
    universe: Option<grammar::SampleInfo<grammar::VAFUniverse>>,
    germline_mutation_rate: grammar::SampleInfo<Option<f64>>,
    somatic_effective_mutation_rate: grammar::SampleInfo<Option<f64>>,
    heterozygosity: Option<LogProb>,
    inheritance: grammar::SampleInfo<Option<Inheritance>>,
    genome_size: u64,
    cache: RefCell<BTreeMap<Vec<likelihood::Event>, LogProb>>,
}

impl Prior {
    fn collect_events(
        &self,
        event: Vec<likelihood::Event>,
        events: &mut Vec<Vec<likelihood::Event>>,
    ) {
        if event.len() < self.universe.unwrap().len() {
            let sample = event.len();
            let new_event = |vaf| likelihood::Event {
                allele_freq: vaf,
                biases: Biases::none(),
            };

            for vaf_spectrum in self.universe.unwrap()[sample].iter() {
                match vaf_spectrum {
                    grammar::formula::VAFSpectrum::Set(vafs) => {
                        for vaf in vafs {
                            let mut next = event.clone();
                            next.push(new_event(*vaf));
                            self.collect_events(next, events);
                        }
                    }
                    grammar::formula::VAFSpectrum::Range(range) => {
                        for vaf in linspace(*range.start, *range.end, 10) {
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

    pub(crate) fn plot(&self, sample_names: &grammar::SampleInfo<String>) -> String {
        use vega_lite_4::*;

        let mut events = Vec::new();
        self.collect_events(Vec::new(), &mut events);

        let data = events
            .iter()
            .map(|event| {
                let prob = self.compute(event);

                event.iter().zip(sample_names.iter()).map(|(e, sample)| {
                json!({
                    "sample": sample.to_owned(),
                    "prob": prob,
                    "vaf": *e.allele_freq,
                    "others": event.iter().zip(sample_names.iter()).filter_map(|(e, other_sample)| {
                        if sample != other_sample {
                            Some(format!("{}={}", other_sample, *e.allele_freq))
                        } else {
                            None
                        }
                    }).join(", ")
                })
            })
            })
            .flatten()
            .collect_vec();

        let chart = VegaliteBuilder::default()
            .data(&data)
            .mark({
                let mut mark = MarkDefClass::default();
                mark.def_type = Some(Mark::Line);
                mark.point = Some(true.into());
                mark
            })
            .facet({
                let mut facet = Facet::default();
                facet.field = Some("sample".into());
                facet
            })
            .encoding(
                EncodingBuilder::default()
                    .x(XClassBuilder::default()
                        .field("vaf")
                        .def_type(StandardType::Quantitative)
                        .build()
                        .unwrap())
                    .y(YClassBuilder::default()
                        .field("prob")
                        .def_type(StandardType::Quantitative)
                        .build()
                        .unwrap())
                    .color({
                        let mut def = DefWithConditionMarkPropFieldDefGradientStringNull::default();
                        def.field = Some("others".into());
                        def
                    })
                    .build()
                    .unwrap(),
            )
            .build()
            .unwrap();

        chart.to_string().unwrap()
    }

    fn is_valid_germline_vaf(&self, sample: usize, vaf: AlleleFreq) -> bool {
        if let Some(ploidy) = self.ploidies.as_ref().unwrap()[sample] {
            let n_alt = ploidy as f64 * *vaf;
            relative_eq!(n_alt, n_alt.round())
        } else {
            true
        }
    }

    fn has_somatic_variation(&self, sample: usize) -> bool {
        self.somatic_effective_mutation_rate[sample].is_some()
    }

    fn calc_prob(
        &self,
        event: &[likelihood::Event],
        somatic_vafs: Vec<AlleleFreq>,
        germline_vafs: Vec<AlleleFreq>,
    ) -> LogProb {
        if somatic_vafs.len() == event.len() {
            // recursion end

            // step 1: population
            let mut prob = if let Some(heterozygosity) = self.heterozygosity {
                // calculate population prior
                let population_samples = self
                    .inheritance
                    .iter()
                    .zip(self.ploidies.as_ref().unwrap().iter())
                    .enumerate()
                    .filter_map(|(sample, (inheritance, ploidy))| {
                        if inheritance.is_none() && ploidy.is_some() {
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

            // step 2: inheritance
            prob += self
                .inheritance
                .iter()
                .enumerate()
                .filter_map(|(sample, inheritance)| match inheritance {
                    Some(Inheritance::Mendelian { from: parents }) => {
                        Some(self.prob_mendelian_inheritance(sample, *parents, &germline_vafs))
                    }
                    Some(Inheritance::Clonal { from: parent }) => {
                        Some(self.prob_clonal_inheritance(sample, *parent, &germline_vafs))
                    }
                    None => None,
                })
                .sum(); // product in log space

            // step 3: somatic mutations
            prob += self
                .somatic_effective_mutation_rate
                .iter()
                .enumerate()
                .filter_map(|(sample, somatic_effective_mutation_rate)| {
                    somatic_effective_mutation_rate
                        .map(|r| self.prob_somatic_mutation(sample, r, &somatic_vafs))
                })
                .sum(); // product in log space

            prob
        } else {
            // recursion

            let sample = somatic_vafs.len();
            let sample_event = &event[sample];
            let sample_ploidy = self.ploidies.as_ref().unwrap()[sample];
            let push_vafs = |somatic, germline| {
                let mut somatic_vafs = somatic_vafs.clone();
                let mut germline_vafs = germline_vafs.clone();
                somatic_vafs.push(somatic);
                germline_vafs.push(germline);

                (somatic_vafs, germline_vafs)
            };

            if self.has_somatic_variation(sample) {
                if let Some(ploidy) = self.ploidies.as_ref().unwrap()[sample] {
                    let mut probs = Vec::with_capacity(ploidy as usize + 1);
                    for n_alt in 0..ploidy + 1 {
                        // for each possible number of germline alt alleles, obtain necessary somatic VAF to get the event VAF.
                        let germline_vaf = AlleleFreq(n_alt as f64 / ploidy as f64);
                        let somatic_vaf =
                            AlleleFreq((germline_vaf - sample_event.allele_freq).abs());
                        let (somatic_vafs, germline_vafs) = push_vafs(somatic_vaf, germline_vaf);
                        probs.push(self.calc_prob(event, somatic_vafs, germline_vafs));
                    }
                    LogProb::ln_sum_exp(&probs)
                } else {
                    unreachable!("bug: sample with somatic mutation rate but no ploidy")
                }
            } else if sample_ploidy.is_some() {
                if self.is_valid_germline_vaf(sample, sample_event.allele_freq) {
                    let (somatic_vafs, germline_vafs) =
                        push_vafs(AlleleFreq(0.0), sample_event.allele_freq);

                    self.calc_prob(event, somatic_vafs, germline_vafs)
                } else {
                    // this vaf is impossible without somatic mutation
                    LogProb::ln_zero()
                }
            } else {
                // sample has a uniform prior
                if self.universe.as_ref().unwrap()[sample].contains(sample_event.allele_freq) {
                    LogProb::ln_one()
                } else {
                    LogProb::ln_zero()
                }
            }
        }
    }

    fn prob_somatic_mutation(
        &self,
        sample: usize,
        somatic_effective_mutation_rate: f64,
        somatic_vafs: &[AlleleFreq],
    ) -> LogProb {
        let density = |vaf: f64| {
            LogProb(
                somatic_effective_mutation_rate.ln()
                    - (2.0 * vaf.ln() + (self.genome_size as f64).ln()),
            )
        };
        if *somatic_vafs[sample] == 0.0 {
            LogProb::ln_simpsons_integrate_exp(|_, vaf| density(vaf), 0.0, 1.0, 5)
                .ln_one_minus_exp()
        } else {
            density(*somatic_vafs[sample])
        }
    }

    fn prob_clonal_inheritance(
        &self,
        sample: usize,
        parent: usize,
        germline_vafs: &[AlleleFreq],
    ) -> LogProb {
        if relative_eq!(*germline_vafs[sample], *germline_vafs[parent]) {
            LogProb::ln_one()
        } else {
            LogProb::ln_zero()
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
            LogProb::ln_sum_exp(&(1..n + 1).into_iter().map(prob_m).collect_vec())
        }
    }

    fn prob_select_alleles(&self, ploidy: u32, source_n: u32, target_n: u32) -> LogProb {
        let choices = ln_binomial(source_n as u64, target_n as u64);
        LogProb(choices + *LogProb::from(Prob(target_n as f64 / source_n as f64)))
    }

    fn prob_select_ref_alt_alleles(
        &self,
        ploidy: u32,
        source_alt: u32,
        target_alt: u32,
        target_ref: u32,
    ) -> LogProb {
        self.prob_select_alleles(ploidy, source_alt, target_alt)
            + self.prob_select_alleles(ploidy, ploidy - source_alt, target_ref)
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
            (0..cmp::min(source_alt.0, first_split_ploidy))
                .into_iter()
                .cartesian_product(0..cmp::min(source_alt.1, second_split_ploidy))
                .map(|(alt_from_first, alt_from_second)| {
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
                    prob + LogProb(germline_mutation_rate.ln() * missing as f64)
                })
        };

        let probs = match (source_ploidy.0, source_ploidy.1, target_ploidy) {
            (p1, p2, c) if p1 % 2 == 0 && p2 % 2 == 0 && c == (p1 / 2 + p2 / 2) => {
                // Default case, normal meiosis (child inherits one half from each parent).
                prob_after_meiotic_split(p1 / 2, p2 / 2).collect_vec()
            }
            (0, p2, c) if p2 == c => {
                // e.g. monosomal inheritance from single parent (e.g. Y chromosome) or no meiotic split from that parent
                prob_after_meiotic_split(0, p2).collect_vec()
            }
            (p1, 0, c) if p1 == c => {
                // e.g. monosomal inheritance from single parent (e.g. Y chromosome) or no meiotic split from that parent
                prob_after_meiotic_split(p1, 0).collect_vec()
            }
            (p1, p2, c) => {
                // something went wrong, there are more chromosomes in the child than in the parents
                // case 1: no separation in the first meiotic split (choose from all chromosomes of that parent)
                // case 2: no separation in the second meiotic split (duplicate a parental chromosome)
                panic!(format!(
                    "ploidies of child and parents do not match ({}, {} => {}) chromosome duplication events \
                     (e.g. trisomy) are not yet supported by the mendelian inheritance model of varlociraptor", 
                    p1, p2, c
                ));
            }
        };

        LogProb::ln_sum_exp(&probs)
    }

    fn prob_mendelian_inheritance(
        &self,
        child: usize,
        parents: (usize, usize),
        germline_vafs: &[AlleleFreq],
    ) -> LogProb {
        let ploidies = self.ploidies.as_ref().unwrap();

        let ploidy = |sample: usize| ploidies[sample].unwrap();
        // we control above that the vafs are valid for the ploidy, but the rounding ensures that there are no numeric glitches
        let n_alt = |sample: usize| (*germline_vafs[sample] * ploidy(sample) as f64).round() as u32;

        self.prob_mendelian_alt_counts(
            (ploidy(parents.0), ploidy(parents.1)),
            ploidy(child),
            (n_alt(parents.0), n_alt(parents.1)),
            n_alt(child),
            self.germline_mutation_rate[child].expect("bug: no germline VAF for child"),
        )
    }
}

impl bayesian::model::Prior for Prior {
    type Event = Vec<likelihood::Event>;

    fn compute(&self, event: &Self::Event) -> LogProb {
        if let Some(prob) = self.cache.borrow().get(event) {
            *prob
        } else {
            let prob = self.calc_prob(
                event,
                Vec::with_capacity(event.len()),
                Vec::with_capacity(event.len()),
            );
            self.cache.borrow_mut().insert(event.to_owned(), prob);

            prob
        }
    }
}

impl UpdatablePrior for Prior {
    fn set_universe(&mut self, universe: grammar::SampleInfo<grammar::VAFUniverse>) {
        self.universe = Some(universe);
    }

    fn set_ploidies(&mut self, ploidies: grammar::SampleInfo<Option<u32>>) {
        self.ploidies = Some(ploidies);
    }
}
