use std::cell::RefCell;
use std::collections::BTreeMap;

use bio::stats::bayesian;
use bio::stats::{LogProb, Prob};
use itertools::Itertools;
use statrs::function::factorial::ln_binomial;

use crate::grammar;
use crate::variants::model::{likelihood, AlleleFreq};

pub(crate) trait UpdatablePrior {
    fn set_universe(&mut self, universe: grammar::SampleInfo<grammar::VAFUniverse>);
    fn set_ploidies(&mut self, ploidies: grammar::SampleInfo<Option<u32>>);
}

pub(crate) enum Inheritance {
    Mendelian { from: (usize, usize) },
    Clonal { from: usize },
}

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
        LogProb(
            somatic_effective_mutation_rate.ln()
                - (2.0 * somatic_vafs[sample].ln() + (self.genome_size as f64).ln()),
        )
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

    fn prob_mendelian_alt_counts(
        &self,
        ploidy: u32,
        source_alt: (u32, u32),
        target_alt: u32,
        germline_mutation_rate: f64,
    ) -> LogProb {
        LogProb::ln_sum_exp(
            &(0..source_alt.0)
                .into_iter()
                .cartesian_product(0..source_alt.1)
                .map(|(alt_from_first, alt_from_second)| {
                    let choices_from_first =
                        ln_binomial(source_alt.0 as u64, alt_from_first as u64);
                    let choices_from_second =
                        ln_binomial(source_alt.1 as u64, alt_from_second as u64);

                    let prob = LogProb(
                        choices_from_first
                            + *LogProb::from(Prob(alt_from_first as f64 / ploidy as f64))
                            + choices_from_second
                            + *LogProb::from(Prob(alt_from_second as f64 / ploidy as f64)),
                    );

                    let missing = target_alt as i32 - (alt_from_first + alt_from_second) as i32;
                    prob + LogProb(germline_mutation_rate.ln() * missing as f64)
                })
                .collect_vec(),
        )
    }

    fn prob_mendelian_inheritance(
        &self,
        child: usize,
        parents: (usize, usize),
        germline_vafs: &[AlleleFreq],
    ) -> LogProb {
        let ploidy = self.ploidies.as_ref().unwrap()[child].unwrap();
        assert_eq!(self.ploidies.as_ref().unwrap()[parents.0].unwrap(), ploidy);
        assert_eq!(self.ploidies.as_ref().unwrap()[parents.1].unwrap(), ploidy);

        // we control above that the vafs are valid for the ploidy, but the rounding ensures that there are no numeric glitches
        let n_alt = |vaf: AlleleFreq| (*vaf * ploidy as f64).round() as u32;

        self.prob_mendelian_alt_counts(
            ploidy,
            (
                n_alt(germline_vafs[parents.0]),
                n_alt(germline_vafs[parents.1]),
            ),
            n_alt(germline_vafs[child]),
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
