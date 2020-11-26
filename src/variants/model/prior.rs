use statrs::function::factorial::ln_binomial;

use crate::model::AlleleFreq;
use crate::utils::{PROB_025, PROB_05};

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
    cache: RefCell<BTreeMap<likelihood::Event, LogProb>>,
}

impl Prior {
    fn is_valid_germline_vaf(&self, sample: usize, vaf: AlleleFreq) -> bool {
        if let Some(ploidy) = self.ploidies[sample] {
            let n_alt = (ploidy as f64 * vaf);
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
        prob_sum_cache: &mut Vec<LogProb>,
    ) -> LogProb {
        if somatic_vaf.len() == event.len() {
            // recursion end

            // step 1: population
            let mut prob = if let Some(heterozygosity) = self.heterozygosity {
                // calculate population prior
                let population_samples = self
                    .inheritance
                    .iter()
                    .zip(self.ploidies.iter())
                    .enumerate()
                    .filter_map(|(sample, (inheritance, ploidy))| {
                        if inheritance.is_none() && ploidy.is_some() {
                            Some(sample)
                        } else {
                            None
                        }
                    })
                    .collect_vec();
                self.prob_population_germline(
                    &self,
                    &population_samples,
                    &germline_vafs,
                    heterozygosity,
                );
            } else {
                LogProb::ln_one()
            };

            // step 2: inheritance
            prob += self
                .inheritance
                .iter()
                .enumerate()
                .map(|(sample, inheritance)| match inheritance {
                    Some(Inheritance::Mendelian { from: parents }) => {
                        self.prob_mendelian_inheritance(sample, parents, &germline_vafs)
                    }
                    Some(Inheritance::Clonal { from: parent }) => {
                        self.prob_clonal_inheritance(sample, parent)
                    }
                })
                .sum(); // product in log space

            // step 3: somatic mutations
            prob += self
                .somatic_mutation_rate
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
            let sample_event = event[sample];
            let sample_ploidy = self.ploidies[sample];
            let push_vafs = |somatic, germline| {
                let mut somatic_vafs = somatic_vafs.clone();
                let mut germline_vafs = germline_vafs.clone();
                somatic_vafs.push(somatic);
                germline_vafs.push(germline);

                (somatic_vafs, germline_vafs)
            };

            if self.has_somatic_variation(sample) {
                prob_sum_cache.clear();

                if let Some(ploidy) = self.ploidies[sample] {
                    for n_alt in 0..ploidy {
                        // for each possible number of germline alt alleles, obtain necessary somatic VAF to get the event VAF.
                        let germline_vaf = AlleleFreq(n_alt as f64 / ploidy as f64);
                        let somatic_vaf = (germline_vaf - sample_event.allele_freq).abs();
                        let (somatic_vafs, germline_vafs) = push_vafs(somatic_vaf, germline_vaf);
                        prob_sum_cache.push(self.calc_prob(event, somatic_vafs, germline_vafs));
                    }
                } else {
                    unreachable!("bug: sample with somatic mutation rate but no ploidy")
                }

                LogProb::ln_sum_exp(prob_sum_cache)
            } else if sample_ploidy.is_some() {
                if self.is_valid_germline_vaf(sample_event.allele_freq) {
                    let (somatic_vafs, germline_vafs) =
                        push_vafs(AlleleFreq(0.0), sample_event.allele_freq);

                    self.calc_prob(event, somatic_vafs, germline_vafs)
                } else {
                    // this vaf is impossible without somatic mutation
                    LogProb::ln_zero()
                }
            } else {
                // sample has a uniform prior
                if self.universe[sample].contains(sample_event.allele_freq) {
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
        somatic_effective_mutation_rate.ln()
            - (2.0 * somatic_vafs[sample].ln() + (self.genome_size as f64).ln())
    }

    fn prob_clonal_inheritance(
        &self,
        sample: usize,
        parent: usize,
        germline_vafs: &[AlleleFreq],
    ) -> LogProb {
        if relative_eq!(germline_vafs[sample], germline_vafs[parent]) {
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
                (self.ploidies[sample].unwrap() * germline_vafs[sample]).round() as u32
            })
            .sum();

        let prob_m = |m| LogProb(heterozygosity - (m as f64).ln());

        if m > 0 {
            // m alt alleles
            prob_m
        } else {
            // no alt alleles
            let n = population_samples
                .iter()
                .map(|sample| self.ploidies[sample].unwrap())
                .sum();
            LogProb::ln_sum_exp(&(1..n + 1).iter().map(prob_m).collect_vec())
        }
    }

    fn prob_mendelian_alt_counts(
        &self,
        ploidy: u32,
        source_alt: (u32, u32),
        target_alt: u32,
        germline_mutation_rate: LogProb,
    ) -> LogProb {
        LogProb::ln_sum_exp(
            &(0..source_alt.0)
                .iter()
                .product(0..source_alt.1)
                .map(|(alt_from_first, alt_from_second)| {
                    let choices_from_first = ln_binomial(source_alt.0, alt_from_first);
                    let choices_from_second = ln_binomial(source_alt.1, alt_from_second);

                    let prob = choices_from_first
                        + LogProb::from(Prob(alt_from_first as f64 / ploidy as f64))
                        + choices_from_second
                        + LogProb::from(Prob(alt_from_second as f64 / ploidy as f64));

                    let missing = target_alt as i32 - (alt_from_first + alt_from_second) as i32;
                    prob + LogProb(germline_mutation_rate * missing)
                })
                .collect_vec(),
        );
    }

    fn prob_mendelian(
        &self,
        child: usize,
        parents: (usize, usize),
        germline_vafs: &[AlleleFreq],
    ) -> LogProb {
        let ploidy = self.ploidies[child];
        for parent in parents {
            assert_eq!(self.ploidies[parent], ploidy);
        }
        // we control above that the vafs are valid for the ploidy, but the rounding ensures that there are no numeric glitches
        let n_alt = |vaf| (vaf * ploidy as f64).round() as u32;

        self.prob_mendelian_alt_counts(
            ploidy,
            (
                n_alt(germline_vafs[parent.0]),
                n_alt(germline_vafs[parent.1]),
            ),
            n_alt(germline_vafs[child]),
            self.germline_mutation_rate[child].expect("bug: no germline VAF for child"),
        )
    }
}

impl Prior for Prior {
    fn compute(&self, event: &Self::Event) -> LogProb {
        if let Some(prob) = self.cache.get().get(event) {
            prob
        } else {
            let prob = self.calc_prob(
                event,
                Vec::with_capacity(event.len()),
                Vec::with_capacity(event.len()),
                Vec::new(),
            );
            self.cache.get_mut().insert(event.to_owned(), prob);

            prob
        }
    }
}

impl model::modes::UpdatablePrior for Prior {
    fn set_universe(&mut self, universe: grammar::SampleInfo<grammar::VAFUniverse>) {
        self.universe = Some(universe);
    }

    fn set_ploidies(&mut self, ploidies: grammar::SampleInfo<Option<u32>>) {
        self.ploidies = Some(ploidies);
    }
}
