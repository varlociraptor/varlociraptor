use statrs::function::factorial::ln_binomial

use crate::model::AlleleFreq;
use crate::utils::{PROB_025, PROB_05};


pub(crate) enum Inheritance {
    Mendelian { from: (usize, usize) },
    Clonal { from: usize },
}


pub(crate) struct Prior {
    ploidies: grammar::SampleInfo<Option<u32>>,
    germline_mutation_rate: grammar::SampleInfo<Option<f64>>,
    somatic_mutation_rate: grammar::SampleInfo<Option<f64>>,
    heterozygosity: Option<LogProb>,
    inheritance: grammar::SampleInfo<Option<Inheritance>>,
    cache: RefCell<BTreeMap<likelihood::Event, LogProb>>,
}


impl Prior {
    fn calc_prob(&self, event: &likelihood::Event, somatic_vafs: Vec<AlleleFreq>, germline_vafs: Vec<AlleleFreq>) -> LogProb {
        if somatic_vaf.len() == event.len() {
            // step 1: population
            let mut prob = if let Some(heterozygosity) = self.heterozygosity {
                // calculate population prior
                let population_samples = self.inheritance.iter().zip(self.ploidies.iter()).enumerate().filter_map(|(sample, (inheritance, ploidy))| if inheritance.is_none() && ploidy.is_some() { Some(sample) } else { None }).collect_vec();
                self.prob_population_germline(&self, &population_samples, &germline_vafs, heterozygosity);
            } else {
                LogProb::ln_one()
            };

            // step 2: inheritance
            prob += self.inheritance.iter().enumerate().map(|(sample, inheritance)| {
                match inheritance {
                    Some(Inheritance::Mendelian { from: parents }) => self.prob_mendelian_inheritance(sample, parents, &germline_vafs),
                    Some(Inheritance::Clonal { from: parent }) => self.prob_clonal_inheritance(sample, parent),
                }
            }).sum(); // product in log space

            // step 3: somatic mutations
            prob += self.somatic_mutation_rate.iter().enumerate().map(|(sample, somatic_mutation_rate)| {
                self.prob_somatic_mutation(sample, somatic_mutation_rate)
            }).sum(); // product in log space

            prob
        }
    }

    fn prob_somatic_mutation(&self, sample: usize, somatic_mutation_rate: LogProb) -> LogProb {
        
    }

    fn prob_clonal_inheritance(&self, sample: usize, parent: usize, germline_vafs: &[AlleleFreq]) -> LogProb {
        if relative_eq!(germline_vafs[sample], germline_vafs[parent]) {
            LogProb::ln_one()
        } else {
            LogProb::ln_zero()
        }
    }

    fn prob_population_germline(&self, population_samples: &[usize], germline_vafs: &[AlleleFreq], heterozygosity: LogProb) -> LogProb {
        let m = population_samples.iter().map(|sample| {
            (self.ploidies[sample].unwrap() * germline_vafs[sample]) as u32
        }).sum();

        let prob_m = |m| LogProb(heterozygosity - (m as f64).ln());

        if m > 0 {
            // m alt alleles
            prob_m
        } else {
            // no alt alleles
            let n = population_samples.iter().map(|sample| {
                self.ploidies[sample].unwrap()
            }).sum();
            LogProb::ln_sum_exp(&(1..n + 1).iter().map(prob_m).collect_vec())
        }
    }

    fn prob_mendelian_alt_counts(&self, ploidy: u32, source_alt: (u32, u32), target_alt: u32, germline_mutation_rate: LogProb) -> LogProb {
        LogProb::ln_sum_exp(&(0..source_alt.0).iter().product(0..source_alt.1).map(|(alt_from_first, alt_from_second)| {
            let choices_from_first = ln_binomial(source_alt.0, alt_from_first);
            let choices_from_second = ln_binomial(source_alt.1, alt_from_second);

            let prob = choices_from_first + LogProb::from(Prob(alt_from_first as f64 / ploidy as f64)) + 
                choices_from_second + LogProb::from(Prob(alt_from_second as f64 / ploidy as f64));
            
            let missing = target_alt as i32 - (alt_from_first + alt_from_second) as i32;
            prob + LogProb(germline_mutation_rate * missing)
        }).collect_vec());
    }

    fn prob_mendelian(&self, child: usize, parents: (usize, usize), germline_vafs: &[AlleleFreq]) -> LogProb {
        let ploidy = self.ploidies[child];
        for parent in parents {
            assert_eq!(self.ploidies[parent], ploidy);
        }
        let n_alt = |vaf| (vaf * ploidy as f64) as u32;

        self.prob_mendelian_alt_counts(ploidy, (n_alt(germline_vafs[parent.0]), n_alt(germline_vafs[parent.1])), n_alt(germline_vafs[child]), self.germline_mutation_rate[child].expect("bug: no germline VAF for child"))
    }
}


impl Prior for Prior {
    fn compute(&self, event: &Self::Event) -> LogProb {
        if let Some(prob) = self.cache.get().get(event) {
            prob
        } else {
            let mut prob = LogProb::ln_one();

            self.cache.get_mut().insert(event.to_owned(), prob);
        }
    }
}


impl model::modes::UpdatablePrior for Prior {
    fn set_universe(&mut self, universe: grammar::SampleInfo<grammar::VAFUniverse>) {
        // universe is ignored here
    }

    fn set_ploidies(&mut self, ploidies: grammar::SampleInfo<Option<u32>>) {
        self.ploidies = ploidies;
    }
}


pub(crate) struct CompositePrior<PP> 
where
    PP: Prior<Event=Vec<likelihood::Event>>
{
    sample_priors: grammar::SampleInfo<Vec<Box<SamplePrior>>>,
    population_prior: Option<PP>,
}

impl Prior for CompositePrior {
    type Event = Vec<likelihood::Event>;

    fn compute(&self, event: &Self::Event) -> LogProb {
        let mut prob = if let Some(population_prior) = self.population_prior {
            population_prior.compute(event)
        } else {
            LogProb::ln_one()
        };

        for (sample, sample_priors) in self.sample_priors.iter().enumerate() {
            if prob == LogProb::ln_zero() {
                return prob;
            }
            let mut sample_vaf = event[sample].allele_freq;
            for sample_prior in sample_priors {
                let res = sample_prior.compute(sample_vaf, event);
                sample_vaf -= res.explained_vaf;
                prob += res.prob;
            }
        }
    }

}


pub(crate) trait SamplePrior {
    fn compute(&self, sample_vaf: AlleleFreq, context_event: &Vec<likelihood::Event>) -> LogProb;
}