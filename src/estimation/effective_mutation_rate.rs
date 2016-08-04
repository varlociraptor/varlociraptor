use std::collections::BTreeMap;
use std::collections::btree_map;

use rusty_machine::learning::lin_reg::LinRegressor;
use rusty_machine::learning::SupModel;
use rusty_machine::linalg::{Matrix, Vector};
use ordered_float::NotNaN;
use itertools::Itertools;


#[derive(Debug)]
pub struct Estimator {
    observations: BTreeMap<NotNaN<f64>, u64>,
    model: LinRegressor
}


impl Estimator {
    pub fn train<F: IntoIterator<Item=f64>>(allele_frequencies: F) -> Self {
        let mut observations = BTreeMap::new();
        for f in allele_frequencies {
            // count occurrences of 1 / f
            *observations.entry(NotNaN::new(1.0 / f).unwrap()).or_insert(0) += 1;
        }
        let freqs = observations.keys().map(|f| **f).collect_vec();
        // calculate the cumulative sum
        let counts = observations.values().scan(0.0, |s, c| {
            *s += *c as f64;
            Some(*s)
        }).collect_vec();

        let freqs = Matrix::new(freqs.len(), 1, freqs);
        let counts = Vector::new(counts);

        let mut lin_mod = LinRegressor::default();
        lin_mod.train(&freqs, &counts);

        Estimator {
            observations: observations,
            model: lin_mod
        }
    }

    pub fn effective_mutation_rate(&self) -> f64 {
        self.slope()
    }

    pub fn slope(&self) -> f64 {
        self.model.parameters().unwrap()[1]
    }

    pub fn intercept(&self) -> f64 {
        self.model.parameters().unwrap()[0]
    }

    pub fn observations(&self) -> btree_map::Iter<NotNaN<f64>, u64> {
        self.observations.iter()
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use itertools::linspace;

    #[test]
    fn test_estimate() {
        // example from Williams et al. Nature Genetics 2016.
        let freqs = linspace(0.12, 0.25, 2539);
        let estimate = estimate(freqs);
        assert_relative_eq!(estimate.effective_mutation_rate, 596.16, epsilon=0.01);
    }
}
