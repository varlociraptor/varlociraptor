use std::collections::BTreeMap;

use rusty_machine::learning::lin_reg::LinRegressor;
use rusty_machine::linalg::{Matrix, Vector};
use rusty_machine::learning::SupModel;
use ordered_float::NotNaN;
use itertools::Itertools;


#[derive(Debug, RustcDecodable, RustcEncodable)]
pub struct Estimate {
    pub observations: Vec<(f64, u64)>,
    pub intercept: f64,
    pub slope: f64
}


impl Estimate {
    pub fn effective_mutation_rate(&self) -> f64 {
        self.slope
    }
}


pub fn estimate<F: IntoIterator<Item=f64>>(allele_frequencies: F) -> Estimate {
    let mut observations = BTreeMap::new();
    for f in allele_frequencies {
        // count occurrences of 1 / f
        *observations.entry(NotNaN::new(1.0 / f).unwrap()).or_insert(0) += 1u64;
    }
    let reciprocal_freqs = observations.keys().map(|f| **f).collect_vec();
    // calculate the cumulative sum
    let _counts = observations.values().scan(0.0, |s, c| {
        *s += *c as f64;
        Some(*s)
    }).collect_vec();

    let observations = reciprocal_freqs.iter().cloned().zip(_counts.iter().map(|c| *c as u64)).collect_vec();

    let freqs = Matrix::new(reciprocal_freqs.len(), 1, reciprocal_freqs);
    let counts = Vector::new(_counts);

    let mut lin_mod = LinRegressor::default();
    lin_mod.train(&freqs, &counts);

    Estimate {
        observations: observations,
        intercept: lin_mod.parameters().unwrap()[0],
        slope: lin_mod.parameters().unwrap()[1]
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
        assert_relative_eq!(estimate.effective_mutation_rate(), 596.16, epsilon=0.01);
    }
}
