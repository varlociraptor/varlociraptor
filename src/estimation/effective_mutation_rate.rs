use std::collections::BTreeMap;

use rusty_machine::learning::lin_reg::LinRegressor;
use rusty_machine::learning::SupModel;
use rusty_machine::linalg::{Matrix, Vector};
use ordered_float::NotNaN;
use itertools::Itertools;


#[derive(Debug)]
pub struct Estimate {
    pub effective_mutation_rate: f64,
    pub observations: BTreeMap<NotNaN<f64>, u64>,
    pub model: LinRegressor
}


/// Estimate the effective mutation rate from given allele frequencies of presumable somatic passenger
/// mutations.
///
/// # Arguments
///
/// * `allele_frequencies` - allele frequencies of somatic passenger mutations (e.g. SNPs)
pub fn estimate<F: IntoIterator<Item=f64>>(allele_frequencies: F) -> Estimate {
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

    let slope = lin_mod.parameters().unwrap()[1];

    Estimate {
        effective_mutation_rate: slope,
        observations: observations,
        model: lin_mod
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
