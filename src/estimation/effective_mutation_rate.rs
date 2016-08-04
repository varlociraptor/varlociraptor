use std::collections::BTreeMap;

use rusty_machine::learning::lin_reg::LinRegressor;
use rusty_machine::learning::SupModel;
use rusty_machine::linalg::{Matrix, Vector};
use ordered_float::NotNaN;
use itertools::Itertools;


pub struct Estimate {
    pub effective_mutation_rate: f64,
    pub observations: BTreeMap<NotNaN<f64>, u64>
}


/// Estimate the effective mutation rate from given allele frequencies of presumable somatic passenger
/// mutations.
///
/// # Arguments
///
/// * `allele_frequencies` - allele frequencies of somatic passenger mutations (e.g. SNPs)
pub fn estimate<F: IntoIterator<Item=NotNaN<f64>>>(allele_frequencies: F) -> Estimate {
    let mut observations = BTreeMap::new();
    for f in allele_frequencies {
        // count occurrences of 1 / f
        *observations.entry(NotNaN::new(1.0 / *f).unwrap()).or_insert(0) += 1;
    }
    let freqs = observations.keys().map(|f| **f).collect_vec();
    let counts = observations.values().map(|c| *c as f64).collect_vec();

    let freqs = Matrix::new(freqs.len(), 1, freqs);
    let counts = Vector::new(counts);

    let mut lin_mod = LinRegressor::default();
    lin_mod.train(&freqs, &counts);

    let slope = lin_mod.parameters().unwrap()[0];

    Estimate {
        effective_mutation_rate: slope,
        observations: observations
    }
}
