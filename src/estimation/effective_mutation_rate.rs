use std::collections::BTreeMap;

use rusty_machine::learning::lin_reg::LinRegressor;
use rusty_machine::learning::SupModel;
use rusty_machine::linalg::{Matrix, Vector};
use ordered_float::NotNaN;

/// Estimate the effective mutation rate from given allele frequencies of presumable somatic passenger
/// mutations.
///
/// # Arguments
///
/// * `allele_frequencies` - allele frequencies of somatic passenger mutations (e.g. SNPs)
pub fn estimate<F: IntoIterator<Item=NotNaN<f64>>>(allele_frequencies: F) {
    let mut counter = BTreeMap::new();
    for f in allele_frequencies {
        // count occurrences of 1 / f
        *counter.entry(NotNaN::new(1.0 / *f).unwrap()).or_insert(0) += 1;
    }
    //let freqs = counter.keys().collect_vec();
    //let freqs = Matrix::new(4, 1, allele_frequencies);
}
