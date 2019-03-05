// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::BTreeMap;

use itertools::Itertools;
use ordered_float::NotNan;
use rusty_machine::learning::lin_reg::LinRegressor;
use rusty_machine::learning::SupModel;
use rusty_machine::linalg::{Matrix, Vector};

use crate::model::AlleleFreq;

#[derive(Debug, Deserialize, Serialize)]
pub struct Estimate {
    pub observations: Vec<(f64, u64)>,
    pub intercept: f64,
    pub slope: f64,
}

impl Estimate {
    pub fn effective_mutation_rate(&self) -> f64 {
        self.slope
    }
}

pub fn estimate<F: IntoIterator<Item = AlleleFreq>>(allele_frequencies: F) -> Estimate {
    let mut observations = BTreeMap::new();
    for f in allele_frequencies {
        // count occurrences of 1 / f
        *observations
            .entry(NotNan::new(1.0).unwrap() / f)
            .or_insert(0) += 1u64;
    }
    let reciprocal_freqs = observations.keys().map(|f| **f).collect_vec();
    // calculate the cumulative sum
    let _counts = observations
        .values()
        .scan(0.0, |s, c| {
            *s += *c as f64;
            Some(*s)
        })
        .collect_vec();

    let observations = reciprocal_freqs
        .iter()
        .cloned()
        .zip(_counts.iter().map(|c| *c as u64))
        .collect_vec();

    let freqs = Matrix::new(reciprocal_freqs.len(), 1, reciprocal_freqs);
    let counts = Vector::new(_counts);

    let mut lin_mod = LinRegressor::default();
    lin_mod.train(&freqs, &counts);

    Estimate {
        observations: observations,
        intercept: lin_mod.parameters().unwrap()[0],
        slope: lin_mod.parameters().unwrap()[1],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::AlleleFreq;
    use itertools_num::linspace;

    #[test]
    fn test_estimate() {
        // example from Williams et al. Nature Genetics 2016.
        let freqs = linspace(0.12, 0.25, 2539).map(|af| AlleleFreq(af));
        let estimate = estimate(freqs);
        assert_relative_eq!(estimate.effective_mutation_rate(), 596.16, epsilon = 0.01);
    }
}
