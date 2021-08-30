// Copyright 2021 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::HashMap;
use std::convert::Into;
use std::hash::Hash;
use std::{
    fmt::Debug,
    ops::{Add, Div, Mul, Sub},
};

use bio::stats::probs::LogProb;
use itertools::Itertools;
use ordered_float::NotNan;

/// Integrate over an interval of type T with a given density function while trying to minimize
/// the number of grid points evaluated and still hit the maximum likelihood point.
/// This is achieved via a binary search over the grid points.
/// The assumption is that the density is unimodal. If that is not the case,
/// the binary search will not find the maximum and the integral can miss features.
pub(crate) fn ln_integrate_exp<T, F>(
    mut density: F,
    min_point: T,
    max_point: T,
    max_resolution: T,
) -> LogProb
where
    T: Copy
        + Add<Output = T>
        + Sub<Output = T>
        + Div<Output = T>
        + Div<NotNan<f64>, Output = T>
        + Mul<Output = T>
        + Into<f64>
        + From<f64>
        + Ord
        + Debug
        + Hash,
    F: FnMut(T) -> LogProb,
    f64: From<T>,
{
    let mut probs = HashMap::new();

    let mut grid_point = |point, probs: &mut HashMap<_, _>| {
        probs.insert(point, density(point));
        point
    };
    let middle_grid_point = |left: T, right: T| (right + left) / NotNan::new(2.0).unwrap();
    // METHOD:
    // Step 1: perform binary search for maximum likelihood point
    // Remember all points.
    let mut left = grid_point(min_point, &mut probs);
    let mut right = grid_point(max_point, &mut probs);

    while ((right - left) >= max_resolution) || left < right {
        let middle = grid_point(middle_grid_point(left, right), &mut probs);

        let left_prob = probs.get(&left).unwrap();
        let right_prob = probs.get(&right).unwrap();

        if left_prob > right_prob {
            // investigate left window more closely
            right = middle;
        } else {
            // investiage right window more closely
            left = middle;
        }
    }
    let sorted_grid_points: Vec<f64> = probs.keys().sorted().map(|point| (*point).into()).collect();

    // METHOD:
    // Step 2: add additional grid points around optimum?
    dbg!(&sorted_grid_points);

    // METHOD:
    // Step 3: integrate over grid points visited during the binary search.
    LogProb::ln_trapezoidal_integrate_grid_exp::<f64, _>(
        |_, g| *probs.get(&T::from(g)).unwrap(),
        &sorted_grid_points,
    )
}
