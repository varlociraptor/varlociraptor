// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::ops::Range;

use bio::stats::LogProb;
use itertools::Itertools;
use rgsl::randist::gaussian::ugaussian_P;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::utils::NUMERICAL_EPSILON;
use crate::variants::sampling_bias::SamplingBias;
use crate::variants::types::Variant;

pub(crate) trait FragmentSamplingBias: Variant + SamplingBias {
    /// Get range of insert sizes with probability above zero.
    /// We use 6 SDs around the mean.
    fn isize_pmf_range(&self, alignment_properties: &AlignmentProperties) -> Range<u64> {
        let m = alignment_properties.insert_size().mean.round() as u64;
        let s = alignment_properties.insert_size().sd.ceil() as u64 * 6;
        m.saturating_sub(s)..m + s
    }

    /// Get probability of given insert size from distribution shifted by the given value.
    fn isize_pmf(
        &self,
        insert_size: u64,
        shift: f64,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        isize_pmf(
            insert_size as f64,
            alignment_properties.insert_size().mean + shift,
            alignment_properties.insert_size().sd,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn expected_prob_sample_alt(
        &self,
        left_read_len: u64,
        right_read_len: u64,
        left_feasible: u64,
        right_feasible: u64,
        delta_ref: u64,
        delta_alt: u64,
        enclose_only: bool,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        let mut infeasible_read_pos_left = left_read_len.saturating_sub(left_feasible);
        let mut infeasible_read_pos_right = right_read_len.saturating_sub(right_feasible);
        if !enclose_only {
            // If we don't need enclosing, the fragment may also overlap on the far right or left.
            // Hence, the are less infeasible positions (by the maximal feasible overlap).
            infeasible_read_pos_left = infeasible_read_pos_left.saturating_sub(left_feasible);
            infeasible_read_pos_right = infeasible_read_pos_right.saturating_sub(right_feasible);
        }
        let infeasible_read_pos = infeasible_read_pos_left + infeasible_read_pos_right;

        // Calculate over each possible true insert size instead of the concrete insert size
        // of this fragment. This way we get a global value. Otherwise, we would obtain a
        // bias for each fragment that is too small to enclose the variant,
        // e.g., a ref fragment.
        let expected_p_alt = LogProb::ln_sum_exp(
            &self
                .isize_pmf_range(alignment_properties)
                .filter_map(|x| {
                    let internal_segment = x
                        .saturating_sub(left_read_len)
                        .saturating_sub(right_read_len);
                    // Number of posititions in the internal segment where variant may not start,
                    // because it would lead to zero overlap.
                    let infeasible_internal_pos_alt =
                        (internal_segment + 1).saturating_sub(delta_alt);
                    let infeasible_pos_alt = infeasible_read_pos + infeasible_internal_pos_alt;
                    let infeasible_pos_ref = (internal_segment + 1).saturating_sub(delta_ref);

                    let valid_pos_alt = x
                        .saturating_sub(delta_alt)
                        .saturating_sub(infeasible_pos_alt);
                    let valid_pos_ref = x.saturating_sub(infeasible_pos_ref);

                    if x <= delta_alt
                        || x <= delta_alt + infeasible_pos_alt
                        || x <= infeasible_pos_ref
                    {
                        // if x is too small to enclose the variant, we skip it as it adds zero to the sum
                        None
                    } else {
                        // probability to sample a valid placement
                        let p = self.isize_pmf(x, 0.0, alignment_properties)
                            + LogProb((valid_pos_alt as f64).ln() - (valid_pos_ref as f64).ln());

                        assert!(
                            p.is_valid(),
                            "bug: invalid probability sampling probability {:?}",
                            p
                        );
                        Some(p)
                    }
                })
                .collect_vec(),
        )
        .cap_numerical_overshoot(NUMERICAL_EPSILON);

        assert!(expected_p_alt.is_valid());
        expected_p_alt
    }

    fn prob_sample_alt_fragment(
        &self,
        left_read_len: u64,
        right_read_len: u64,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        let left_feasible = self.feasible_bases(left_read_len, alignment_properties);
        let right_feasible = self.feasible_bases(right_read_len, alignment_properties);

        let delta_ref = self.len();
        let delta_alt = 0;
        self.expected_prob_sample_alt(
            left_read_len,
            right_read_len,
            left_feasible,
            right_feasible,
            delta_ref,
            delta_alt,
            true,
            alignment_properties,
        )
    }

    fn is_within_sd(
        &self,
        insert_size: u64,
        shift: f64,
        alignment_properties: &AlignmentProperties,
    ) -> bool {
        let m = alignment_properties.insert_size().mean + shift;
        (insert_size as f64 - m).abs() <= alignment_properties.insert_size().sd
    }
}

/// as shown in http://www.milefoot.com/math/stat/pdfc-normaldisc.htm
pub(crate) fn isize_pmf(value: f64, mean: f64, sd: f64) -> LogProb {
    LogProb((ugaussian_P((value + 0.5 - mean) / sd) - ugaussian_P((value - 0.5 - mean) / sd)).ln())
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     fn _test_n_fragment_positions(insert_size: u32) -> u32 {
//         // window to the left and right of the variant
//         let window = 364;
//         let read_len = 100;

//         n_fragment_positions(insert_size, read_len, read_len, window)
//     }

//     #[test]
//     fn test_n_fragment_positions_too_small() {
//         let n = _test_n_fragment_positions(150);
//         // Not enough space for placements around a centerpoint.
//         assert_eq!(n, 0);

//         let n = _test_n_fragment_positions(200);
//         // Not enough space for placements around a centerpoint.
//         assert_eq!(n, 0);
//     }

//     #[test]
//     fn test_n_fragment_positions_exact() {
//         let n = _test_n_fragment_positions(201);
//         // Enough space for 1 placement around a centerpoint.
//         assert_eq!(n, 1);
//     }

//     #[test]
//     fn test_n_fragment_positions() {
//         let n = _test_n_fragment_positions(202);
//         // Enough space for 2 placements around a centerpoint.
//         assert_eq!(n, 2);
//     }

//     #[test]
//     fn test_n_fragment_positions_too_large() {
//         let n = _test_n_fragment_positions(800);
//         assert_eq!(n, 0);
//     }

//     #[test]
//     fn test_isize_pmf() {
//         isize_pmf(300.0, 312.0, 15.0);
//     }
// }
