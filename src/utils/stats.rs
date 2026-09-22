//! stats.rs
//!
//! Statistical utility functions for `Varlociraptor`.
//!
//! This module provides precise statistical calculations including:
//! 1. PHRED score conversion to linear probability
//! 2. Checked usize-to-f64 conversion (panics on precision loss)
//! 3. Percentage computation

use bio::stats::{PHREDProb, Prob};

/// Convert PHRED-scaled probability to linear probability.
///
/// PHRED scores encode error probabilities on a logarithmic scale where:
/// - PHRED = -10 × log₁₀(P)
/// - P = 10^(-PHRED/10)
///
/// Common PHRED values:
/// - PHRED 0 = P(1.0) = 100% probability
/// - PHRED 10 = P(0.1) = 10% probability
///
/// # Arguments
/// * `phred` - PHRED-scaled probability score
///
/// # Returns
/// Mathematically in `(0.0, 1.0]`, linear probability for valid (non-negative)
/// PHRED input - it only approaches but never reaches 0 as PHRED grows.
///
/// # Panics
/// Panics if `phred` is negative - a negative PHRED score is not meaningful
/// (PHRED = -10·log₁₀(P) is always >= 0 for P in a valid [0, 1] range).
/// Also panics on NaN and on infinite input: a real quality/probability
/// computation never legitimately produces literal infinity, so its presence
/// signals upstream corruption (e.g. a division by zero), not a valid
/// score.
pub(crate) fn phred_to_prob(phred: f64) -> f64 {
    assert!(
        phred.is_finite() && phred >= 0.0,
        "PHRED score must be finite and non-negative, got {}",
        phred
    );
    *Prob::from(PHREDProb(phred))
}

/// Convert a `usize` to `f64`, panicking if the value can't be represented exactly.
///
/// `f64` represents integers exactly only up to 2^53. This codebase's counts
/// (e.g region totals) never approach that scale in practice, but this
/// guards against silently producing a corrupted percentage/statistic if
/// that assumption is ever violated, instead of failing loudly.
///
/// # Arguments
/// * `value` - Count to convert (e.g. a region or allele total)
///
/// # Returns
/// `value` as `f64`, exact for any input within range.
///
/// # Panics
/// Panics if `value > 2^53` (9_007_199_254_740_992), where the conversion
/// would lose precision.
pub(crate) fn usize_to_f64_exact(value: usize) -> f64 {
    const MAX_EXACT_F64_INT: usize = 1 << 53;

    assert!(
        value <= MAX_EXACT_F64_INT,
        "usize value {} exceeds f64 exact-integer range (2^53); conversion would lose precision",
        value
    );

    value as f64
}

/// Calculate percentage as a simple f32 ratio.
///
/// # Arguments
/// * `numerator` - Count value (e.g., number of unstable regions)
/// * `denominator` - Total value (e.g., total number of regions)
///
/// # Returns
/// * `numerator / denominator × 100`, as f32
///
/// # Panics
/// Panics if `denominator == 0` - division by zero is an error, not
/// a value this function can meaningfully return.
///
/// # Note
/// usize_to_f64_exact (not a plain `as f64` cast) is used for the raw values so a
/// count large enough to lose precision panics loudly here instead of silently
/// erroring the percentage.
pub(crate) fn calculate_percentage(numerator: usize, denominator: usize) -> f32 {
    assert!(denominator != 0, "denominator must not be zero");

    ((usize_to_f64_exact(numerator) / usize_to_f64_exact(denominator)) * 100.0) as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::constants::test_constants::TEST_EPSILON;

    /* ========== phred_to_prob tests ================ */

    #[test]
    fn test_phred_to_prob() {
        // Key boundary: PHRED 0 = probability 1.0
        assert!((phred_to_prob(0.0) - 1.0).abs() < TEST_EPSILON);
        // Common value: PHRED 10 = probability 0.1
        assert!((phred_to_prob(10.0) - 0.1).abs() < TEST_EPSILON);
    }

    #[test]
    #[should_panic(expected = "PHRED score must be finite and non-negative")]
    fn test_phred_to_prob_negative_panics() {
        phred_to_prob(-10.0);
    }

    #[test]
    #[should_panic(expected = "PHRED score must be finite and non-negative")]
    fn test_phred_to_prob_nan_panics() {
        phred_to_prob(f64::NAN);
    }

    #[test]
    #[should_panic(expected = "PHRED score must be finite and non-negative")]
    fn test_phred_to_prob_infinity_panics() {
        phred_to_prob(f64::INFINITY);
    }

    #[test]
    #[should_panic(expected = "PHRED score must be finite and non-negative")]
    fn test_phred_to_prob_negative_infinity_panics() {
        phred_to_prob(f64::NEG_INFINITY);
    }

    /* ========== usize_to_f64_exact tests =========== */

    #[test]
    fn test_usize_to_f64_exact_normal_value() {
        assert_eq!(usize_to_f64_exact(1_000_000), 1_000_000.0);
    }

    #[test]
    fn test_usize_to_f64_exact_at_boundary() {
        let boundary = 1usize << 53;
        assert_eq!(usize_to_f64_exact(boundary), boundary as f64);
    }

    #[test]
    #[should_panic(expected = "exceeds f64 exact-integer range")]
    fn test_usize_to_f64_exact_above_boundary_panics() {
        let over = (1usize << 53) + 1;
        usize_to_f64_exact(over);
    }

    /* ========== calculate_percentage tests ========= */

    #[test]
    fn test_calculate_percentage() {
        assert_eq!(calculate_percentage(1, 4), 25.0);
        assert_eq!(calculate_percentage(0, 100), 0.0);
        assert_eq!(calculate_percentage(200, 100), 200.0);
    }

    #[test]
    #[should_panic(expected = "denominator must not be zero")]
    fn test_calculate_percentage_zero_denominator_panics() {
        calculate_percentage(5, 0);
    }
}
