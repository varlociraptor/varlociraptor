//! stats.rs
//!
//! Statistical utility functions for MSI analysis.
//!
//! This module provides precise statistical calculations including:
//! 1. PHRED score conversion
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
/// Linear probability in range [0.0, 1.0]
///
/// # Examples:
/// assert!((p0 - 1.0).abs() < 1e-6);
pub(crate) fn phred_to_prob(phred: f64) -> f64 {
    *Prob::from(PHREDProb(phred))
}

/// Convert a `usize` to `f64`, panicking if the value can't be represented exactly.
///
/// `f64` represents integers exactly only up to 2^53. This codebase's counts
/// (e.g region totals) never approach that scale in practice, but this
/// guards against silently producing a corrupted percentage/statistic if
/// that assumption is ever violated, instead of failing loudly.
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
/// * Percentage as f32 in range [0.0, 100.0]
/// * Returns 0.0 if denominator is zero (avoiding division by zero)
///
/// # Examples
/// assert_eq!(calculate_percentage(5, 100), 5.0);`
pub fn calculate_percentage(numerator: usize, denominator: usize) -> f32 {
    if denominator == 0 {
        return 0.0;
    }

    ((usize_to_f64_exact(numerator) / usize_to_f64_exact(denominator)) * 100.0) as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::constants::test_constants::TEST_EPSILON;

    #[test]
    fn test_phred_to_prob() {
        // Key boundary: PHRED 0 = probability 1.0
        assert!((phred_to_prob(0.0) - 1.0).abs() < TEST_EPSILON);
        // Common value: PHRED 10 = probability 0.1
        assert!((phred_to_prob(10.0) - 0.1).abs() < TEST_EPSILON);
    }

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
    fn test_usize_to_f64_exact_panics_above_boundary() {
        let over = (1usize << 53) + 1;
        usize_to_f64_exact(over);
    }

    #[test]
    fn test_calculate_percentage() {
        assert_eq!(calculate_percentage(1, 4), 25.0);
        assert_eq!(calculate_percentage(0, 100), 0.0);
        assert_eq!(calculate_percentage(5, 0), 0.0);
    }
}
