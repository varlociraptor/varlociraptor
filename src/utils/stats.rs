//! stats.rs
//!
//! Statistics utilities for precise calculations.

use bio::stats::{PHREDProb, Prob};
use rust_decimal::prelude::ToPrimitive;
use rust_decimal::Decimal;

/// Convert PHRED-scaled probability to linear probability.
/// Uses rust-bio's `PHREDProb` and `Prob` types for precise, idiomatic conversion.
pub(crate) fn phred_to_prob(phred: f64) -> f64 {
    *Prob::from(PHREDProb(phred))
}

/// Calculate percentage using exact decimal arithmetic
///
/// Prevents precision loss when converting large usize values to percentages.
///
/// # Arguments
/// * `numerator` - (e.g., Count: number of unstable regions)
/// * `denominator` - (e.g., Total: total regions)
///
/// # Returns
/// Percentage as f64 (0.0 to 100.0), or 0.0 if denominator is 0`
pub fn calculate_percentage_exact(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        return 0.0;
    }

    let num = Decimal::from(numerator);
    let den = Decimal::from(denominator);

    let ratio = num / den;

    let percentage = ratio * Decimal::from(100);

    percentage.to_f64().unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_phred_to_prob_basic() {
        // PHRED 0 = 1.0, PHRED 10 ≈ 0.1, PHRED 20 ≈ 0.01
        let p0 = phred_to_prob(0.0);
        let p10 = phred_to_prob(10.0);
        let p20 = phred_to_prob(20.0);

        assert!((p0 - 1.0).abs() < 1e-6);
        assert!((p10 - 0.1).abs() < 1e-3);
        assert!((p20 - 0.01).abs() < 1e-4);
    }

    #[test]
    fn test_basic_percentage() {
        assert_eq!(calculate_percentage_exact(5, 100), 5.0);
        assert_eq!(calculate_percentage_exact(1, 4), 25.0);
        assert_eq!(calculate_percentage_exact(0, 100), 0.0);
    }

    #[test]
    fn test_zero_denominator() {
        assert_eq!(calculate_percentage_exact(5, 0), 0.0);
        assert_eq!(calculate_percentage_exact(0, 0), 0.0);
    }

    #[test]
    fn test_large_values() {
        // Values exceeding f64 safe integer range (2^53)
        let pct = calculate_percentage_exact(
            10_000_000_000_000,  // 10 trillion
            100_000_000_000_000, // 100 trillion
        );
        assert!((pct - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_precision() {
        // Test that precision is maintained
        let pct = calculate_percentage_exact(1, 3);
        assert!((pct - 33.333333333333336).abs() < 1e-10);
    }
}
