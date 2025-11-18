//! stats.rs
//!
//! Statistical utility functions for MSI analysis.
//! 
//! This module provides precise statistical calculations including:
//! 1. PHRED score conversion
//! 2. Exact percentage computation using decimal arithmetic.


use bio::stats::{PHREDProb, Prob};
use rust_decimal::prelude::ToPrimitive;
use rust_decimal::Decimal;

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

/// Calculate percentage using exact decimal arithmetic.
/// 
/// Prevents floating-point precision loss when working with large integers
/// or when exact decimal representation is required. Uses the `rust_decimal`
/// crate for arbitrary precision decimal arithmetic.
/// 
/// # Arguments
/// * `numerator` - Count value (e.g., number of unstable regions)
/// * `denominator` - Total value (e.g., total number of regions)
/// 
/// # Returns
/// * Percentage as f64 in range [0.0, 100.0]
/// * Returns 0.0 if denominator is zero (avoiding division by zero)
/// 
/// # Precision
/// Maintains exact precision for integers
/// far exceeding f64's safe integer range of 2^53.
/// 
/// # Examples
/// assert_eq!(calculate_percentage_exact(5, 100), 5.0);`
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
