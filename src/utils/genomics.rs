//! genomics.rs
//!
//! Genomics utility functions.
//!
//! This module provides utilities for:
//! 1. Sequence analysis (Svlen calculation);
//! 2. MSI status classification.

/// Calculate structural variant length (SVLEN) from reference and alternate sequences.
///
/// Implements anchor-aware SVLEN calculation by identifying the common prefix
/// between reference and alternate alleles (the "anchor"), then computing the
/// difference in non-anchor sequence lengths.
///
/// # Algorithm
/// 1. Find longest common prefix (anchor) between REF and ALT
/// 2. Calculate tail lengths after anchor for both sequences
/// 3. Return: alt_tail - ref_tail = alt_len - ref_len
///
/// # Arguments
/// * `ref_seq` - Reference allele sequence
/// * `alt_seq` - Alternate allele sequence
///
/// # Returns
/// * Positive value - Insertion (ALT longer than REF)
/// * Negative value - Deletion (REF longer than ALT)
/// * Zero - Same length (likely SNV or MNV)
///
/// # Note
/// Anchor-aware: calculates length difference ignoring common prefix (anchor).
/// Handles empty sequences (start or end of sequence) gracefully.
///
/// # Examples
/// Insertion: REF=ACAG, ALT=ACAGCAG : +3
/// assert_eq!(calculate_dynamic_svlen(b"ACAG", b"ACAGCAG"), 3);
/// Deletion: REF=ACAGT, ALT=AC : -3  
/// assert_eq!(calculate_dynamic_svlen(b"ACAGT", b"AC"), -3);
pub(crate) fn calculate_dynamic_svlen(ref_seq: &[u8], alt_seq: &[u8]) -> i32 {
    // Find anchor length (longest common prefix)
    let min_len = ref_seq.len().min(alt_seq.len());
    let anchor_len = (0..min_len)
        .take_while(|&i| ref_seq[i].eq_ignore_ascii_case(&alt_seq[i]))
        .count();

    // Calculate length difference after anchor
    let ref_tail = ref_seq.len() - anchor_len;
    let alt_tail = alt_seq.len() - anchor_len;

    alt_tail as i32 - ref_tail as i32
}

/// Classify MSI status based on score and threshold.
///
/// Binary classification of microsatellite instability status:
/// - MSI-High: Score ≥ threshold
/// - MSS (Microsatellite Stable): Score < threshold
///
/// # Arguments
/// * `msi_score` - Calculated MSI score (percentage)
/// * `threshold` - Classification threshold (default 3.5%)
///
/// # Returns
/// * `"MSI-High"` - High microsatellite instability
/// * `"MSS"` - Microsatellite stable
///
/// # Examples
/// assert_eq!(classify_msi_status(5.0, 3.5), "MSI-High");
pub fn classify_msi_status(msi_score: f64, threshold: f64) -> &'static str {
    if msi_score >= threshold {
        "MSI-High"
    } else {
        "MSS"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /* ======= calculate_dynamic_svlen tests ========= */

    #[test]
    fn test_calculate_dynamic_svlen_insertions() {
        assert_eq!(calculate_dynamic_svlen(b"ACAG", b"ACAGCAG"), 3); // Simple insertion
        assert_eq!(calculate_dynamic_svlen(b"AT", b"ATATAT"), 4); // Multiple unit insertion
        assert_eq!(calculate_dynamic_svlen(b"A", b"AT"), 1); // Single base insertion
        assert_eq!(calculate_dynamic_svlen(b"", b"CAG"), 3); // No anchor insertion
        assert_eq!(calculate_dynamic_svlen(b"AAT", b"AACAG"), 2); // (Special Case)
    }

    #[test]
    fn test_calculate_dynamic_svlen_deletions() {
        assert_eq!(calculate_dynamic_svlen(b"ACAGT", b"AC"), -3); // Simple deletion
        assert_eq!(calculate_dynamic_svlen(b"ATCG", b"A"), -3); // Complete deletion after anchor
        assert_eq!(calculate_dynamic_svlen(b"AT", b"A"), -1); // Single base deletion
        assert_eq!(calculate_dynamic_svlen(b"AACAG", b"AAT"), -2); // (Special Case)
    }

    #[test]
    fn test_calculate_dynamic_svlen_substitutions() {
        assert_eq!(calculate_dynamic_svlen(b"A", b"T"), 0); // SNV
        assert_eq!(calculate_dynamic_svlen(b"ACG", b"TGC"), 0); // MNV (multiple nucleotide variant)
        assert_eq!(calculate_dynamic_svlen(b"ATCG", b"ATCG"), 0); // Same sequences
    }

    #[test]
    fn test_calculate_dynamic_svlen_case_insensitive() {
        assert_eq!(calculate_dynamic_svlen(b"acag", b"ACAGCAG"), 3);
        assert_eq!(calculate_dynamic_svlen(b"AcAgCaG", b"aCaG"), -3);
    }

    #[test]
    fn test_calculate_dynamic_svlen_edge_cases() {
        // Empty sequences
        assert_eq!(calculate_dynamic_svlen(b"", b""), 0);
        assert_eq!(calculate_dynamic_svlen(b"ATG", b""), -3);
        assert_eq!(calculate_dynamic_svlen(b"", b"ATG"), 3);

        // No common anchor
        assert_eq!(calculate_dynamic_svlen(b"AAA", b"TTT"), 0);
    }

    /* ========= classify_msi_status tests =========== */

    #[test]
    fn test_classify_msi_status() {
        assert_eq!(classify_msi_status(2.0, 3.5), "MSS"); // Below threshold
        assert_eq!(classify_msi_status(3.5, 3.5), "MSI-High"); // At threshold (inclusive)
        assert_eq!(classify_msi_status(5.0, 3.5), "MSI-High"); // Above threshold
    }
}
