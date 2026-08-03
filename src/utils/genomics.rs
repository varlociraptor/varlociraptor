//! genomics.rs
//!
//! Genomics utility functions.
//!
//! This module provides utilities for:
//! 1. Allele type classification (indel detection);
//! 2. Anchor length calculation (shared prefix between sequences);
//! 3. Clean indel detection (pure insertion/deletion without complex variants);
//! 4. Indel position calculation (adjusting for anchor to find true indel location);
//! 5. Sequence analysis (Svlen calculation);
//! 6. MSI status classification.

/// Check if two sequences represent an indel (different lengths).
///
/// An indel (insertion or deletion) is indicated by different
/// sequence lengths between reference and alternate alleles.
///
/// # Arguments
/// * `ref_seq` - Reference sequence
/// * `alt_seq` - Alternate sequence
///
/// # Returns
/// * `true` if lengths differ (indel)
/// * `false` if same length (SNV, MNV, or identical)
///
/// # Example
/// Insertion:
/// assert!(is_indel(b"ACAG", b"ACAGCAG"));
///
/// SNV (not an indel):
/// assert!(!is_indel(b"A", b"T"));
pub(crate) fn is_indel(ref_seq: &[u8], alt_seq: &[u8]) -> bool {
    ref_seq.len() != alt_seq.len()
}

/// Calculate anchor length (shared prefix) between two sequences.
///
/// The anchor is the longest common prefix between REF and ALT alleles
/// in VCF format. VCF indels always include at least one anchor base.
///
/// # Algorithm
/// Compares sequences byte-by-byte (case-insensitive) until mismatch.
///
/// # Arguments
/// * `ref_seq` - Reference allele sequence
/// * `alt_seq` - Alternate allele sequence
///
/// # Returns
/// Length of shared prefix in bytes
///
/// # Examples
/// Deletion: GCCT -> G, anchor = G (1 base)
/// assert_eq!(calculate_anchor_length(b"GCCT", b"G"), 1);
///
/// Insertion: G -> GCCT, anchor = G (1 base)
/// assert_eq!(calculate_anchor_length(b"G", b"GCCT"), 1);
///
/// Two anchors: TGCCT -> TG, anchor = TG (2 bases)
/// assert_eq!(calculate_anchor_length(b"TGCCT", b"TG"), 2);
pub(crate) fn calculate_anchor_length(ref_seq: &[u8], alt_seq: &[u8]) -> usize {
    let min_len = ref_seq.len().min(alt_seq.len());

    (0..min_len)
        .take_while(|&i| ref_seq[i].eq_ignore_ascii_case(&alt_seq[i]))
        .count()
}

/// Check if an indel is "clean" (pure insertion or deletion).
///
/// A clean indel must:
/// 1. Be an actual indel (different lengths)
/// 2. Have exactly ONE tail empty after removing anchor (XOR condition)
///
/// This rejects:
/// - SNVs/MNVs (same length)
/// - Complex variants (both REF and ALT have non-anchor sequence)
///
/// # Algorithm
/// 1. Check if it's an indel (different lengths)
/// 2. Find anchor (shared prefix)
/// 3. Check exactly one tail is empty (XOR)
///
/// # Arguments
/// * `ref_seq` - Reference sequence
/// * `alt_seq` - Alternate sequence
///
/// # Returns
/// * `true` if clean indel (one tail empty, lengths differ)
/// * `false` if SNV, complex variant, or identical sequences
///
/// # Examples
/// Clean deletion: GCCT -> G (anchor=G, ref_tail=CCT, alt_tail=empty)
/// assert!(is_clean_indel(b"GCCT", b"G"));
///
/// Clean insertion: G -> GCCT (anchor=G, ref_tail=empty, alt_tail=CCT)
/// assert!(is_clean_indel(b"G", b"GCCT"));
///
/// Complex: ATT -> AG (anchor=A, ref_tail=TT, alt_tail=G - BOTH non-empty)
/// assert!(!is_clean_indel(b"ATT", b"AG"));
///
/// SNV: A -> T (same length, not an indel)
/// assert!(!is_clean_indel(b"A", b"T"));
pub(crate) fn is_clean_indel(ref_seq: &[u8], alt_seq: &[u8]) -> bool {
    if !is_indel(ref_seq, alt_seq) {
        return false;
    }

    let anchor_len = calculate_anchor_length(ref_seq, alt_seq);

    let ref_tail = ref_seq.len() - anchor_len;
    let alt_tail = alt_seq.len() - anchor_len;

    (ref_tail == 0) != (alt_tail == 0)
}

/// Calculate the genomic position where an indel actually occurs.
///
/// In variant representation (VCF, etc.), the position field often points to
/// anchor base(s), not where the insertion/deletion actually happens. This
/// function calculates the position where the sequence change begins.
///
/// Returns None if the variant is not a clean indel (complex variant or SNV).
///
/// # Coordinate System
/// **This function is coordinate-system agnostic:**
/// - Input 0-based position: Output 0-based position
/// - Input 1-based position: Output 1-based position
/// - The function preserves whatever coordinate system you provide
///
/// # Algorithm
/// 1. Validate it's a clean indel (one tail empty)
/// 2. Find anchor length (shared prefix)
/// 3. Add anchor length to input position
///
/// # Arguments
/// * `pos` - Genomic position (in ANY coordinate system)
/// * `ref_seq` - Reference sequence
/// * `alt_seq` - Alternate sequence
///
/// # Returns
/// * `Some(position)` - Position where clean indel starts (same coordinate system as input)
/// * `None` - Not a clean indel (complex variant, SNV, or identical sequences)
///
/// # Examples
/// 0-based coordinates (BED):
/// Position 18630802 (0-based) with anchor G
/// Indel starts at: 18630802 + 1 = 18630803 (0-based)
/// assert_eq!(calculate_indel_position(18630802, b"GCCT", b"G"), Some(18630803));
///
/// 1-based coordinates (VCF standard):
/// Position 18630803 (1-based) with anchor G
/// Indel starts at: 18630803 + 1 = 18630804 (1-based)
/// assert_eq!(calculate_indel_position(18630803, b"GCCT", b"G"), Some(18630804));
///
/// Complex variant: returns None
/// assert_eq!(calculate_indel_position(100, b"ATT", b"AG"), None);
///
/// SNV: returns None
/// assert_eq!(calculate_indel_position(100, b"A", b"T"), None);
pub(crate) fn calculate_indel_position(pos: u64, ref_seq: &[u8], alt_seq: &[u8]) -> Option<u64> {
    if !is_clean_indel(ref_seq, alt_seq) {
        return None;
    }

    let anchor_len = calculate_anchor_length(ref_seq, alt_seq);
    Some(pos + anchor_len as u64)
}

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
    let anchor_len = calculate_anchor_length(ref_seq, alt_seq);

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
pub(crate) fn classify_msi_status(msi_score: f64, threshold: f64) -> &'static str {
    if msi_score >= threshold {
        "MSI-High"
    } else {
        "MSS"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /* ========== is_indel tests =================== */

    #[test]
    fn test_is_indel_insertions() {
        assert!(is_indel(b"A", b"ATT"));
        assert!(is_indel(b"ACAG", b"ACAGCAG"));
    }

    #[test]
    fn test_is_indel_deletions() {
        assert!(is_indel(b"GCCT", b"G"));
        assert!(is_indel(b"ATCG", b"A"));
    }

    #[test]
    fn test_is_indel_not_indels() {
        assert!(!is_indel(b"A", b"T")); // SNV
        assert!(!is_indel(b"ACG", b"TGC")); // MNV
        assert!(!is_indel(b"ACGT", b"ACGT")); // Identical
    }

    /* ====== calculate_anchor_length tests ========== */

    #[test]
    fn test_calculate_anchor_length_single() {
        assert_eq!(calculate_anchor_length(b"GCCT", b"G"), 1);
        assert_eq!(calculate_anchor_length(b"G", b"GCCT"), 1);
    }

    #[test]
    fn test_calculate_anchor_length_multiple() {
        assert_eq!(calculate_anchor_length(b"TGCCT", b"TG"), 2);
        assert_eq!(calculate_anchor_length(b"ATGCCT", b"ATG"), 3);
    }

    #[test]
    fn test_calculate_anchor_length_no_anchor() {
        assert_eq!(calculate_anchor_length(b"A", b"T"), 0);
        assert_eq!(calculate_anchor_length(b"AAA", b"TTT"), 0);
    }

    #[test]
    fn test_calculate_anchor_length_complete_match() {
        assert_eq!(calculate_anchor_length(b"ACGT", b"ACGT"), 4);
    }

    #[test]
    fn test_calculate_anchor_length_case_insensitive() {
        assert_eq!(calculate_anchor_length(b"AcGt", b"ACGT"), 4);
        assert_eq!(calculate_anchor_length(b"gcct", b"GCCT"), 4);
    }

    #[test]
    fn test_calculate_anchor_length_empty() {
        assert_eq!(calculate_anchor_length(b"", b""), 0);
        assert_eq!(calculate_anchor_length(b"ACGT", b""), 0);
        assert_eq!(calculate_anchor_length(b"", b"ACGT"), 0);
    }

    /* ======== is_clean_indel tests ================= */

    #[test]
    fn test_is_clean_indel_clean_deletion() {
        assert!(is_clean_indel(b"GCCT", b"G"));
        assert!(is_clean_indel(b"ATCG", b"A"));
    }

    #[test]
    fn test_is_clean_indel_clean_insertion() {
        assert!(is_clean_indel(b"G", b"GCCT"));
        assert!(is_clean_indel(b"A", b"ATCG"));
    }

    #[test]
    fn test_is_clean_indel_complex_variants() {
        // Both have tails - complex variant
        assert!(!is_clean_indel(b"ATT", b"AG"));
        assert!(!is_clean_indel(b"AAAGAGAGAGA", b"AAAT"));
    }

    #[test]
    fn test_is_clean_indel_same_sequence() {
        // SNV (same length) - rejected by first check
        assert!(!is_clean_indel(b"A", b"T"));
        assert!(!is_clean_indel(b"ACGT", b"ACGT"));
    }

    #[test]
    fn test_is_clean_indel_edge_cases() {
        // Empty sequences
        assert!(!is_clean_indel(b"", b""));

        // One empty (clean insertion/deletion from nothing)
        assert!(is_clean_indel(b"", b"ACGT"));
        assert!(is_clean_indel(b"ACGT", b""));
    }

    /* ==== calculate_indel_position tests =========== */

    #[test]
    fn test_calculate_indel_position_clean_deletion_0based() {
        assert_eq!(
            calculate_indel_position(18630802, b"GCCT", b"G"),
            Some(18630803)
        );
    }

    #[test]
    fn test_calculate_indel_position_clean_deletion_1based() {
        assert_eq!(
            calculate_indel_position(18630803, b"GCCT", b"G"),
            Some(18630804)
        );
    }

    #[test]
    fn test_calculate_indel_position_clean_insertion() {
        assert_eq!(calculate_indel_position(100, b"A", b"ATT"), Some(101));
    }

    #[test]
    fn test_calculate_indel_position_multiple_anchors() {
        assert_eq!(calculate_indel_position(200, b"TGCCT", b"TG"), Some(202));
    }

    #[test]
    fn test_calculate_indel_position_complex_variant_or_snv() {
        assert_eq!(calculate_indel_position(100, b"ATT", b"AG"), None);
        assert_eq!(calculate_indel_position(100, b"AAAGAGAGAGA", b"AAAT"), None);
        assert_eq!(calculate_indel_position(300, b"A", b"T"), None);
    }

    #[test]
    fn test_calculate_indel_position_no_anchor() {
        assert_eq!(calculate_indel_position(100, b"A", b"TGC"), None);
    }

    #[test]
    fn test_calculate_indel_position_chromosome_start() {
        // Edge case: Position 0 (0-based) = first chromosome base

        // Insertion at start
        assert_eq!(calculate_indel_position(0, b"A", b"ATT"), Some(1));

        // Deletion at start
        assert_eq!(calculate_indel_position(0, b"ACGT", b"A"), Some(1));
    }

    #[test]
    fn test_calculate_indel_position_no_empty_ref_seq() {
        assert_eq!(calculate_indel_position(100, b"", b"TGC"), Some(100));
    }

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
