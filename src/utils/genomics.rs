//! genomics.rs
//!
//! Genomics utility functions.
//!
//! This module provides utilities for:
//! 1. `MsiStatus` - the result type for MSI status classification
//! 2. Basic indel detection (length-based: are REF and ALT different lengths?)
//! 3. Anchor length calculation (shared prefix between sequences)
//! 4. Clean vs. complex indel classification (pure insertion/deletion vs. mixed)
//! 5. `reverse_alleles` - helper for reversing REF and ALT byte sequences
//! 6. Indel position calculation (adjusting for anchor to find true indel location)
//! 7. Sequence analysis (Svlen calculation)
//! 8. MSI status classification, via the `MsiStatus` result type
//!
//! Note:
//! These are generic byte-comparison utilities: they operate on any two byte
//! slices regardless of source, and make no assumptions about non-empty
//! input. Deciding what counts as meaningful data for a given caller is that
//! caller's responsibility, not this module's.
//!
//! These utilities cover the normal case throughout; unusual edge cases
//! (e.g. an indel at a contig's very first position, where VCF's own anchor
//! convention flips) are left entirely to the caller to detect and handle.

/* ============ Data Structures =================== */

/// Binary classification of microsatellite instability status.
///
/// - `High`   - displays as `"MSI-High"`
/// - `Stable` - displays as `"MSS"`
#[derive(Debug, PartialEq)]
pub enum MsiStatus {
    Stable,
    High,
}

impl MsiStatus {
    /// Returns the canonical text form of MSI status.
    pub fn as_str(&self) -> &'static str {
        match self {
            MsiStatus::Stable => "MSS",
            MsiStatus::High => "MSI-High",
        }
    }
}

impl std::fmt::Display for MsiStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/* ================================================ */

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
pub(crate) fn is_indel(ref_seq: &[u8], alt_seq: &[u8]) -> bool {
    ref_seq.len() != alt_seq.len()
}

/// Calculate anchor length (shared prefix) between two sequences.
///
/// The anchor is the longest common prefix between REF and ALT alleles.
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
pub(crate) fn is_clean_indel(ref_seq: &[u8], alt_seq: &[u8]) -> bool {
    if !is_indel(ref_seq, alt_seq) {
        return false;
    }

    let anchor_len = calculate_anchor_length(ref_seq, alt_seq);

    let ref_tail = ref_seq.len() - anchor_len;
    let alt_tail = alt_seq.len() - anchor_len;

    (ref_tail == 0) != (alt_tail == 0)
}

/// Reverse both allele sequences.
///
/// # Note
/// A caller handling the position-1 trailing-anchor edge case (see this
/// module's docs) can run this and then use `calculate_anchor_length`/
/// `is_clean_indel` unmodified on the result: a trailing anchor becomes a
/// leading one under reversal. Sufficient for structural checks (which only
/// depend on tail lengths, not byte order); not sufficient if a caller needs
/// the changed content in correct real-world order.
///
/// # Arguments
/// * `ref_seq` - Reference sequence
/// * `alt_seq` - Alternate sequence
///
/// # Returns
/// `(reversed ref_seq, reversed alt_seq)`
pub(crate) fn reverse_alleles(ref_seq: &[u8], alt_seq: &[u8]) -> (Vec<u8>, Vec<u8>) {
    (
        ref_seq.iter().rev().copied().collect(),
        alt_seq.iter().rev().copied().collect(),
    )
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
pub(crate) fn calculate_indel_position(pos: u64, ref_seq: &[u8], alt_seq: &[u8]) -> Option<u64> {
    if !is_clean_indel(ref_seq, alt_seq) {
        return None;
    }

    let anchor_len = calculate_anchor_length(ref_seq, alt_seq);
    Some(pos + anchor_len as u64)
}

/// Calculate structural variant length (SVLEN) from reference and alternate sequences.
///
/// Returns the length difference between the alternate and reference alleles.
/// The shared prefix (anchor) cancels out of the subtraction, so the result
/// is simply `alt_len - ref_len`.
///
/// # Arguments
/// * `ref_seq` - Reference allele sequence
/// * `alt_seq` - Alternate allele sequence
///
/// # Returns
/// * Positive value - Insertion (ALT longer than REF)
/// * Negative value - Deletion (REF longer than ALT)
/// * Zero - Same length (likely SNV or MNV)
pub(crate) fn calculate_dynamic_svlen(ref_seq: &[u8], alt_seq: &[u8]) -> i32 {
    alt_seq.len() as i32 - ref_seq.len() as i32
}

/// Classify MSI status based on score and threshold.
///
/// Binary classification of microsatellite instability status:
/// - `MsiStatus::High`: Score ≥ threshold (displays as "MSI-High")
/// - `MsiStatus::Stable`: Score < threshold (displays as "MSS")
///
/// # Arguments
/// * `msi_score` - Calculated MSI score (percentage)
/// * `threshold` - Classification threshold
///
/// # Returns
/// * `MsiStatus::High` - High microsatellite instability
/// * `MsiStatus::Stable` - Microsatellite stable
pub(crate) fn classify_msi_status(msi_score: f32, threshold: f32) -> MsiStatus {
    if msi_score >= threshold {
        MsiStatus::High
    } else {
        MsiStatus::Stable
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /* ========== MsiStatus tests ================== */

    #[test]
    fn test_msi_status_as_str() {
        assert_eq!(MsiStatus::Stable.as_str(), "MSS");
        assert_eq!(MsiStatus::High.as_str(), "MSI-High");
    }

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

    /* ========== reverse_alleles tests ============== */

    #[test]
    fn test_reverse_alleles() {
        // Non-palindromic: verifies bytes are actually reversed, not a no-op.
        assert_eq!(
            reverse_alleles(b"ACGT", b"TT"),
            (b"TGCA".to_vec(), b"TT".to_vec())
        );
        // Palindromic: the realistic indel edge case (position-1 trailing anchor).
        assert_eq!(
            reverse_alleles(b"GAGAGAG", b"G"),
            (b"GAGAGAG".to_vec(), b"G".to_vec())
        );
    }

    #[test]
    fn test_reverse_alleles_empty() {
        assert_eq!(reverse_alleles(b"", b""), (vec![], vec![]));
    }

    /* ==== calculate_indel_position tests =========== */

    #[test]
    fn test_calculate_indel_position_clean_deletion() {
        assert_eq!(
            calculate_indel_position(18630802, b"GCCT", b"G"),
            Some(18630803)
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
    fn test_calculate_indel_position_empty_ref_seq_yields_input_pos() {
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
        assert_eq!(calculate_dynamic_svlen(b"AAA", b"TTTTT"), 2);
    }

    /* ========= classify_msi_status tests =========== */

    #[test]
    fn test_classify_msi_status() {
        assert_eq!(classify_msi_status(2.0, 3.5), MsiStatus::Stable); // Below threshold
        assert_eq!(classify_msi_status(3.5, 3.5), MsiStatus::High); // At threshold (inclusive)
        assert_eq!(classify_msi_status(5.0, 3.5), MsiStatus::High); // Above threshold
    }
}
