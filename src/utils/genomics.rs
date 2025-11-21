//! genomics.rs
//!
//! Genomics utility functions.
//!
//! This module provides utilities for:
//! 1. Chromosome manipulation(chromosome normalization, chromosome to rank mapping);
//! 2. Sequence analysis(Svlen calculation);
//! 3. MSI status classification.

/// Normalize chromosome names by removing "chr" prefix.
///
/// Provides consistent chromosome naming across different reference formats.
/// Both "chr1" and "1" normalize to "1" for uniform processing.
///
/// # Arguments
/// * `chrom` - Chromosome name with or without "chr" prefix
///
/// # Returns
/// Normalized chromosome name without "chr" prefix
///
/// # Note
/// Case-insensitive removal of "chr" prefix.
/// Preserves original case for the rest of the name.
///
/// # Example
/// assert_eq!(normalize_chrom("chr1"), "1");
pub(crate) fn normalize_chrom(chrom: &str) -> String {
    if chrom.len() >= 3 && chrom[..3].eq_ignore_ascii_case("chr") {
        chrom[3..].to_string()
    } else {
        chrom.to_string()
    }
}

/// Convert chromosome name to sortable rank for natural ordering.
///
/// Maps chromosome names to numeric ranks that preserve biological ordering:
/// - Autosomes: 1-22 → ranks 1-22
/// - Sex chromosomes: X → 23, Y → 24  
/// - Mitochondrial: M/MT → 25
///
/// Returns `None` for unrecognized chromosomes (e.g., scaffolds, decoys).
/// Currently supports human chromosomes only; extend for other species as needed.
///
/// # Arguments
/// * `chrom` - Chromosome name (with or without "chr" prefix)
///
/// # Returns
/// * `Some(rank)` - Numeric rank for recognized chromosomes
/// * `None` - For unrecognized chromosomes
///
/// # Note
/// Currently supports only human chromosomes (1-22, X, Y, M/MT).
/// Other species will return `None`. Extend with custom logic for other genomes.
///
/// # Example
/// assert_eq!(chrom_rank_checked("chr1"), Some(1));
/// assert_eq!(chrom_rank_checked("GL000192.1"), None);
pub(crate) fn chrom_rank_checked(chrom: &str) -> Option<u32> {
    let normalized = normalize_chrom(chrom);
    match normalized.parse::<u32>() {
        Ok(n) if (1..=22).contains(&n) => Some(n),
        _ => match normalized.as_str() {
            "X" | "x" => Some(23),
            "Y" | "y" => Some(24),
            "M" | "m" | "MT" | "Mt" | "mT" | "mt" => Some(25),
            _ => None,
        },
    }
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
pub fn classify_msi_status(msi_score: f64, threshold: f64) -> String {
    if msi_score >= threshold {
        "MSI-High".to_string()
    } else {
        "MSS".to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /* ========== normalize_chrom tests ============== */

    #[test]
    fn test_normalize_chrom() {
        assert_eq!(normalize_chrom("chr1"), "1");
        assert_eq!(normalize_chrom("1"), "1");
        assert_eq!(normalize_chrom("chrX"), "X");
        assert_eq!(normalize_chrom("MT"), "MT");
        assert_eq!(normalize_chrom("chrMT"), "MT");
        assert_eq!(normalize_chrom("Chr1"), "1");
        assert_eq!(normalize_chrom("CHR2"), "2");

        /* Upstream should be aware of these: */
        assert_eq!(normalize_chrom("GL000192.1"), "GL000192.1");
        assert_eq!(normalize_chrom(""), "");
        assert_eq!(normalize_chrom("chr"), "");
        assert_eq!(normalize_chrom("chromosome1"), "omosome1");
    }

    /* ======= chrom_rank_checked tests ============== */

    #[test]
    fn test_chrom_rank_checked_valid() {
        assert_eq!(chrom_rank_checked("22"), Some(22)); // Autosomes
        assert_eq!(chrom_rank_checked("X"), Some(23)); // Sex chromosomes
        assert_eq!(chrom_rank_checked("M"), Some(25)); // Mitochondrial

        assert_eq!(chrom_rank_checked("chr1"), Some(1)); // With chr prefix
        assert_eq!(chrom_rank_checked("chry"), Some(24)); // With chr prefix lowercase
    }

    #[test]
    fn test_chrom_rank_checked_invalid() {
        // Out of range autosomes
        assert_eq!(chrom_rank_checked("0"), None);
        assert_eq!(chrom_rank_checked("23"), None);

        // Scaffolds and decoys
        assert_eq!(chrom_rank_checked("GL000192.1"), None);

        // Other invalid inputs
        assert_eq!(chrom_rank_checked(""), None);
        assert_eq!(chrom_rank_checked("chr"), None);
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
