//! genomics.rs
//!
//! Utility functions for genomics-related tasks

/// Strip "chr" prefix, so both "chr1" and "1" normalize to "1".
pub(crate) fn normalize_chrom(chrom: &str) -> String {
    chrom.trim_start_matches("chr").to_string()
}

/// Convert chromosome name to a sortable rank (natural human order)
/// This caters for chromosomes 1-22, X, Y, M/MT only at the moment,
/// as these are the only ones relevant for Human MSI analysis.
/// Note: This can be extended later if needed.
pub(crate) fn chrom_rank_checked(chrom: &str) -> Option<u32> {
    let s = normalize_chrom(chrom);
    match s.parse::<u32>() {
        Ok(n @ 1..=22) => Some(n),
        Ok(_) => None,
        Err(_) => match s.as_str() {
            "X" => Some(23),
            "Y" => Some(24),
            "M" | "MT" => Some(25),
            _ => None,
        },
    }
}

/// Calculate SVLEN dynamically using anchor detection
/// Positive = insertion, Negative = deletion
/// Example: REF=ACAG, ALT=ACAGCAG -> SVLEN=3
pub(crate) fn calculate_dynamic_svlen(ref_seq: &[u8], alt_seq: &[u8]) -> i32 {
    let mut anchor_len = 0;
    let min_len = ref_seq.len().min(alt_seq.len());

    for i in 0..min_len {
        if ref_seq[i].eq_ignore_ascii_case(&alt_seq[i]) {
            anchor_len += 1;
        } else {
            break;
        }
    }

    if alt_seq.len() > ref_seq.len() {
        (alt_seq.len() - anchor_len) as i32
    } else {
        -((ref_seq.len() - anchor_len) as i32)
    }
}

/// Classify MSI status based on score and threshold
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

    #[test]
    fn test_normalize_chrom() {
        assert_eq!(normalize_chrom("chr1"), "1");
        assert_eq!(normalize_chrom("1"), "1");
        assert_eq!(normalize_chrom("chrX"), "X");
        assert_eq!(normalize_chrom("MT"), "MT");
    }

    #[test]
    fn test_chrom_rank_checked_valid() {
        assert_eq!(chrom_rank_checked("1"), Some(1));
        assert_eq!(chrom_rank_checked("22"), Some(22));
        assert_eq!(chrom_rank_checked("chrX"), Some(23));
        assert_eq!(chrom_rank_checked("chrY"), Some(24));
        assert_eq!(chrom_rank_checked("MT"), Some(25));
    }

    #[test]
    fn test_chrom_rank_checked_invalid() {
        assert_eq!(chrom_rank_checked("chr23"), None);
        assert_eq!(chrom_rank_checked("GL000192.1"), None);
        assert_eq!(chrom_rank_checked("random"), None);
    }

    #[test]
    fn test_calculate_dynamic_svlen_insertion() {
        let ref_seq = b"ACAG";
        let alt_seq = b"ACAGCAG";
        assert_eq!(calculate_dynamic_svlen(ref_seq, alt_seq), 3);
    }

    #[test]
    fn test_calculate_dynamic_svlen_deletion() {
        let ref_seq = b"ACAGT";
        let alt_seq = b"AC";
        assert_eq!(calculate_dynamic_svlen(ref_seq, alt_seq), -3);
    }

    #[test]
    fn test_calculate_dynamic_svlen_case_insensitive() {
        let ref_seq = b"acag";
        let alt_seq = b"ACAGCAG";
        assert_eq!(calculate_dynamic_svlen(ref_seq, alt_seq), 3);
    }

    #[test]
    fn test_classify_msi_status() {
        assert_eq!(classify_msi_status(5.0, 3.5), "MSI-High");
        assert_eq!(classify_msi_status(2.0, 3.5), "MSS");
        assert_eq!(classify_msi_status(3.5, 3.5), "MSI-High");
    }
}
