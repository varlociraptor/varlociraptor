//! variant_analysis.rs
//!
//! Variant Analysis Utilities for MSI preprocessing.
//!
//! This module provides:
//! 1. `is_perfect_repeat` function to determine if an indel is a perfect tandem repeat of a microsatellite motif.
//! 2. `should_include_variant` function to analyze a variant-allele pair and determine if it should be included
//!    in the preprocessed output based on MSI relevance and perfect repeat status.
//!

use anyhow::Result;
use log::debug;
use rust_htslib::bcf;

use crate::utils::bcf_utils::{
    get_chrom, is_breakend, is_reference_allele, is_spanning_deletion, is_symbolic,
};
use crate::utils::genomics::{calculate_anchor_length, calculate_dynamic_svlen, is_indel};
use crate::utils::ms_bed::BedRegion;

/* ============ Data Structures =================== */

/// Classification of indel repeat pattern relative to microsatellite motif.
#[derive(Debug, PartialEq)]
enum RepeatStatus {
    /// Variant indel is a perfect tandem repeat of the motif.
    Perfect,
    /// Variant indel does not match motif pattern or fails validation.
    NA,
}

/* ================================================ */

/// Check if an indel is a perfect tandem repeat of a microsatellite motif.
///
/// Determines if the changed sequence (insertion or deletion) consists
/// entirely of complete motif units.
///
/// # Algorithm
/// 1. Find anchor (common prefix between REF and ALT)
/// 2. Extract changed sequence after anchor
/// 3. Verify changed sequence length matches |SVLEN|
/// 4. Check if changed sequence is exact motif repeats
///
/// # Arguments
/// * `alt_seq` - Alternate allele sequence
/// * `svlen` - Structural variant length (positive=insertion, negative=deletion)
/// * `motif` - Microsatellite motif to match
/// * `ref_seq` - Reference allele sequence
///
/// # Returns
/// * `RepeatStatus::Perfect` - Changed sequence is exact motif repeats
/// * `RepeatStatus::NA` - Does not match or fails validation
///
/// # Note
/// Small indels (even a single-base insertion) count as MS indels.
/// Slippage acts per repeat unit; for mononucleotide repeats that's 1bp,
/// and it's a common real MSI event.
/// https://www.sciencedirect.com/topics/neuroscience/slipped-strand-mispairing
/// https://en.wikipedia.org/wiki/Microsatellite
///
/// No upper size limit is applied to the changed sequence here - any
/// perfect-motif match is valid regardless of length, no literature
/// we know establishes a cutoff. Large deletions are already excluded
/// before reaching this function (should_include_variant filters symbolic
/// alleles via is_symbolic; Varlociraptor converts deletions >50bp to
/// symbolic <DEL> upstream), so no separate cap is needed here.
///
/// We discard insertions that are not perfect multiples of the motif, and
/// this is the right approach. Replication slippage, the mechanism
/// underlying MSI, produces exclusively whole repeat-unit changes -
/// Brinkmann et al. 1998 observed zero incomplete-repeat events across
/// several germline STR mutations, concluding replication slippage was
/// the sole mechanism. A CAGCAA insertion therefore cannot arise from
/// slippage alone; it would require either a subsequent point mutation
/// (a separate mutational process) or a different insertion mechanism
/// altogether. Both scenarios are computationally indistinguishable,
/// introducing unresolvable ambiguity into the MSI score. This could be
/// an interesting future direction.
/// https://pmc.ncbi.nlm.nih.gov/articles/instance/1377148/pdf/9585597.pdf
///
/// We do not track net tract-length change after an indel (i.e. whether
/// the region "stops being" a microsatellite post-variant). MSI is defined
/// as a repeat-count change relative to reference, not a constraint on
/// resulting tract length - the BED regions are already MS loci by
/// definition from their source annotation.
///
/// # Example
/// assert_eq!(is_perfect_repeat(b"ACAGCAG", 3, "CAG", b"ACAG"), RepeatStatus::Perfect);
fn is_perfect_repeat(alt_seq: &[u8], svlen: i32, motif: &str, ref_seq: &[u8]) -> RepeatStatus {
    // 0. Handling Edge Cases
    // svlen == 0 is not checked separately - is_clean_indel (guaranteed
    // upstream) already ensures ref_seq.len() != alt_seq.len(), and svlen
    // is always self-computed as that exact difference, so svlen != 0
    // whenever this function is reached.
    if ref_seq.is_empty() || alt_seq.is_empty() || !ref_seq[0].eq_ignore_ascii_case(&alt_seq[0]) {
        return RepeatStatus::NA;
    }

    // 1. Use genomics utility to check if clean indel
    // Relies on caller guarantee: is_clean_indel already checked upstream
    // in variant_overlaps_region (intersection.rs), same (ref_seq, alt_seq).

    // 2. Finding the anchor length and absolute SVLEN
    // Note: Anchor length 0 is not errored as a valid indel
    let abs_svlen = svlen.unsigned_abs() as usize;
    let anchor_len = calculate_anchor_length(ref_seq, alt_seq);

    // 3. Extracting the changed sequence
    // Relies on caller guarantee: svlen is always correctly-computed,
    // see should_include_variant. Combined with is_clean_indel check upstream,
    // as mentioned in (1), this guarantees anchor_len is always < the relevant seq length,
    // so the else branches for that are unreachable and therefore not included.
    let changed_seq = if svlen > 0 {
        &alt_seq[anchor_len..]
    } else {
        &ref_seq[anchor_len..]
    };

    // Validate: changed sequence length should match SVLEN
    if changed_seq.len() != abs_svlen {
        return RepeatStatus::NA;
    }

    // 4. Check if changed sequence is a perfect repeat of the motif
    let motif_bytes: Vec<u8> = motif.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let motif_len = motif_bytes.len();

    // Guaranteed non-empty by ms_bed.rs's parse_motif_from_name. Re-enable if that changes.
    // if motif_len == 0 {
    //     return RepeatStatus::NA
    // }

    if changed_seq.len() % motif_len != 0 {
        return RepeatStatus::NA;
    }

    for (i, &base) in changed_seq.iter().enumerate() {
        let expected_base = motif_bytes[i % motif_len];
        if base.to_ascii_uppercase() != expected_base {
            return RepeatStatus::NA;
        }
    }

    RepeatStatus::Perfect
}

/// Analyze a variant if it should be included in preprocessed output.
///
/// This function processes a single variant-allele pair to determine if it's
/// relevant for preprocessed output. It performs the following steps:
///
/// # Algorithm
/// 1. **Filtering**: Skip variants that aren't relevant for MSI:
///    - Reference alleles (ALT=<REF>)
///    - Symbolic alleles (<DEL>, <INS>)
///    - Breakends (complex structural variants)
///    - Spanning deletions (*)
///    - Non-indel variants
/// 2. **Check Perfect Repeat Status**:
///    - Calculate SVLEN (indel length)
///    - Verify indel is perfect tandem repeat of motif
///
/// # Arguments
/// * `record` - BCF record representing the variant
/// * `alt_idx` - Index of the alternate allele to analyze
/// * `region` - BED region context for motif information
///
/// # Returns
/// * `Ok(true)` - Variant is a perfect MS indel at this locus
/// * `Ok(false)` - Variant should be skipped (not relevant for MSI quantification)
/// * `Err` - Error reading variant data
///
/// # Example
/// assert!(should_include_variant(&record, 0, &region).unwrap());
pub(super) fn should_include_variant(
    record: &bcf::Record,
    alt_idx: usize,
    region: &BedRegion,
) -> Result<bool> {
    let alleles = record.alleles();
    let ref_allele = alleles[0];
    let alt_allele = alleles[alt_idx + 1]; // +1 because alleles[0] is REF

    /* 1. Filter non indel variants */
    if is_reference_allele(alt_allele)
        || is_symbolic(alt_allele)
        || is_breakend(alt_allele)
        || is_spanning_deletion(alt_allele)
        || !is_indel(ref_allele, alt_allele)
    {
        debug!(
            "Filtering non-indel variant at {}:{}",
            get_chrom(record)?,
            record.pos() + 1
        );
        return Ok(false);
    }

    /* 2. Check Perfect Repeat Status */
    let svlen = calculate_dynamic_svlen(ref_allele, alt_allele);
    let repeat_status = is_perfect_repeat(alt_allele, svlen, &region.motif, ref_allele);

    Ok(repeat_status == RepeatStatus::Perfect)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::bcf_utils::tests::{create_minimal_vcf, read_first_record};
    use crate::utils::ms_bed::BedRegion;

    /* ========== is_perfect_repeat tests ============ */

    #[test]
    fn test_is_perfect_repeat_insertion() {
        assert_eq!(
            is_perfect_repeat(b"ACAGCAG", 3, "CAG", b"ACAG"),
            RepeatStatus::Perfect
        );
        assert_eq!(
            is_perfect_repeat(b"ACAGCAGCAG", 6, "CAG", b"ACAG"),
            RepeatStatus::Perfect
        );
        assert_eq!(
            is_perfect_repeat(b"ATT", 1, "T", b"AT"),
            RepeatStatus::Perfect
        );
        // Small indels (1bp mononucleotide) are valid MSI events - see is_perfect_repeat's doc comment.
    }

    #[test]
    fn test_is_perfect_repeat_deletion() {
        assert_eq!(
            is_perfect_repeat(b"ACAG", -3, "CAG", b"ACAGCAG"),
            RepeatStatus::Perfect
        );
    }

    #[test]
    fn test_is_perfect_repeat_case_insensitive() {
        assert_eq!(
            is_perfect_repeat(b"acagCAG", 3, "cag", b"acag"),
            RepeatStatus::Perfect
        );
    }

    #[test]
    fn test_is_perfect_repeat_case_insensitive_first_byte() {
        assert_eq!(
            is_perfect_repeat(b"AcagCAG", 3, "cag", b"acag"),
            RepeatStatus::Perfect
        );
    }

    #[test]
    fn test_is_perfect_repeat_not_perfect() {
        assert_eq!(
            is_perfect_repeat(b"ACAGCAT", 3, "CAG", b"ACAG"),
            RepeatStatus::NA
        );
        assert_eq!(
            is_perfect_repeat(b"ACAGCA", 2, "CAG", b"ACAG"),
            RepeatStatus::NA
        );
        assert_eq!(
            is_perfect_repeat(b"AAAGAGAGAGA", 7, "GA", b"AAAT"),
            RepeatStatus::NA
        );
        // Caught by the length check - is_clean_indel is guaranteed upstream, not re-checked here.
    }

    #[test]
    fn test_is_perfect_repeat_edge_cases() {
        assert_eq!(is_perfect_repeat(b"CAG", 3, "CAG", b""), RepeatStatus::NA);
        assert_eq!(is_perfect_repeat(b"", 0, "CAG", b"CAG"), RepeatStatus::NA);
        assert_eq!(
            is_perfect_repeat(b"ACAT", 0, "CAG", b"ACAG"),
            RepeatStatus::NA
        );
        assert_eq!(
            is_perfect_repeat(b"TCAG", 3, "TCAG", b"A"),
            RepeatStatus::NA
        );
    }

    /* ======= should_include_variant tests ========== */

    #[test]
    fn test_should_include_variant_filters_snv() {
        let tmp_vcf = create_minimal_vcf(
            &[br"##contig=<ID=chr1,length=1000000>"],
            &[(0, 99, b"A", &[b"T"])],
        );

        let record = read_first_record(tmp_vcf.path());

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 200,
            motif: "A".to_string(),
        };

        let result = should_include_variant(&record, 0, &region).unwrap();
        assert!(!result); // SNV should be filtered out
    }

    #[test]
    fn test_should_include_variant_perfect_indel() {
        let tmp_vcf = create_minimal_vcf(
            &[br"##contig=<ID=chr1,length=1000000>"],
            &[(0, 99, b"ACAG", &[b"ACAGCAG"])],
        );

        let record = read_first_record(tmp_vcf.path());

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 100,
            end: 121,
            motif: "CAG".to_string(),
        };

        let result = should_include_variant(&record, 0, &region).unwrap();

        assert!(result); // Perfect indel should be included
    }

    #[test]
    fn test_should_include_variant_multi_allelic() {
        let tmp_vcf = create_minimal_vcf(
            &[br"##contig=<ID=chr1,length=1000000>"],
            &[(0, 99, b"A", &[b"T", b"ATG"])],
        );

        let record = read_first_record(tmp_vcf.path());

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 200,
            motif: "TG".to_string(),
        };

        // alt_idx=0 is SNV (T) - should filter out
        assert!(!should_include_variant(&record, 0, &region).unwrap());

        // alt_idx=1 is indel (ATG) - should include
        assert!(should_include_variant(&record, 1, &region).unwrap());
    }

    #[test]
    fn test_should_include_variant_filters_symbolic() {
        let tmp_vcf = create_minimal_vcf(
            &[br"##contig=<ID=chr1,length=1000000>"],
            &[(0, 99, b"A", &[b"<DEL>"])],
        );

        let record = read_first_record(tmp_vcf.path());

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 200,
            motif: "A".to_string(),
        };

        let result = should_include_variant(&record, 0, &region).unwrap();
        assert!(!result); // Symbolic allele should be filtered out
    }

    #[test]
    fn test_should_include_variant_imperfect_repeat() {
        let tmp_vcf = create_minimal_vcf(
            &[br"##contig=<ID=chr1,length=1000000>"],
            &[(0, 99, b"ACAG", &[b"ACAGCAT"])], // Not perfect CAG repeat
        );

        let record = read_first_record(tmp_vcf.path());

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 100,
            end: 121,
            motif: "CAG".to_string(),
        };

        let result = should_include_variant(&record, 0, &region).unwrap();
        assert!(!result); // Imperfect repeat should be filtered
    }
}
