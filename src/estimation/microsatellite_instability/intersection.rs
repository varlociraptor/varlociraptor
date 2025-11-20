//! intersection.rs
//!
//! Streaming intersection and variant analysis for MSI detection.
//!
//! This module provides:
//! 1. Streaming intersection of BED regions with VCF variants
//! 2. Perfect microsatellite repeat detection
//! 3. Variant filtering and probability extraction
//!
//! The streaming approach processes BED regions sequentially while maintaining
//! a sliding window of VCF variants, enabling memory-efficient analysis of
//! large datasets.

use std::collections::{HashMap, VecDeque};
use std::path::PathBuf;

use anyhow::{Context, Result};
use bio::io::bed;
use log::{debug, info};
use rust_htslib::bcf::{self, header::HeaderView, Read};

use crate::errors::Error;
use crate::utils::bcf_utils::{
    get_chrom, get_prob_absent, get_sample_afs, get_svlen, is_breakend, is_indel,
    is_reference_allele, is_spanning_deletion, is_symbolic,
};
use crate::utils::genomics::chrom_rank_checked;
use crate::utils::is_phred_scaled;
use crate::utils::ms_bed::{parse_bed_record, BedRegion};

/* ============ Data Structures =================== */

/// Classification of indel repeat pattern relative to microsatellite motif.
#[derive(Debug, Clone, PartialEq)]
enum RepeatStatus {
    /// Variant indel is a perfect tandem repeat of the motif.
    Perfect,
    /// Variant indel does not match motif pattern or fails validation.
    NA,
}

/// Analyzed variant passing all filters for MSI analysis.
///
/// Contains probability and allele frequency data extracted
/// from a VCF record that passed perfect repeat validation.
#[derive(Debug, Clone)]
pub(super) struct Variant {
    /// Combined probability variant is absent (PROB_ABSENT + PROB_ARTIFACT)
    pub(super) prob_absent: f64,
    /// Per-sample allele frequencies from FORMAT:AF
    pub(super) sample_afs: HashMap<String, f64>,
}

/// Summary of variants found within a single BED region.
#[derive(Debug)]
pub(super) struct RegionSummary {
    /// Variants passing all filters at intersection step within this region
    pub(super) variants: Vec<Variant>,
}

impl RegionSummary {
    /// Add a variant to the region summary.
    fn add_variant(&mut self, variant: Variant) {
        self.variants.push(variant);
    }
}

/// Statistics collected during intersection process.
#[derive(Debug, Clone)]
pub(super) struct IntersectionStats {
    /// Total number of regions processed from BED file.
    pub total_regions: usize,
    /// Number of regions skipped due to invalid motif length.
    pub skipped_invalid_motif: usize,
}

impl IntersectionStats {
    /// Calculate number of valid regions processed.
    fn valid_regions(&self) -> usize {
        self.total_regions - self.skipped_invalid_motif
    }

    /// Log a summary of intersection statistics.
    fn log_summary(&self) {
        info!("Intersection complete:");
        info!("    - Total regions in BED: {}", self.total_regions);
        info!(
            "    Valid regions (1-6 motif length): {}",
            self.valid_regions()
        );
        info!(
            "    - Skipped (invalid motif): {}",
            self.skipped_invalid_motif
        );
    }
}

/* ================================================ */

/* ======== Variant Analysis Functions ============ */

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
/// # Example
/// assert_eq!(is_perfect_repeat(b"ACAGCAG", 3, "CAG", b"ACAG"), RepeatStatus::Perfect);
fn is_perfect_repeat(alt_seq: &[u8], svlen: i32, motif: &str, ref_seq: &[u8]) -> RepeatStatus {
    // 0. Handling Edge Cases
    if ref_seq.is_empty() || alt_seq.is_empty() || ref_seq[0] != alt_seq[0] || svlen == 0 {
        return RepeatStatus::NA;
    }

    // 1. Finding the anchor
    // Note: Anchor length 0 is not errored as a valid indel
    let abs_svlen = svlen.unsigned_abs() as usize;
    let mut anchor_len = 0;
    let min_len = ref_seq.len().min(alt_seq.len());

    for i in 0..min_len {
        if ref_seq[i].eq_ignore_ascii_case(&alt_seq[i]) {
            anchor_len += 1;
        } else {
            break;
        }
    }

    // 2. Making sure is a clean indel, ignoring complex indels
    let ref_tail = ref_seq.len() - anchor_len;
    let alt_tail = alt_seq.len() - anchor_len;

    // Require a clean indel: one tail must be empty
    // This rejects complex variants like AAT -> AAGAGAGAGA
    if ref_tail != 0 && alt_tail != 0 {
        return RepeatStatus::NA;
    }

    // 3. Extracting the changed sequence
    let changed_seq = if svlen > 0 {
        if anchor_len < alt_seq.len() {
            &alt_seq[anchor_len..]
        } else {
            return RepeatStatus::NA;
        }
    } else if anchor_len < ref_seq.len() {
        &ref_seq[anchor_len..]
    } else {
        return RepeatStatus::NA;
    };

    // Validate: changed sequence length should match SVLEN
    if changed_seq.len() != abs_svlen || changed_seq.is_empty() {
        return RepeatStatus::NA;
    }

    // 4. Check if changed sequence is a perfect repeat of the motif
    let motif_bytes: Vec<u8> = motif.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let motif_len = motif_bytes.len();

    /* NOTE: In case we toggle the Error on motif.len() in bed parsing off, turn this onn.*/
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

/// Analyze a variant within a microsatellite region
///
/// This function processes a single variant-allele pair to determine if it's
/// relevant for MSI analysis. It performs the following steps:
///
/// # Algorithm
/// 1. **Filtering**: Skip variants that aren't relevant for MSI:
///    - SNVs (single nucleotide variants)
///    - Symbolic alleles (<DEL>, <INS>)
///    - Breakends (complex structural variants)
///    - Spanning deletions (*)
///    - Non-indel variants
/// 2. **SVLEN Calculation**: Determine the length of insertion/deletion
/// 3. **Perfect Repeat Check**: Determine if indel is a perfect repeat
/// 4. **Probability Extraction**: Get PROB_ABSENT from INFO field
/// 5. **Allele Frequency Extraction**: Get per-sample AFs from FORMAT:AF
/// 6. **Final Validation**: Ensure all required data is present for perfect repeats
///    - If perfect repeat but missing data(prob_absent/sample_afs) : (skip variant)
///
/// # Arguments
/// * `record` - BCF record representing the variant
/// * `header` - BCF header for metadata access
/// * `alt_idx` - Index of the alternate allele to analyze
/// * `region` - BED region context for motif information
/// * `samples_index_map` - Map of sample names to their indices in the VCF
/// * `is_phred` - Whether probabilities are PHRED-scaled
///
/// # Returns
/// * `Ok(None)` - Variant should be skipped (SNV, symbolic, etc.)
/// * `Ok(Some(VariantAnalysis))` - Valid indel with extracted data
/// * `Err(...)` - Error occurred (invalid data, I/O error, etc.)
///
/// # Example
/// assert!(analyze_variant(&record, &header, 0, &region, &samples_index_map, false).is_ok());
fn analyze_variant(
    record: &bcf::Record,
    header: &HeaderView,
    alt_idx: usize,
    region: &BedRegion,
    samples_index_map: &HashMap<String, usize>,
    is_phred: bool,
) -> Result<Option<Variant>> {
    let alleles = record.alleles();
    let ref_allele = alleles[0];
    let alt_allele = alleles[alt_idx + 1]; // +1 because alleles[0] is REF

    if is_reference_allele(alt_allele)
        || is_symbolic(alt_allele)
        || is_breakend(alt_allele)
        || is_spanning_deletion(alt_allele)
        || !is_indel(ref_allele, alt_allele)
    {
        debug!(
            "Skipping filter variant at {}:{}",
            get_chrom(record, header)?,
            record.pos() + 1
        );
        return Ok(None);
    }

    let svlen = get_svlen(record, alt_idx, ref_allele, alt_allele)?;
    let repeat_status = is_perfect_repeat(alt_allele, svlen, &region.motif, ref_allele);
    let prob_absent = get_prob_absent(record, header, alt_idx, is_phred)?;
    let sample_afs = get_sample_afs(record, header, samples_index_map, alt_idx)?;

    if repeat_status == RepeatStatus::Perfect {
        if prob_absent.is_none() || sample_afs.is_empty() {
            debug!(
                "Perfect repeat at {}:{} missing required probability or all_afs - downgrading to NA",
                get_chrom(record, header)?,
                record.pos() + 1
            );
            return Ok(None);
        }
    } else {
        return Ok(None);
    }

    Ok(Some(Variant {
        prob_absent: prob_absent.unwrap(),
        sample_afs,
    }))
}

/* ================================================ */

/* ============ Streaming Intersection ============ */

/// Check if a variant position is within a BED region
/// Uses point-based overlap (variant position only)
///
/// # Arguments
/// * `record` - BCF record representing the variant
/// * `region` - BED region to check overlap against
///
/// # Returns
/// * `true` if variant position is within region, else `false`
///
/// # Note: Point-based overlap means we only consider the variant position,
///
/// # Example
/// assert!(variant_overlaps_region(&record, &region));
#[inline]
fn variant_overlaps_region(record: &bcf::Record, region: &BedRegion) -> bool {
    // Variant position (0-based)
    let pos = record.pos() as u64;
    // BED region (0-based, half-open [start, end))
    pos >= region.start && pos < region.end
}

/// Perform streaming intersection of BED regions with VCF variants.
///
/// Uses a sliding window approach to efficiently match variants to regions
/// without loading entire files into memory.
///
/// # Algorithm
/// 1. For each BED region (in genomic order):
///    - Remove variants before this region from window
///    - Load variants until past this region
///    - Analyze overlapping variants for perfect repeats
///
/// # Arguments
/// * `bed_path` - Path to BED file with microsatellite regions
/// * `vcf_path` - Path to VCF/BCF file with variant calls
/// * `samples_index_map` - Map of sample names to VCF indices
///
/// # Returns
/// * `Vec<RegionSummary>` - Regions with valid variant
/// * `usize` - Total number of valid BED regions processed
///
/// # Errors
/// Returns error if:
/// - Files cannot be opened
/// - No chromosome overlap between BED and VCF
/// - Record parsing fails
///
/// # Example
/// assert!(intersect_streaming(&bed_path, &vcf_path, &samples_index_map).is_ok());
pub(super) fn intersect_streaming(
    bed_path: &PathBuf,
    vcf_path: &PathBuf,
    samples_index_map: &HashMap<String, usize>,
) -> Result<(Vec<RegionSummary>, usize)> {
    /* ========== Setup: Initialize readers and state ========== */
    let mut bed_reader = bed::Reader::from_file(bed_path).context("Failed to open BED file")?;
    let mut vcf = bcf::Reader::from_path(vcf_path).context("Failed to open VCF file")?;
    let header: HeaderView = vcf.header().clone();

    let is_phred = is_phred_scaled(&vcf);
    info!(
        "Probabilities in VCF/BCF are {} scaled",
        if is_phred { "PHRED" } else { "linear" }
    );

    // Results & Stats
    let mut results = Vec::new();
    let mut total_regions = 0;
    let mut skipped_invalid_region = 0;
    let mut variant_window: VecDeque<(bcf::Record, String)> = VecDeque::new();
    let mut seen_any_chrom_overlap = false;

    /* ======================================================== */
    /* ========== Main Loop: Process each BED region ========== */
    for (line_num, bed_result) in bed_reader.records().enumerate() {
        let bed_record = bed_result.map_err(|e| Error::BedRecordReadFailed {
            line: line_num + 1,
            details: e.to_string(),
        })?;
        let region = parse_bed_record(&bed_record)?;
        let region_rank = chrom_rank_checked(&region.chrom);

        total_regions += 1;

        if !region.is_valid_motif() || region_rank.is_none() {
            skipped_invalid_region += 1;
            debug!(
                " Skipping region {} with invalid motif length {} or chromosome: {}",
                region.region_id(),
                region.motif_length(),
                region.chrom
            );
            continue;
        }
        let region_rank = region_rank.unwrap();

        /* ------------------------------------------------------------- */
        /* - 1. Window management: Remove variants before this region -- */
        while let Some((record, chrom)) = variant_window.front() {
            let chrom_rank = chrom_rank_checked(chrom);
            if chrom_rank.is_none() {
                debug!(
                    " Skipping VCF record with unsupported chromosome: {}",
                    chrom
                );
                variant_window.pop_front();
                continue;
            }
            let chrom_rank = chrom_rank.unwrap();

            if chrom_rank < region_rank {
                variant_window.pop_front();
            } else if chrom_rank == region_rank {
                let pos = record.pos() as u64;
                if pos < region.start {
                    variant_window.pop_front();
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        /* ------------------------------------------------------------- */
        /* ------------------------------------------------------------- */

        /* ------------------------------------------------------------- */
        /* -2. Window management: Load variants until we pass region --- */
        loop {
            if let Some((record, chrom)) = variant_window.back() {
                let chrom_rank = chrom_rank_checked(chrom);
                if chrom_rank.is_none() {
                    debug!(
                        " Skipping VCF record with unsupported chromosome: {}",
                        chrom
                    );
                    variant_window.pop_back();
                    continue;
                }
                let chrom_rank = chrom_rank.unwrap();

                if chrom_rank > region_rank {
                    break;
                } else if chrom_rank == region_rank {
                    let pos = record.pos() as u64;
                    if pos >= region.end {
                        break;
                    }
                }
            }

            let mut next_record = vcf.empty_record();
            match vcf.read(&mut next_record) {
                None => break, // EOF
                Some(Err(e)) => {
                    return Err(Error::VcfRecordReadFailed {
                        details: e.to_string(),
                    }
                    .into());
                }
                Some(Ok(())) => {
                    let chrom = get_chrom(&next_record, &header)?;
                    let pos = next_record.pos();

                    if pos < 0 {
                        debug!(
                            " Skipping malformed VCF record with invalid position: {}:{}",
                            chrom, pos
                        );
                        continue;
                    }

                    variant_window.push_back((next_record, chrom));
                }
            }
        }
        /* ------------------------------------------------------------- */
        /* ------------------------------------------------------------- */

        /* ------------------------------------------------------------- */
        /* - 3. Process variants overlapping this region --------------- */
        let mut region_summary: Option<RegionSummary> = None;

        for (record, chrom) in &variant_window {
            let Some(chrom_rank) = chrom_rank_checked(chrom) else {
                continue;
            };

            if chrom_rank != region_rank {
                continue;
            }

            if !seen_any_chrom_overlap {
                seen_any_chrom_overlap = true;
            }

            if variant_overlaps_region(record, &region) {
                // Analyze each alt allele
                let allele_count = record.allele_count() as usize;
                for alt_idx in 0..(allele_count - 1) {
                    match analyze_variant(
                        record,
                        &header,
                        alt_idx,
                        &region,
                        samples_index_map,
                        is_phred,
                    )? {
                        Some(analysis) => {
                            let summary = region_summary.get_or_insert_with(|| RegionSummary {
                                variants: Vec::new(),
                            });
                            summary.add_variant(analysis);
                        }
                        None => {
                            // Filtered
                        }
                    }
                }
            }
        }
        /* ------------------------------------------------------------- */
        /* ------------------------------------------------------------- */

        if let Some(summary) = region_summary {
            results.push(summary);
        }
    }
    /* ======================================================== */
    /* ====== End Main Loop: Process each BED region ========== */

    /* ======================================================== */
    /* === Finalization: Validate and return results ========== */
    if !seen_any_chrom_overlap {
        return Err(Error::MsiVcfChromMismatch.into());
    }

    let stats = IntersectionStats {
        total_regions,
        skipped_invalid_motif: skipped_invalid_region,
    };
    stats.log_summary();
    info!("Number of intersected regions: {}", results.len());

    Ok((results, stats.valid_regions()))
}

/* ================================================ */

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use std::io::Write;

    use rust_htslib::bcf::{self, Read};
    use tempfile::NamedTempFile;

    use crate::utils::bcf_utils::tests::{
        create_multi_chromosome_vcf, create_test_vcf, TestVcfConfig,
    };
    use crate::utils::stats::TEST_EPSILON;

    /* ============ Bed Helper  ====================== */

    fn create_multi_region_bed() -> NamedTempFile {
        let tmp_bed = NamedTempFile::new().unwrap();
        writeln!(tmp_bed.as_file(), "chr1\t97\t104\t7xT").unwrap(); // pos 100 inside [97, 104)
        writeln!(tmp_bed.as_file(), "chr1\t197\t204\t7xT").unwrap(); // pos 200 inside [197, 204)
        writeln!(tmp_bed.as_file(), "chr2\t147\t154\t7xT").unwrap(); // pos 150 inside [147, 154)
        writeln!(tmp_bed.as_file(), "chrX\t172\t179\t7xT").unwrap(); // pos 175 inside [172, 179)
        tmp_bed
    }

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
        /* We are considering even low indels like in the last case. TODO: Check if discarding is better behaviour. */
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
        /* Last test: Special Case when tail remains in ref/alt apart from anchor.
            Here we consider it NA as it's not a clean indel.
        */
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

    /* ========== analyze_variant tests ============== */

    #[test]
    fn test_analyze_variant_filters_snv() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"A",
            alt_alleles: vec![b"T"],
            af_values: Some(vec![0.5, 0.5]),
            ..Default::default()
        });

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 200,
            motif: "A".to_string(),
        };

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);

        let result =
            analyze_variant(&record, &header, 0, &region, &samples_index_map, false).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_analyze_variant_perfect_indel() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"ACAG",
            alt_alleles: vec![b"ACAGCAG"],
            af_values: Some(vec![0.5, 0.8]),
            prob_absent: Some(vec![0.01]),
            prob_artifact: Some(vec![0.005]),
            ..Default::default()
        });

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 200,
            motif: "CAG".to_string(),
        };

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);
        samples_index_map.insert("sample2".to_string(), 1);

        let result =
            analyze_variant(&record, &header, 0, &region, &samples_index_map, false).unwrap();

        assert!(result.is_some());
        let variant = result.unwrap();
        assert!((variant.prob_absent - 0.015).abs() < TEST_EPSILON);
        assert_eq!(variant.sample_afs.len(), 2);
        assert!((variant.sample_afs["sample1"] - 0.5).abs() < TEST_EPSILON);
        assert!((variant.sample_afs["sample2"] - 0.8).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_analyze_variant_phred_probabilities() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"ACAG",
            alt_alleles: vec![b"ACAGCAG"],
            af_values: Some(vec![0.5, 0.8]),
            use_phred: true,
            ..Default::default()
        });

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 200,
            motif: "CAG".to_string(),
        };

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);

        let result =
            analyze_variant(&record, &header, 0, &region, &samples_index_map, true).unwrap();

        assert!(result.is_some());
        let variant = result.unwrap();
        assert!((variant.prob_absent - 0.11).abs() < 0.01);
    }

    #[test]
    fn test_analyze_variant_multi_allelic() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"A",
            alt_alleles: vec![b"T", b"ATG"],
            af_values: Some(vec![0.3, 0.3, 0.6, 0.6]),
            ..Default::default()
        });

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 200,
            motif: "TG".to_string(),
        };

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);

        assert!(
            analyze_variant(&record, &header, 0, &region, &samples_index_map, false)
                .unwrap()
                .is_none()
        );
        assert!(
            analyze_variant(&record, &header, 1, &region, &samples_index_map, false)
                .unwrap()
                .is_some()
        );
    }

    #[test]
    fn test_analyze_variant_filters_symbolic() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"A",
            alt_alleles: vec![b"<DEL>"],
            ..Default::default()
        });

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 200,
            motif: "A".to_string(),
        };

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);

        let result =
            analyze_variant(&record, &header, 0, &region, &samples_index_map, false).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_analyze_variant_downgrades_without_prob() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();

        let mut header = rust_htslib::bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");
        header.push_record(br##"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">"##);
        // Headers defined but won't add values to the record for PROB_ABSENT and PROB_ARTIFACT
        header.push_record(
            br##"##INFO=<ID=PROB_ABSENT,Number=A,Type=Float,Description="Probability absent">"##,
        );
        header.push_record(br##"##INFO=<ID=PROB_ARTIFACT,Number=A,Type=Float,Description="Probability artifact">"##);
        header.push_sample(b"sample1");
        header.push_record(
            br##"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">"##,
        );

        let mut wtr = rust_htslib::bcf::Writer::from_path(
            path,
            &header,
            false,
            rust_htslib::bcf::Format::Vcf,
        )
        .unwrap();

        let mut rec = wtr.empty_record();
        rec.set_rid(Some(0));
        rec.set_pos(99);
        rec.set_alleles(&[b"ACAG", b"ACAGCAG"]).unwrap();
        rec.push_info_integer(b"SVLEN", &[3]).unwrap();
        rec.push_format_float(b"AF", &[0.5]).unwrap();
        wtr.write(&rec).unwrap();
        drop(wtr);

        let mut reader = bcf::Reader::from_path(path).unwrap();
        let header_view = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 200,
            motif: "CAG".to_string(),
        };

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);

        let result =
            analyze_variant(&record, &header_view, 0, &region, &samples_index_map, false).unwrap();
        assert!(result.is_none());
    }

    /* ====== variant_overlaps_region tests ========== */

    #[test]
    fn test_variant_overlaps_region() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig::default());
        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let record = reader.records().next().unwrap().unwrap();

        // Inside
        assert!(variant_overlaps_region(
            &record,
            &BedRegion {
                chrom: "chr1".to_string(),
                start: 0,
                end: 200,
                motif: "A".to_string(),
            }
        ));

        // At start (inclusive)
        assert!(variant_overlaps_region(
            &record,
            &BedRegion {
                chrom: "chr1".to_string(),
                start: 99,
                end: 200,
                motif: "A".to_string(),
            }
        ));

        // Before region
        assert!(!variant_overlaps_region(
            &record,
            &BedRegion {
                chrom: "chr1".to_string(),
                start: 100,
                end: 200,
                motif: "A".to_string(),
            }
        ));

        // At end (exclusive)
        assert!(!variant_overlaps_region(
            &record,
            &BedRegion {
                chrom: "chr1".to_string(),
                start: 0,
                end: 99,
                motif: "A".to_string(),
            }
        ));
    }

    /* ========== intersect_streaming tests ========== */

    #[test]
    fn test_intersect_streaming_basic() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"ACAG",
            alt_alleles: vec![b"ACAGCAG"],
            af_values: Some(vec![0.5, 0.8]),
            prob_absent: Some(vec![0.01]),
            prob_artifact: Some(vec![0.005]),
            ..Default::default()
        });

        let tmp_bed = NamedTempFile::new().unwrap();
        writeln!(tmp_bed.as_file(), "chr1\t96\t105\t3xCAG").unwrap();

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);
        samples_index_map.insert("sample2".to_string(), 1);

        let (results, total_regions) = intersect_streaming(
            &tmp_bed.path().to_path_buf(),
            &tmp_vcf.path().to_path_buf(),
            &samples_index_map,
        )
        .unwrap();

        assert_eq!(total_regions, 1);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].variants.len(), 1);
    }

    #[test]
    fn test_intersect_streaming_no_overlap() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            af_values: Some(vec![0.5, 0.8]),
            ..Default::default()
        });

        let tmp_bed = NamedTempFile::new().unwrap();
        writeln!(tmp_bed.as_file(), "chr1\t200\t209\t3xCAG").unwrap();

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);

        let (results, total_regions) = intersect_streaming(
            &tmp_bed.path().to_path_buf(),
            &tmp_vcf.path().to_path_buf(),
            &samples_index_map,
        )
        .unwrap();

        assert_eq!(total_regions, 1);
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_intersect_streaming_multi_chromosome() {
        let tmp_vcf = create_multi_chromosome_vcf();
        let tmp_bed = create_multi_region_bed();

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);

        let (results, total_regions) = intersect_streaming(
            &tmp_bed.path().to_path_buf(),
            &tmp_vcf.path().to_path_buf(),
            &samples_index_map,
        )
        .unwrap();

        assert_eq!(total_regions, 4);
        assert!(results.len() == 4);
        assert!(results[0].variants.len() == 1);
        assert!(results[1].variants.len() == 1);
        assert!(results[2].variants.len() == 1);
        assert!(results[3].variants.len() == 1);
    }
}
