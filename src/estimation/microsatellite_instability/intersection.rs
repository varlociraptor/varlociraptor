//! intersection.rs
//!
//! Streaming intersection and variant analysis
//! Processes BED regions and VCF variants to identify perfect microsatellite repeats.
//! Prepares variants for DP analysis by extracting probability data and filtering not
//! required variants(probability missing, no af value in all samples, not a perfect repeat).

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
#[derive(Debug, Clone, PartialEq)]
enum RepeatStatus {
    Perfect,
    NA,
}

#[derive(Debug, Clone)]
pub(super) struct Variant {
    pub(super) prob_absent: f64,
    pub(super) sample_afs: HashMap<String, f64>,
}

#[derive(Debug)]
pub(super) struct RegionSummary {
    pub(super) variants: Vec<Variant>,
}

impl RegionSummary {
    fn add_variant(&mut self, variant: Variant) {
        self.variants.push(variant);
    }
}

#[derive(Debug, Clone)]
pub(super) struct IntersectionStats {
    pub total_regions: usize,
    pub skipped_invalid_motif: usize,
}

impl IntersectionStats {
    fn valid_regions(&self) -> usize {
        self.total_regions - self.skipped_invalid_motif
    }

    fn log_summary(&self) {
        info!("  Intersection complete:");
        info!("    Total regions in BED: {}", self.total_regions);
        info!(
            "    Valid regions (1-6 motif length): {}",
            self.valid_regions()
        );
        info!(
            "    Skipped (invalid motif): {}",
            self.skipped_invalid_motif
        );
    }
}
/* ================================================ */

/* ======== Variant Analysis Functions ============ */
/// Check if a sequence is a perfect repeat of a motif
///
/// This function determines if an indel represents a perfect tandem repeat
/// of the microsatellite motif by:
/// 1. Finding the anchor (full matching prefix between REF and ALT)
/// 2. Extracting the changed sequence (non-matching portion after anchor)
/// 3. Checking if changed sequence is a perfect repeat of the motif
///
/// Example:
///   REF:     ACAG        (anchor: ACAG, changed: none)
///   ALT:     ACAGCAG     (anchor: ACAG, changed: CAG)
///   Motif:   CAG
///   Anchor:  ACAG (all 4 bases match)
///   Changed: CAG (only the non-matching part)
///   Result:  CAG = 1× CAG -> Perfect repeat!
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

    // 2. Extracting the changed sequence
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

    // 3. Check if changed sequence is a perfect repeat of the motif
    let motif_bytes: Vec<u8> = motif.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let motif_len = motif_bytes.len();

    /* NOTE: In case we toggle the Error on motif.len() in bed parsing off, turn this onn.*/
    // if motif_len == 0 {
    //     return RepeatStatus::NA
    // }

    if changed_seq.len() % motif_len != 0 {
        // "Changed sequence length ({}) is not a multiple of motif length ({})"
        return RepeatStatus::NA;
    }

    for (i, &base) in changed_seq.iter().enumerate() {
        let expected_base = motif_bytes[i % motif_len];
        if base.to_ascii_uppercase() != expected_base {
            // Not a perfect repeat: position {} has {} but expected {}
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
/// 1. **Filtering**: Skip variants that aren't relevant for MSI:
///    - SNVs (single nucleotide variants)
///    - Symbolic alleles (<DEL>, <INS>)
///    - Breakends (complex structural variants)
///    - Spanning deletions (*)
///    - Non-indel variants
/// 2. **SVLEN Calculation**: Determine the length of insertion/deletion///
/// 3. **Perfect Repeat Check**: Determine if indel is a perfect repeat
/// 4. **Probability Extraction**: Get PROB_ABSENT from INFO field
/// 5. **Allele Frequency Extraction**: Get per-sample AFs from FORMAT:AF
///
/// # Returns
/// * `Ok(None)` - Variant should be skipped (SNV, symbolic, etc.)
/// * `Ok(Some(VariantAnalysis))` - Valid indel with extracted data
/// * `Err(...)` - Error occurred (invalid data, I/O error, etc.)
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
            // reason: "Perfect repeat but missing required probability (PROB_ABSENT) or All sample Afs"
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
#[inline]
fn variant_overlaps_region(record: &bcf::Record, region: &BedRegion) -> bool {
    // Variant position (0-based)
    let pos = record.pos() as u64;
    // BED region (0-based, half-open [start, end))
    pos >= region.start && pos < region.end
}

/// Perform streaming intersection of BED regions with VCF variants
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
        "  Probabilities in VCF/BCF are {} scaled",
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
        let bed_record = bed_result.map_err(|e| Error::BedRecordRead {
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
                    return Err(Error::VcfRecordRead {
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
        return Err(Error::NoChromosomeMatch.into());
    }

    let stats = IntersectionStats {
        total_regions,
        skipped_invalid_motif: skipped_invalid_region,
    };
    stats.log_summary();
    info!(" Number of intersected regions: {}", results.len());

    Ok((results, stats.valid_regions()))
}

/* ================================================ */

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_is_perfect_repeat_simple_insertion() {
        let ref_seq = b"ACAG";
        let alt_seq = b"ACAGCAG";
        let motif = "CAG";
        let svlen = 3;

        let status = is_perfect_repeat(alt_seq, svlen, motif, ref_seq);
        assert_eq!(status, RepeatStatus::Perfect);
    }

    #[test]
    fn test_is_perfect_repeat_simple_deletion() {
        let ref_seq = b"ACAGCAG";
        let alt_seq = b"ACAG";
        let motif = "CAG";
        let svlen = -3;

        let status = is_perfect_repeat(alt_seq, svlen, motif, ref_seq);
        assert_eq!(status, RepeatStatus::Perfect);
    }

    #[test]
    fn test_is_perfect_repeat_multiple_units() {
        let ref_seq = b"ACAG";
        let alt_seq = b"ACAGCAGCAG"; // 2 CAG units inserted
        let motif = "CAG";
        let svlen = 6;

        let status = is_perfect_repeat(alt_seq, svlen, motif, ref_seq);
        assert_eq!(status, RepeatStatus::Perfect);
    }

    #[test]
    fn test_is_perfect_repeat_not_perfect() {
        let ref_seq = b"ACAG";
        let alt_seq = b"ACAGCAT"; // CAT, not CAG
        let motif = "CAG";
        let svlen = 3;

        let status = is_perfect_repeat(alt_seq, svlen, motif, ref_seq);
        assert_eq!(status, RepeatStatus::NA);
    }

    #[test]
    fn test_is_perfect_repeat_case_insensitive() {
        let ref_seq = b"acag";
        let alt_seq = b"acagCAG"; // Mixed case
        let motif = "cag";
        let svlen = 3;

        let status = is_perfect_repeat(alt_seq, svlen, motif, ref_seq);
        assert_eq!(status, RepeatStatus::Perfect);
    }

    #[test]
    fn test_is_perfect_repeat_partial_motif() {
        let ref_seq = b"ACAG";
        let alt_seq = b"ACAGCA"; // Only "CA" added, not full "CAG"
        let motif = "CAG";
        let svlen = 2;

        let status = is_perfect_repeat(alt_seq, svlen, motif, ref_seq);
        assert_eq!(status, RepeatStatus::NA);
    }

    #[test]
    fn test_is_perfect_repeat_dinucleotide() {
        let ref_seq = b"ACA";
        let alt_seq = b"ACACA"; // "CA" dinucleotide
        let motif = "CA";
        let svlen = 2;

        let status = is_perfect_repeat(alt_seq, svlen, motif, ref_seq);
        assert_eq!(status, RepeatStatus::Perfect);
    }

    #[test]
    fn test_is_perfect_repeat_mononucleotide() {
        let ref_seq = b"AT";
        let alt_seq = b"ATT"; // "T" mononucleotide
        let motif = "T";
        let svlen = 1;

        let status = is_perfect_repeat(alt_seq, svlen, motif, ref_seq);
        assert_eq!(status, RepeatStatus::Perfect);
    }

    #[test]
    fn test_is_perfect_repeat_svlen_zero() {
        let ref_seq = b"ACAG";
        let alt_seq = b"ACAT"; // SNV, not indel
        let motif = "CAG";
        let svlen = 0;

        let status = is_perfect_repeat(alt_seq, svlen, motif, ref_seq);
        assert_eq!(status, RepeatStatus::NA);
    }

    #[test]
    fn test_is_perfect_repeat_anchor_zero() {
        let ref_seq = b"A";
        let alt_seq = b"TCAG";
        let motif = "TCAG";
        let svlen = 4; // No Anchor match

        let status = is_perfect_repeat(alt_seq, svlen, motif, ref_seq);
        assert_eq!(status, RepeatStatus::NA);
    }

    #[test]
    fn test_is_perfect_repeat_length_mismatch() {
        let ref_seq = b"ACAG";
        let alt_seq = b"ACAGCAG";
        let motif = "CAG";
        let svlen = 5; // Wrong! Should be 3

        let status = is_perfect_repeat(alt_seq, svlen, motif, ref_seq);
        assert_eq!(status, RepeatStatus::NA);
    }

    /// Helper to encode genotype alleles for VCF FORMAT:GT
    fn encode_genotype_allele(allele_index: i32, phased: bool) -> i32 {
        let phased_flag = if phased { 1 } else { 0 };
        (allele_index + 1) * 2 | phased_flag
    }

    /// Create a test VCF file with configurable alleles
    /// Returns (temp_file, sample_names)
    fn create_test_vcf_with_alleles(
        ref_allele: &[u8],
        alt_alleles: Vec<&[u8]>,
        af_values: Option<Vec<f32>>,
        use_phred: bool,
    ) -> (NamedTempFile, Vec<String>) {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();

        let mut header = rust_htslib::bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");
        header.push_record(br##"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">"##);

        if use_phred {
            header.push_record(br##"##INFO=<ID=PROB_ABSENT,Number=A,Type=Float,Description="Probability absent (PHRED)">"##);
            header.push_record(br##"##INFO=<ID=PROB_ARTIFACT,Number=A,Type=Float,Description="Probability artifact (PHRED)">"##);
        } else {
            header.push_record(br##"##INFO=<ID=PROB_ABSENT,Number=A,Type=Float,Description="Probability absent (linear)">"##);
            header.push_record(br##"##INFO=<ID=PROB_ARTIFACT,Number=A,Type=Float,Description="Probability artifact (linear)">"##);
        }

        header.push_sample(b"sample1");
        header.push_sample(b"sample2");
        header.push_record(br##"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"##);
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

        // Build alleles array: [REF, ALT1, ALT2, ...]
        let mut alleles_vec = vec![ref_allele];
        alleles_vec.extend(alt_alleles.iter());
        rec.set_alleles(&alleles_vec).unwrap();

        // Add SVLEN for each ALT
        let svlen_values: Vec<i32> = alt_alleles
            .iter()
            .map(|alt| (alt.len() as i32) - (ref_allele.len() as i32))
            .collect();
        rec.push_info_integer(b"SVLEN", &svlen_values).unwrap();

        let (prob_absent_values, prob_artifact_values) = if use_phred {
            // For Phred: use values that convert to reasonable probabilities
            // Phred 20 = 0.01 probability, Phred 10 = 0.1 probability
            (vec![20.0; alt_alleles.len()], vec![10.0; alt_alleles.len()])
        } else {
            // For linear: use direct probabilities that sum appropriately
            (
                vec![0.005; alt_alleles.len()],
                vec![0.005; alt_alleles.len()],
            )
        };

        // // Add PROB_ABSENT for each ALT
        // let prob_absent_values: Vec<f32> = vec![0.01; alt_alleles.len()];
        rec.push_info_float(b"PROB_ABSENT", &prob_absent_values)
            .unwrap();

        // // Add PROB_ARTIFACT for each ALT
        // let prob_artifact_values: Vec<f32> = vec![0.01; alt_alleles.len()];
        rec.push_info_float(b"PROB_ARTIFACT", &prob_artifact_values)
            .unwrap();

        // Add GT (0/1 for sample1, 1/1 for sample2)
        let genotypes_data = [
            encode_genotype_allele(0, false),
            encode_genotype_allele(1, false),
            encode_genotype_allele(1, false),
            encode_genotype_allele(1, false),
        ];
        rec.push_format_integer(b"GT", &genotypes_data).unwrap();

        // Add AF values if provided
        if let Some(afs) = af_values {
            rec.push_format_float(b"AF", &afs).unwrap();
        }

        wtr.write(&rec).unwrap();

        (tmp, vec!["sample1".to_string(), "sample2".to_string()])
    }

    #[test]
    fn test_analyze_variant_filters_snv() {
        let (tmp_vcf, _) = create_test_vcf_with_alleles(
            b"A",
            vec![b"T"], // SNV
            // Some(vec![0.5]),
            Some(vec![0.5, 0.5]),
            false,
        );

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
    fn test_analyze_variant_analyzes_indel() {
        let (tmp_vcf, _) = create_test_vcf_with_alleles(
            b"ACAG",
            vec![b"ACAGCAG"], // Insertion
            Some(vec![0.5, 0.8]),
            false,
        );

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
        assert!((variant.prob_absent - 0.01).abs() < 1e-6);
        assert_eq!(variant.sample_afs.len(), 2);
        assert!((variant.sample_afs["sample1"] - 0.5).abs() < 1e-6);
        assert!((variant.sample_afs["sample2"] - 0.8).abs() < 1e-6);
    }

    #[test]
    fn test_analyze_variant_multi_allelic_site() {
        let (tmp_vcf, _) = create_test_vcf_with_alleles(
            b"A",
            vec![b"T", b"ATG"],             // SNV + Indel
            Some(vec![0.3, 0.3, 0.6, 0.6]), // 2 samples × 2 alts = 4 values
            false,
        );

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

        // alt_idx=0 (T, SNV) should be filtered
        let result_snv =
            analyze_variant(&record, &header, 0, &region, &samples_index_map, false).unwrap();
        assert!(result_snv.is_none());

        // alt_idx=1 (ATG, indel) should be analyzed
        let result_indel =
            analyze_variant(&record, &header, 1, &region, &samples_index_map, false).unwrap();
        assert!(result_indel.is_some());
    }

    #[test]
    fn test_analyze_variant_filters_symbolic() {
        let (tmp_vcf, _) = create_test_vcf_with_alleles(b"A", vec![b"<DEL>"], None, false);

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

    // #[test]
    // fn test_analyze_variant_downgrades_perfect_without_prob() {
    //     // Create VCF without PROB_ABSENT field
    //     let tmp = NamedTempFile::new().unwrap();
    //     let path = tmp.path();

    //     let mut header = rust_htslib::bcf::Header::new();
    //     header.push_record(br"##fileformat=VCFv4.2");
    //     header.push_record(br"##contig=<ID=chr1,length=1000000>");
    //     header.push_record(br##"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">"##);
    //     // NO PROB_ABSENT field!
    //     header.push_sample(b"sample1");
    //     header.push_record(br##"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"##);
    //     header.push_record(br##"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">"##);
    //     header.push_record(br##"##INFO=<ID=PROB_ABSENT,Number=A,Type=Float,Description="Probability absent">"##);

    //     let mut wtr = rust_htslib::bcf::Writer::from_path(
    //         path, &header, false, rust_htslib::bcf::Format::Vcf
    //     ).unwrap();

    //     let mut rec = wtr.empty_record();
    //     rec.set_rid(Some(0));
    //     rec.set_pos(99);
    //     rec.set_alleles(&[b"ACAG", b"ACAGCAG"]).unwrap();
    //     rec.push_info_integer(b"SVLEN", &[3]).unwrap();
    //     rec.push_format_float(b"AF", &[0.5]).unwrap();
    //     wtr.write(&rec).unwrap();

    //     let mut reader = bcf::Reader::from_path(path).unwrap();
    //     let header_view = reader.header().clone();
    //     let record = reader.records().next().unwrap().unwrap();

    //     let region = BedRegion {
    //         chrom: "chr1".to_string(),
    //         start: 0,
    //         end: 200,
    //         motif: "CAG".to_string(),
    //     };

    //     let mut samples_index_map = HashMap::new();
    //     samples_index_map.insert("sample1".to_string(), 0);

    //     let result = analyze_variant(&record, &header_view, 0, &region, &samples_index_map, false).unwrap();

    //     // Should still return Some, but downgraded to NA
    //     assert!(result.is_some());
    //     let analysis = result.unwrap();
    //     assert!(matches!(analysis.repeat_status, RepeatStatus::NA));
    //     assert!(analysis.prob_absent.is_none());
    // }

    #[test]
    fn test_variant_overlaps_region_inside() {
        let (tmp_vcf, _) = create_test_vcf_with_alleles(b"A", vec![b"T"], None, false);
        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let record = reader.records().next().unwrap().unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 200,
            motif: "A".to_string(),
        };

        // Variant at position 99, region [0, 200)
        assert!(variant_overlaps_region(&record, &region));
    }

    #[test]
    fn test_variant_overlaps_region_before() {
        let (tmp_vcf, _) = create_test_vcf_with_alleles(b"A", vec![b"T"], None, false);
        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let record = reader.records().next().unwrap().unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 100, // Variant at 99, region starts at 100
            end: 200,
            motif: "A".to_string(),
        };

        // Variant at 99 < start 100
        assert!(!variant_overlaps_region(&record, &region));
    }

    #[test]
    fn test_variant_overlaps_region_at_start() {
        let (tmp_vcf, _) = create_test_vcf_with_alleles(b"A", vec![b"T"], None, false);
        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let record = reader.records().next().unwrap().unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 99, // Variant at 99, region starts at 99
            end: 200,
            motif: "A".to_string(),
        };

        // Variant at 99 >= start 99
        assert!(variant_overlaps_region(&record, &region));
    }

    #[test]
    fn test_variant_overlaps_region_at_end() {
        let (tmp_vcf, _) = create_test_vcf_with_alleles(b"A", vec![b"T"], None, false);
        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let record = reader.records().next().unwrap().unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 99, // Variant at 99, region ends at 99 (exclusive)
            motif: "A".to_string(),
        };

        // Variant at 99 >= end 99 (half-open interval)
        assert!(!variant_overlaps_region(&record, &region));
    }

    #[test]
    fn test_variant_overlaps_region_after() {
        let (tmp_vcf, _) = create_test_vcf_with_alleles(b"A", vec![b"T"], None, false);
        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let record = reader.records().next().unwrap().unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 50, // Variant at 99, region ends at 50
            motif: "A".to_string(),
        };

        // Variant at 99 >= end 50
        assert!(!variant_overlaps_region(&record, &region));
    }

    #[test]
    fn test_intersect_streaming_basic() {
        // Create test VCF
        let (tmp_vcf, _) = create_test_vcf_with_alleles(
            b"ACAG",
            vec![b"ACAGCAG"], // Perfect repeat insertion
            Some(vec![0.5, 0.8]),
            false,
        );

        // Create test BED
        let tmp_bed = NamedTempFile::new().unwrap();
        writeln!(tmp_bed.as_file(), "chr1\t0\t200\t3xCAG").unwrap();

        // Run intersection
        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);
        samples_index_map.insert("sample2".to_string(), 1);

        let (results, total_regions) = intersect_streaming(
            &tmp_bed.path().to_path_buf(),
            &tmp_vcf.path().to_path_buf(),
            &samples_index_map,
        )
        .unwrap();

        // Verify results
        assert_eq!(total_regions, 1);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].variants.len(), 1);
        assert!(
            results[0].variants[0].prob_absent >= 0.0 && results[0].variants[0].prob_absent <= 1.0
        ); //TODO:
    }

    #[test]
    fn test_intersect_streaming_no_overlap() {
        // Create test VCF at position 99
        let (tmp_vcf, _) =
            create_test_vcf_with_alleles(b"A", vec![b"AT"], Some(vec![0.5, 0.8]), false);

        // Create test BED far from variant
        let tmp_bed = NamedTempFile::new().unwrap();
        writeln!(tmp_bed.as_file(), "chr1\t200\t300\t3xCAG").unwrap();

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);

        let (results, total_regions) = intersect_streaming(
            &tmp_bed.path().to_path_buf(),
            &tmp_vcf.path().to_path_buf(),
            &samples_index_map,
        )
        .unwrap();

        // Verify no overlap
        assert_eq!(total_regions, 1);
        assert_eq!(results.len(), 0);
    }
}
