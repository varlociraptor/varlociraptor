//! extraction.rs
//!
//! VCF streaming and variant extraction for MSI calling.
//!
//! # This module provides:
//! 1. Data Structures (`Variant`, `RegionSummary`, `ExtractionStats`)
//! 2. Main Extraction (`extract_regions`)
//!
//! # Input Contract
//! Input VCF must be the output of `varlociraptor preprocess msi` followed by
//! `varlociraptor call variants`. This guarantees:
//! - Records with `INFO/REGION_ID` are MS indels (real or dummy)
//! - `INFO/MSI_DUMMY` marks dummy records (no real indel observed in the region)
//! - `INFO/REGION_ID` is `Number=A` - one value per ALT allele
//! - `INFO/PROB_{EVENT}` fields are expected for all user-specified events
//!   (missing values are skipped and counted in `ExtractionStats`)
//! - `FORMAT/AF` is expected for the sample of interest
//!   (missing values are skipped and counted in `ExtractionStats`)
//!
//! # Design
//! Records without `INFO/REGION_ID` are skipped - non-MS variants that
//! passed through preprocessing unchanged. Per-ALT REGION_ID values of
//! `"."` indicate unmatched ALTs in multi-allelic records and are skipped.
//! Regions are grouped via `HashMap` for O(1) lookup and sorted by
//! `(chrom, start)` before return.
//!

use std::collections::HashMap;

use anyhow::Result;
use log::{debug, info, warn};
use rust_htslib::bcf::{self, Read};

use crate::constants::{MSI_DUMMY_TAG, MSI_REGION_ID_TAG};
use crate::errors::Error;
use crate::utils::bcf_utils::{
    get_chrom, get_events_probability, get_info_strings, get_sample_af, record_has_info_flag,
    record_has_info_string,
};

/* ============ Data Structures =================== */

/// A single variant contributing to MSI analysis for one region.
///
/// Note: In multi-allelic records each ALT allele yields a separate Variant.
#[derive(Debug)]
pub(super) struct Variant {
    /// P(variant absent) = 1.0 - P(at least one specified event).
    pub prob_absent: f64,
    /// Allele frequency for the sample (FORMAT:AF).
    pub af: f32,
}

/// All variants observed within one microsatellite region.
#[derive(Debug)]
pub(super) struct RegionSummary {
    /// Chromosome — used for window assignment in heatmap analysis.
    pub chrom: String,
    /// Region start (0-based) — used for window assignment in heatmap analysis.
    pub start: u64,
    /// Variants in this region (real or dummy, ordered by encounter).
    pub variants: Vec<Variant>,
    /// True if at least one non-dummy indel record was encountered,
    /// regardless of whether prob or AF extraction succeeded.
    /// Used to report the "real indel only" MSI score alongside
    /// the full (real + dummy) score.
    pub has_real_indel: bool,
}

impl RegionSummary {
    /// Parse a REGION_ID string ("chrom:start-end") into a RegionSummary.
    ///
    /// Returns Err if the format is invalid so the caller can warn and skip.
    ///
    /// # Errors
    /// Returns `Error::MsiRegionIdMalformed` if:
    /// - No `:` separator between chrom and coordinates
    /// - No `-` separator between start and end
    /// - Start coordinate is not a valid integer
    fn from_region_id(region_id: &str) -> Result<Self> {
        let (chrom, coords) =
            region_id
                .split_once(':')
                .ok_or_else(|| Error::MsiRegionIdMalformed {
                    region_id: region_id.to_string(),
                    details: "missing ':' separator".to_string(),
                })?;

        let (start_str, _end_str) =
            coords
                .split_once('-')
                .ok_or_else(|| Error::MsiRegionIdMalformed {
                    region_id: region_id.to_string(),
                    details: "missing '-' in coordinate part".to_string(),
                })?;

        let start: u64 = start_str.parse().map_err(|_| Error::MsiRegionIdMalformed {
            region_id: region_id.to_string(),
            details: format!("start '{}' is not a valid integer", start_str),
        })?;

        Ok(RegionSummary {
            chrom: chrom.to_string(),
            start,
            variants: Vec::new(),
            has_real_indel: false,
        })
    }
}

/// Statistics collected during extraction.
#[derive(Debug, Default)]
pub(super) struct ExtractionStats {
    /// Total VCF records read (MS and non-MS: one per line, regardless of ALT count).
    pub total_records: usize,
    /// Records skipped — no REGION_ID, not MS-relevant.
    /// Includes SNVs and non-perfect indels that passed through preprocessing unchanged.
    /// Note: Counted per record, not per allele.
    pub skipped_non_ms: usize,
    /// Unique MS regions encountered — used as the MSI score denominator.
    pub total_ms_regions: usize,
    /// Regions with at least one non-dummy indel record, counted the moment
    /// the region is first marked real - independent of whether it later
    /// retains any variants after AF/prob filtering or heatmap-gated pruning.
    pub regions_with_real_indel: usize,
    /// Dummy indel records processed (MSI_DUMMY flag set), counted per record.
    /// Each corresponds to a region where no perfect indel was observed in reads.
    /// Note: Each individual dummy record has one ALT by construction; since
    /// inject_dummy_indels writes a deletion+insertion pair per region, a
    /// region with no perfect indel contributes 2 to this count, not 1.
    pub dummy_records: usize,
    /// ALT Alleles skipped — FORMAT:AF absent for the sample.
    pub skipped_missing_af: usize,
    /// ALT Alleles skipped — event probability field absent.
    pub skipped_missing_prob: usize,
}

impl ExtractionStats {
    /// Log extraction statistics alongside region-level counts.
    pub fn log_stats(&self) {
        info!("Extraction statistics:");
        info!("  - Total records read:             {}", self.total_records);
        info!(
            "  - Non-MS records skipped:         {}",
            self.skipped_non_ms
        );
        info!(
            "  - Total MS regions:               {}",
            self.total_ms_regions
        );
        info!(
            "  - Regions with real indel:        {}",
            self.regions_with_real_indel
        );
        info!(
            "  - Regions needing a dummy indel:  {}",
            self.dummy_records / 2
        );
        info!(
            "  - ALT alleles skipped (no AF):    {}",
            self.skipped_missing_af
        );
        info!(
            "  - ALT alleles skipped (no prob):  {}",
            self.skipped_missing_prob
        );
    }
}

/// Stream a preprocessed+called VCF and extract per-region variant summaries.
///
/// Processes only records with `INFO/REGION_ID`, grouping variants by region
/// via a HashMap index into a Vec.
///
/// `has_real_indel` is set for regions where at least one non-dummy record
/// was encountered, regardless of whether prob or AF extraction succeeded.
/// This allows reporting a real-indel-only MSI score alongside the full score.
///
/// # Arguments
/// * `vcf`        - Reader for the called VCF
/// * `sample`     - Sample name (used in log messages)
/// * `sample_idx` - Index of this sample in FORMAT columns
/// * `events`     - Event names to combine (e.g., `["somatic"]`)
/// * `is_phred`   - Whether `INFO/PROB_*` values are PHRED-scaled
/// * `needs_heatmap` - Whether heatmap output was requested; if `false`,
///   regions with zero usable variants are pruned before returning, since
///   heatmap's per-window denominator requires every region present, and
///   window generation itself depends on every region's position - a
///   pruned variant-less region could otherwise cause its whole window
///   (or chromosome) to never appear on the heatmap at all.
///   (distribution/pseudotime are unaffected either way - see
///   filter_regions_by_af in dp_analysis.rs).
///
/// # Returns
/// `(Vec<RegionSummary>, ExtractionStats)`. `stats.total_ms_regions` remains
/// the MSI score denominator regardless of pruning. Vec length equals
/// `stats.total_ms_regions` only when `needs_heatmap` is `true`; otherwise
/// it may be smaller, since variant-less regions are pruned.
///
/// # Example
/// assert!(extract_regions(&mut vcf, "sample1", 0, &["somatic".to_string()], true).is_ok());
pub(super) fn extract_regions(
    vcf: &mut bcf::Reader,
    sample: &str,
    sample_idx: usize,
    events: &[String],
    is_phred: bool,
    needs_heatmap: bool,
) -> Result<(Vec<RegionSummary>, ExtractionStats)> {
    let mut regions: Vec<RegionSummary> = Vec::new();
    let mut region_index: HashMap<String, usize> = HashMap::new();
    let mut stats = ExtractionStats::default();
    let mut record = vcf.empty_record();

    loop {
        match vcf.read(&mut record) {
            None => break,
            Some(Err(e)) => {
                return Err(Error::VcfRecordReadFailed {
                    details: e.to_string(),
                }
                .into())
            }
            Some(Ok(())) => {
                stats.total_records += 1;
            }
        }

        // Preprocessing guarantees REGION_ID is only written when at least one
        // ALT matched a BED region - all-dot records are never written.
        // So presence of the field means at least one non-dot value exists.
        if !record_has_info_string(&record, MSI_REGION_ID_TAG) {
            stats.skipped_non_ms += 1;
            continue;
        }

        let is_dummy = record_has_info_flag(&record, MSI_DUMMY_TAG);
        if is_dummy {
            stats.dummy_records += 1;
        }

        // Extract one Variant per ALT allele.
        let allele_count = record.allele_count() as usize;
        let region_ids = get_info_strings(&record, MSI_REGION_ID_TAG).unwrap_or_default();
        for alt_idx in 0..allele_count.saturating_sub(1) {
            let region_id = match region_ids.get(alt_idx).filter(|s| *s != ".").cloned() {
                Some(id) => id,
                None => {
                    debug!(
                        "ALT {} at {}:{} has no region - skipping",
                        alt_idx,
                        get_chrom(&record).unwrap_or_default(),
                        record.pos() + 1
                    );
                    continue;
                }
            };

            // Register region on first encounter, look up existing index otherwise.
            let region_idx = match region_index.get(&region_id) {
                Some(&idx) => idx,
                None => match RegionSummary::from_region_id(&region_id) {
                    Ok(summary) => {
                        let idx = regions.len();
                        regions.push(summary);
                        region_index.insert(region_id, idx);
                        stats.total_ms_regions += 1;
                        idx
                    }
                    Err(e) => {
                        warn!("Skipping record with malformed REGION_ID: {}", e);
                        continue;
                    }
                },
            };

            // Mark region as having a real (non-dummy) indel for this alt if applicable.
            // This is set regardless of whether prob/AF extraction succeeds -
            // the indel structurally existed in the reads.
            if !is_dummy && !regions[region_idx].has_real_indel {
                regions[region_idx].has_real_indel = true;
                stats.regions_with_real_indel += 1;
            }

            // Combine event probabilities, P(at least one event)
            let prob_events = match get_events_probability(&record, alt_idx, events, is_phred)? {
                Some(p) => p,
                None => {
                    stats.skipped_missing_prob += 1;
                    debug!(
                        "Event prob missing at {}:{} alt={} — skipping allele",
                        get_chrom(&record).unwrap_or_default(),
                        record.pos() + 1,
                        alt_idx
                    );
                    continue;
                }
            };

            // FORMAT:AF for this sample and ALT allele
            let af = match get_sample_af(&record, sample_idx, alt_idx)? {
                Some(a) => a,
                None => {
                    stats.skipped_missing_af += 1;
                    debug!(
                        "FORMAT:AF missing for '{}' at {}:{} alt={} — skipping allele",
                        sample,
                        get_chrom(&record).unwrap_or_default(),
                        record.pos() + 1,
                        alt_idx
                    );
                    continue;
                }
            };

            // prob_absent = 1 - P(events).
            regions[region_idx].variants.push(Variant {
                prob_absent: 1.0 - prob_events,
                af,
            });
        }
    }

    regions.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));

    if !needs_heatmap {
        regions.retain(|r| !r.variants.is_empty());
    }

    Ok((regions, stats))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bcf::{self, record::Numeric};
    use tempfile::NamedTempFile;

    use crate::constants::test_constants::{TEST_EPSILON_LOOSE, TEST_EPSILON_LOOSE_F32};
    use crate::constants::{MSI_DUMMY_HEADER, MSI_REGION_ID_HEADER};

    /* ====== Shared header builder =============== */

    /// Build a minimal BCF header for extraction tests.
    /// Includes REGION_ID (Number=A), MSI_DUMMY, PROB_SOMATIC, FORMAT:AF,
    /// and a single sample "tumor" on contigs chr1 and chr2.
    fn make_header() -> bcf::Header {
        let mut header = bcf::Header::new();
        header.push_record(br#"##fileformat=VCFv4.2"#);
        header.push_record(br#"##contig=<ID=chr1,length=10000000>"#);
        header.push_record(br#"##contig=<ID=chr2,length=10000000>"#);
        header.push_record(MSI_REGION_ID_HEADER.as_bytes());
        header.push_record(MSI_DUMMY_HEADER.as_bytes());
        header.push_record(
            br##"##INFO=<ID=PROB_SOMATIC,Number=A,Type=Float,Description="P(somatic)">"##,
        );
        header.push_record(br##"##FORMAT=<ID=AF,Number=A,Type=Float,Description="AF">"##);
        header.push_sample(b"tumor");
        header
    }

    /// Run extraction with arbitrary event list against a test VCF.
    fn extract_with_events(
        tmp: &NamedTempFile,
        events: &[&str],
        needs_heatmap: bool,
    ) -> (Vec<RegionSummary>, ExtractionStats) {
        let mut vcf = bcf::Reader::from_path(tmp.path()).unwrap();
        extract_regions(
            &mut vcf,
            "tumor",
            0,
            &events.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
            false,
            needs_heatmap,
        )
        .unwrap()
    }

    /// Run extraction with the single "somatic" event - covers the common case.
    /// Thin wrapper around [`extract_with_events`].
    fn extract_somatic(tmp: &NamedTempFile) -> (Vec<RegionSummary>, ExtractionStats) {
        extract_with_events(tmp, &["somatic"], true)
    }

    /* ====== RegionSummary::from_region_id ======= */

    #[test]
    fn test_from_region_id_valid() {
        let r = RegionSummary::from_region_id("chr1:1000-2000").unwrap();
        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.start, 1000);
        assert!(!r.has_real_indel);
        assert!(r.variants.is_empty());
    }

    #[test]
    fn test_from_region_id_missing_colon() {
        let err = RegionSummary::from_region_id("chr1_1000-2000").unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("chr1_1000-2000"));
        assert!(msg.contains("missing ':'"));
    }

    #[test]
    fn test_from_region_id_missing_dash() {
        let err = RegionSummary::from_region_id("chr1:10002000").unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("chr1:10002000"));
        assert!(msg.contains("missing '-'"));
    }

    #[test]
    fn test_from_region_id_non_numeric_start() {
        let err = RegionSummary::from_region_id("chr1:abc-2000").unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("chr1:abc-2000"));
        assert!(msg.contains("not a valid integer"));
    }

    /* ====== extract_regions tests =============== */

    #[test]
    fn test_non_ms_record_skipped() {
        let tmp = NamedTempFile::new().unwrap();
        let mut w =
            bcf::Writer::from_path(tmp.path(), &make_header(), true, bcf::Format::Vcf).unwrap();

        // Non-MS — no REGION_ID
        let mut r = w.empty_record();
        r.set_rid(Some(0));
        r.set_pos(500);
        r.set_alleles(&[b"A", b"T"]).unwrap();
        r.push_info_float(b"PROB_SOMATIC", &[0.5_f32]).unwrap();
        r.push_format_float(b"AF", &[0.3_f32]).unwrap();
        w.write(&r).unwrap();

        // MS record
        let mut r = w.empty_record();
        r.set_rid(Some(0));
        r.set_pos(999);
        r.set_alleles(&[b"GCAG", b"GCAGCAG"]).unwrap();
        r.push_info_string(MSI_REGION_ID_TAG, &[b"chr1:1000-1030"])
            .unwrap();
        r.push_info_float(b"PROB_SOMATIC", &[0.9_f32]).unwrap();
        r.push_format_float(b"AF", &[0.8_f32]).unwrap();
        w.write(&r).unwrap();
        drop(w);

        let (regions, stats) = extract_somatic(&tmp);
        assert_eq!(stats.total_records, 2);
        assert_eq!(stats.skipped_non_ms, 1);
        assert_eq!(stats.total_ms_regions, 1);
        assert_eq!(regions.len(), 1);
    }

    #[test]
    fn test_two_records_same_region_grouped() {
        // Two records with same REGION_ID - one RegionSummary, two Variants
        let tmp = NamedTempFile::new().unwrap();
        let mut w =
            bcf::Writer::from_path(tmp.path(), &make_header(), true, bcf::Format::Vcf).unwrap();

        for (pos, prob, af) in [(999, 0.9_f32, 0.8_f32), (1002, 0.7, 0.5)] {
            let mut r = w.empty_record();
            r.set_rid(Some(0));
            r.set_pos(pos);
            r.set_alleles(&[b"GCAG", b"GCAGCAG"]).unwrap();
            r.push_info_string(MSI_REGION_ID_TAG, &[b"chr1:1000-1030"])
                .unwrap();
            r.push_info_float(b"PROB_SOMATIC", &[prob]).unwrap();
            r.push_format_float(b"AF", &[af]).unwrap();
            w.write(&r).unwrap();
        }
        drop(w);

        let (regions, stats) = extract_somatic(&tmp);
        assert_eq!(stats.total_ms_regions, 1);
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].variants.len(), 2);
    }

    #[test]
    fn test_prob_absent_computed_correctly() {
        // prob_somatic=0.9 -> prob_absent=0.1
        let tmp = NamedTempFile::new().unwrap();
        let mut w =
            bcf::Writer::from_path(tmp.path(), &make_header(), true, bcf::Format::Vcf).unwrap();

        let mut r = w.empty_record();
        r.set_rid(Some(0));
        r.set_pos(999);
        r.set_alleles(&[b"GCAG", b"GCAGCAG"]).unwrap();
        r.push_info_string(MSI_REGION_ID_TAG, &[b"chr1:1000-1030"])
            .unwrap();
        r.push_info_float(b"PROB_SOMATIC", &[0.9_f32]).unwrap();
        r.push_format_float(b"AF", &[0.8_f32]).unwrap();
        w.write(&r).unwrap();
        drop(w);

        let (regions, _) = extract_somatic(&tmp);
        let v = &regions[0].variants[0];
        assert!((v.prob_absent - 0.1).abs() < TEST_EPSILON_LOOSE);
        assert!((v.af - 0.8).abs() < TEST_EPSILON_LOOSE_F32);
    }

    #[test]
    fn test_dummy_counted_has_real_indel_false() {
        let tmp = NamedTempFile::new().unwrap();
        let mut w =
            bcf::Writer::from_path(tmp.path(), &make_header(), true, bcf::Format::Vcf).unwrap();

        // Deletion half of the dummy pair
        let mut r1 = w.empty_record();
        r1.set_rid(Some(0));
        r1.set_pos(1002);
        r1.set_alleles(&[b"GCAG", b"G"]).unwrap();
        r1.push_info_string(MSI_REGION_ID_TAG, &[b"chr1:1000-1030"])
            .unwrap();
        r1.push_info_flag(MSI_DUMMY_TAG).unwrap();
        r1.push_info_float(b"PROB_SOMATIC", &[0.001_f32]).unwrap();
        r1.push_format_float(b"AF", &[0.0_f32]).unwrap();
        w.write(&r1).unwrap();

        // Insertion half of the dummy pair (swapped REF/ALT, same position/region)
        let mut r2 = w.empty_record();
        r2.set_rid(Some(0));
        r2.set_pos(1002);
        r2.set_alleles(&[b"G", b"GCAG"]).unwrap();
        r2.push_info_string(MSI_REGION_ID_TAG, &[b"chr1:1000-1030"])
            .unwrap();
        r2.push_info_flag(MSI_DUMMY_TAG).unwrap();
        r2.push_info_float(b"PROB_SOMATIC", &[0.001_f32]).unwrap();
        r2.push_format_float(b"AF", &[0.0_f32]).unwrap();
        w.write(&r2).unwrap();

        drop(w);

        let (regions, stats) = extract_somatic(&tmp);
        assert_eq!(stats.dummy_records, 2);
        assert_eq!(regions.len(), 1);
        assert!(!regions[0].has_real_indel);
        assert_eq!(
            regions[0].variants.len(),
            2,
            "both dummy variants should be present"
        );
        assert!((regions[0].variants[0].prob_absent - 0.999).abs() < TEST_EPSILON_LOOSE);
        assert!((regions[0].variants[1].prob_absent - 0.999).abs() < TEST_EPSILON_LOOSE);
    }

    #[test]
    fn test_real_record_sets_has_real_indel() {
        let tmp = NamedTempFile::new().unwrap();
        let mut w =
            bcf::Writer::from_path(tmp.path(), &make_header(), true, bcf::Format::Vcf).unwrap();

        let mut r = w.empty_record();
        r.set_rid(Some(0));
        r.set_pos(999);
        r.set_alleles(&[b"GCAG", b"GCAGCAG"]).unwrap();
        r.push_info_string(MSI_REGION_ID_TAG, &[b"chr1:1000-1030"])
            .unwrap();
        r.push_info_float(b"PROB_SOMATIC", &[0.9_f32]).unwrap();
        r.push_format_float(b"AF", &[0.8_f32]).unwrap();
        w.write(&r).unwrap();
        drop(w);

        let (regions, _) = extract_somatic(&tmp);
        assert!(regions[0].has_real_indel);
    }

    #[test]
    fn test_multiple_regions_sorted_by_chrom_start() {
        let tmp = NamedTempFile::new().unwrap();
        let mut w =
            bcf::Writer::from_path(tmp.path(), &make_header(), true, bcf::Format::Vcf).unwrap();

        for (rid, pos, region_id) in [
            (0, 2002, b"chr1:2000-2030" as &[u8]),
            (0, 999, b"chr1:1000-1030"),
            (1, 499, b"chr2:500-530"),
        ] {
            let mut r = w.empty_record();
            r.set_rid(Some(rid));
            r.set_pos(pos);
            r.set_alleles(&[b"GCAG", b"GCAGCAG"]).unwrap();
            r.push_info_string(MSI_REGION_ID_TAG, &[region_id]).unwrap();
            r.push_info_float(b"PROB_SOMATIC", &[0.9_f32]).unwrap();
            r.push_format_float(b"AF", &[0.8_f32]).unwrap();
            w.write(&r).unwrap();
        }
        drop(w);

        let (regions, stats) = extract_somatic(&tmp);
        assert_eq!(stats.total_ms_regions, 3);
        assert_eq!(regions[0].chrom, "chr1");
        assert_eq!(regions[0].start, 1000);
        assert_eq!(regions[1].chrom, "chr1");
        assert_eq!(regions[1].start, 2000);
        assert_eq!(regions[2].chrom, "chr2");
        assert_eq!(regions[2].start, 500);
    }

    #[test]
    fn test_non_contiguous_same_region_grouped_via_hashmap() {
        // region1 / region2 / region1 - HashMap handles non-contiguous correctly
        let tmp = NamedTempFile::new().unwrap();
        let mut w =
            bcf::Writer::from_path(tmp.path(), &make_header(), true, bcf::Format::Vcf).unwrap();

        for (pos, region_id) in [
            (999, b"chr1:1000-1030" as &[u8]),
            (1999, b"chr1:2000-2030"),
            (1002, b"chr1:1000-1030"), // same as first - non-contiguous
        ] {
            let mut r = w.empty_record();
            r.set_rid(Some(0));
            r.set_pos(pos);
            r.set_alleles(&[b"GCAG", b"GCAGCAG"]).unwrap();
            r.push_info_string(MSI_REGION_ID_TAG, &[region_id]).unwrap();
            r.push_info_float(b"PROB_SOMATIC", &[0.9_f32]).unwrap();
            r.push_format_float(b"AF", &[0.8_f32]).unwrap();
            w.write(&r).unwrap();
        }
        drop(w);

        let (regions, stats) = extract_somatic(&tmp);
        assert_eq!(stats.total_ms_regions, 2);
        assert_eq!(regions[0].chrom, "chr1");
        assert_eq!(regions[0].start, 1000);
        assert_eq!(regions[0].variants.len(), 2);
        assert_eq!(regions[1].chrom, "chr1");
        assert_eq!(regions[1].start, 2000);
        assert_eq!(regions[1].variants.len(), 1);
    }

    #[test]
    fn test_two_alts_different_regions_each_registered() {
        // One record, two ALTs, each matching a different region via Number=A
        // ALT0: ACAGCAG -> indel_pos=1001 -> chr1:1001-1030
        // ALT1: A       -> indel_pos=998  -> chr1:992-1000
        let tmp = NamedTempFile::new().unwrap();
        let mut w =
            bcf::Writer::from_path(tmp.path(), &make_header(), true, bcf::Format::Vcf).unwrap();

        let mut r = w.empty_record();
        r.set_rid(Some(0));
        r.set_pos(997);
        r.set_alleles(&[b"GCAG", b"GCAGCAG", b"G"]).unwrap();
        r.push_info_string(MSI_REGION_ID_TAG, &[b"chr1:1001-1030", b"chr1:992-1000"])
            .unwrap();
        r.push_info_float(b"PROB_SOMATIC", &[0.9_f32, 0.7_f32])
            .unwrap();
        r.push_format_float(b"AF", &[0.6_f32, 0.4_f32]).unwrap();
        w.write(&r).unwrap();
        drop(w);

        let (regions, stats) = extract_somatic(&tmp);

        assert_eq!(stats.total_ms_regions, 2);
        assert_eq!(regions.len(), 2);

        // After sort by (chrom, start): chr1:992-1000 first
        assert_eq!(regions[0].chrom, "chr1");
        assert_eq!(regions[0].start, 992);
        assert_eq!(regions[1].chrom, "chr1");
        assert_eq!(regions[1].start, 1001);

        assert_eq!(regions[0].variants.len(), 1);
        assert_eq!(regions[1].variants.len(), 1);

        assert!(regions[0].has_real_indel);
        assert!(regions[1].has_real_indel);

        // ALT1 -> chr1:992-1000: prob_somatic=0.7 -> prob_absent≈0.3
        assert!((regions[0].variants[0].prob_absent - 0.3).abs() < TEST_EPSILON_LOOSE);
        // ALT0 -> chr1:1001-1030: prob_somatic=0.9 -> prob_absent≈0.1
        assert!((regions[1].variants[0].prob_absent - 0.1).abs() < TEST_EPSILON_LOOSE);
    }

    #[test]
    fn test_alt_with_dot_region_skipped() {
        // Two ALTs: ALT0 has region, ALT1 has "." - ALT1 skipped
        let tmp = NamedTempFile::new().unwrap();
        let mut w =
            bcf::Writer::from_path(tmp.path(), &make_header(), true, bcf::Format::Vcf).unwrap();

        let mut r = w.empty_record();
        r.set_rid(Some(0));
        r.set_pos(997);
        r.set_alleles(&[b"ACAG", b"ACAGCAG", b"A"]).unwrap();
        r.push_info_string(MSI_REGION_ID_TAG, &[b"chr1:1001-1030", b"."])
            .unwrap();
        r.push_info_float(b"PROB_SOMATIC", &[0.9_f32, 0.7_f32])
            .unwrap();
        r.push_format_float(b"AF", &[0.6_f32, 0.4_f32]).unwrap();
        w.write(&r).unwrap();
        drop(w);

        let (regions, stats) = extract_somatic(&tmp);
        assert_eq!(stats.total_ms_regions, 1);
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].chrom, "chr1");
        assert_eq!(regions[0].start, 1001);
        assert_eq!(regions[0].variants.len(), 1);
        assert!((regions[0].variants[0].prob_absent - 0.1).abs() < TEST_EPSILON_LOOSE);
    }

    #[test]
    fn test_missing_af_skipped_counted_in_stats() {
        // Record has REGION_ID and PROB_SOMATIC but no FORMAT:AF
        // Region is registered but variant is skipped - skipped_missing_af increments
        let tmp = NamedTempFile::new().unwrap();
        let mut w =
            bcf::Writer::from_path(tmp.path(), &make_header(), true, bcf::Format::Vcf).unwrap();

        let mut r = w.empty_record();
        r.set_rid(Some(0));
        r.set_pos(999);
        r.set_alleles(&[b"GCAG", b"GCAGCAG"]).unwrap();
        r.push_info_string(MSI_REGION_ID_TAG, &[b"chr1:1000-1030"])
            .unwrap();
        r.push_info_float(b"PROB_SOMATIC", &[0.9_f32]).unwrap();
        r.push_format_float(b"AF", &[f32::missing()]).unwrap(); // missing AF
        w.write(&r).unwrap();
        drop(w);

        let (regions, stats) = extract_somatic(&tmp);

        assert_eq!(stats.skipped_missing_af, 1);
        // Region is registered (REGION_ID was present) but has no variants
        assert_eq!(stats.total_ms_regions, 1);
        assert_eq!(regions[0].variants.len(), 0);
        // has_real_indel still set - the record existed regardless of AF
        assert!(regions[0].has_real_indel);
    }

    #[test]
    fn test_missing_prob_skipped_counted_in_stats() {
        let tmp = NamedTempFile::new().unwrap();
        let mut w =
            bcf::Writer::from_path(tmp.path(), &make_header(), true, bcf::Format::Vcf).unwrap();

        let mut r = w.empty_record();
        r.set_rid(Some(0));
        r.set_pos(999);
        r.set_alleles(&[b"GCAG", b"GCAGCAG"]).unwrap();
        r.push_info_string(MSI_REGION_ID_TAG, &[b"chr1:1000-1030"])
            .unwrap();
        r.push_info_float(b"PROB_SOMATIC", &[f32::missing()])
            .unwrap();
        r.push_format_float(b"AF", &[0.8_f32]).unwrap();
        w.write(&r).unwrap();
        drop(w);

        let (regions, stats) = extract_somatic(&tmp);

        assert_eq!(stats.skipped_missing_prob, 1);
        assert_eq!(stats.total_ms_regions, 1);
        assert_eq!(regions[0].variants.len(), 0);
        assert!(regions[0].has_real_indel);
    }

    #[test]
    fn test_two_event_prob_combined_correctly() {
        let tmp = NamedTempFile::new().unwrap();
        let mut header = make_header();
        header.push_record(
            br##"##INFO=<ID=PROB_HIGH_VAF,Number=A,Type=Float,Description="P(high_vaf)">"##,
        );
        let mut w = bcf::Writer::from_path(tmp.path(), &header, true, bcf::Format::Vcf).unwrap();

        let mut r = w.empty_record();
        r.set_rid(Some(0));
        r.set_pos(999);
        r.set_alleles(&[b"GCAG", b"GCAGCAG"]).unwrap();
        r.push_info_string(MSI_REGION_ID_TAG, &[b"chr1:1000-1030"])
            .unwrap();
        r.push_info_float(b"PROB_SOMATIC", &[0.3_f32]).unwrap();
        r.push_info_float(b"PROB_HIGH_VAF", &[0.4_f32]).unwrap();
        r.push_format_float(b"AF", &[0.6_f32]).unwrap();
        w.write(&r).unwrap();
        drop(w);

        let (regions, _) = extract_with_events(&tmp, &["somatic", "high_vaf"], true);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].variants.len(), 1);

        // ln_sum_exp(ln(0.3), ln(0.4)) -> prob_events ≈ 0.7, 1 - prob_events ≈ 0.3
        assert!(
            (regions[0].variants[0].prob_absent - 0.3).abs() < TEST_EPSILON_LOOSE,
            "expected prob_absent ≈ 0.3, got {}",
            regions[0].variants[0].prob_absent
        );
    }

    #[test]
    fn test_variant_less_region_pruned_when_heatmap_not_needed() {
        let tmp = NamedTempFile::new().unwrap();
        let mut w =
            bcf::Writer::from_path(tmp.path(), &make_header(), true, bcf::Format::Vcf).unwrap();

        let mut r = w.empty_record();
        r.set_rid(Some(0));
        r.set_pos(999);
        r.set_alleles(&[b"GCAG", b"GCAGCAG"]).unwrap();
        r.push_info_string(MSI_REGION_ID_TAG, &[b"chr1:1000-1030"])
            .unwrap();
        r.push_info_float(b"PROB_SOMATIC", &[0.9_f32]).unwrap();
        r.push_format_float(b"AF", &[f32::missing()]).unwrap(); // missing AF -> zero variants
        w.write(&r).unwrap();
        drop(w);

        let (regions, stats) = extract_with_events(&tmp, &["somatic"], false);

        assert_eq!(
            stats.total_ms_regions, 1,
            "denominator unaffected by pruning"
        );
        assert_eq!(
            regions.len(),
            0,
            "variant-less region pruned when heatmap not requested"
        );
    }

    #[test]
    fn test_regions_with_real_indel_counted_even_when_pruned() {
        let tmp = NamedTempFile::new().unwrap();
        let mut w =
            bcf::Writer::from_path(tmp.path(), &make_header(), true, bcf::Format::Vcf).unwrap();

        let mut r = w.empty_record();
        r.set_rid(Some(0));
        r.set_pos(999);
        r.set_alleles(&[b"GCAG", b"GCAGCAG"]).unwrap();
        r.push_info_string(MSI_REGION_ID_TAG, &[b"chr1:1000-1030"])
            .unwrap();
        r.push_info_float(b"PROB_SOMATIC", &[0.9_f32]).unwrap();
        r.push_format_float(b"AF", &[f32::missing()]).unwrap();
        w.write(&r).unwrap();
        drop(w);

        let (regions, stats) = extract_with_events(&tmp, &["somatic"], false);

        assert!(regions.is_empty(), "region was pruned");
        assert_eq!(
            stats.regions_with_real_indel, 1,
            "diagnostic count survives pruning"
        );
    }
}
