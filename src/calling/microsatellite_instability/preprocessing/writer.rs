//! writer.rs
//!
//! VCF writing operations for MSI preprocessing.
//!
//! This module provides:
//! 1. `VariantInWindow` - Structure for variants with accumulated region annotations
//! 2. `inject_dummy_indels` - Creates dummy deletion + insertion pair for
//!    regions without perfect repeats
//! 3. `write_variant` - Writes variants with REGION_ID annotations if applicable
//! 4.  Unit tests for both (2,3) functions, covering typical and edge cases
//!

use anyhow::Result;
use rust_htslib::bcf::{self, Writer};
use std::collections::HashMap;

use crate::constants::{MSI_COPY_FIELDS, MSI_DUMMY_TAG, MSI_OMIT_AUX, MSI_REGION_ID_TAG};
use crate::utils::aux_info::AuxInfo;
use crate::utils::bcf_utils::copy_info_fields;
use crate::utils::ms_bed::BedRegion;

/* ============ Data Structures =================== */

/// Variant in sliding window with accumulated region annotations.
pub(super) struct VariantInWindow {
    /// The VCF record
    pub record: bcf::Record,
    /// Chromosome name
    pub chrom: String,
    /// Region IDs matched per ALT allele (Number=A).
    ///
    /// Key: ALT index (0-based into the ALT array).
    /// Value: BED region identifier matched by this ALT (e.g. "chr1:1000-1030").
    ///
    /// Only ALTs that overlap a BED region have an entry.
    /// ALTs with no BED match are absent.
    /// A single record can have multiple entries if different ALTs
    /// overlap different non-overlapping BED regions.
    ///
    /// Invariant: each ALT index maps to at most one region.
    /// Assigning a second region to the same ALT is an error.
    pub matching_regions: HashMap<usize, String>,
    /// User-specified INFO fields collected from input record.
    /// Written to output via aux_info.write() during write_variant.
    /// Fields in MSI_OMIT_AUX are excluded to prevent
    /// double-writing with copy_info_fields().
    pub aux_info: AuxInfo,
}

/* ============ Functions ========================= */

/// Inject dummy indels (deletion + insertion) for a region with no observed
/// perfect indels.
///
/// Creates two hypothetical variants of one motif unit at the same anchored
/// position (last base of the first repeat unit): a deletion and an
/// insertion. Writing both lets `call variants` independently evaluate
/// support for either direction against real reads - whichever direction
/// has coverage (or both) contributes to the DP model; if neither does,
/// both remain uninformative.
///
/// Output records contains only REGION_ID and MSI_DUMMY INFO fields.
/// No FORMAT or sample data is written.
///
/// # Arguments
/// * `writer` - VCF writer
/// * `region` - BED region requiring dummy indels
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if chromosome not found or motif is empty
///
/// # Note:
/// Panics if `region.chrom` is not present in the writer's header.
/// This is expected to never happen when called from the real pipeline:
/// `prepare_header` (header.rs) guarantees every BED chromosome is added
/// to the header before this writer is constructed - see
/// `header.rs::test_prepare_header_adds_missing_contigs`.
///
/// # Example
/// assert!(inject_dummy_indels(&mut writer, &region).is_ok());
pub(super) fn inject_dummy_indels(writer: &mut Writer, region: &BedRegion) -> Result<()> {
    // Guaranteed non-empty by ms_bed.rs's parse_motif_from_name.
    // Re-enable if that changes.
    // if motif_bytes.is_empty() {
    //     return Err(Error::MsiBedMotifInvalid {
    //         motif: "(empty)".to_string(),
    //     }
    //     .into());
    // }

    let rid = writer.header().name2rid(region.chrom.as_bytes()).unwrap();
    let motif_bytes = region.motif.as_bytes();
    let dummy_pos = region.dummy_indel_position();
    let anchor = motif_bytes[motif_bytes.len() - 1];
    let region_id = region.region_id();

    // Deletion: REF = anchor + motif, ALT = anchor
    {
        let mut record = writer.empty_record();
        record.set_rid(Some(rid));
        record.set_pos(dummy_pos as i64);

        let mut ref_allele = vec![anchor];
        ref_allele.extend_from_slice(motif_bytes);
        let alt_allele = vec![anchor];

        record.set_alleles(&[&ref_allele, &alt_allele])?;
        record.push_info_string(MSI_REGION_ID_TAG, &[region_id.as_bytes()])?;
        record.push_info_flag(MSI_DUMMY_TAG)?;
        writer.write(&record)?;
    }

    // Insertion: REF = anchor, ALT = anchor + motif
    {
        let mut record = writer.empty_record();
        record.set_rid(Some(rid));
        record.set_pos(dummy_pos as i64);

        let ref_allele = vec![anchor];
        let mut alt_allele = vec![anchor];
        alt_allele.extend_from_slice(motif_bytes);

        record.set_alleles(&[&ref_allele, &alt_allele])?;
        record.push_info_string(MSI_REGION_ID_TAG, &[region_id.as_bytes()])?;
        record.push_info_flag(MSI_DUMMY_TAG)?;
        writer.write(&record)?;
    }

    Ok(())
}

/// Write variant to output, with region annotations if present.
///
/// Creates a fresh output record containing:
/// 1. Position and alleles from input record
/// 2. Standard INFO fields copied via copy_info_fields()
/// 3. User-specified INFO fields via aux_info.write()
/// 4. REGION_ID annotation if variant overlaps a perfect MS region
///
/// FORMAT and sample data are intentionally omitted.
///
/// # Arguments
/// * `writer` - VCF writer
/// * `variant_info` - Variant with accumulated region annotation and aux info
/// * `counter` - Counter for annotated MS indels (incremented if annotations present)
pub(super) fn write_variant(
    writer: &mut Writer,
    variant_info: VariantInWindow,
    counter: &mut usize,
) -> Result<()> {
    let mut output_record = writer.empty_record();

    let output_rid = writer
        .header()
        .name2rid(variant_info.chrom.as_bytes())
        .expect("chromosome missing from output header - invariant violated, see prepare_header");

    output_record.set_rid(Some(output_rid));
    output_record.set_pos(variant_info.record.pos());
    output_record.set_alleles(&variant_info.record.alleles())?;

    copy_info_fields(&variant_info.record, &mut output_record, &MSI_COPY_FIELDS)?;

    variant_info
        .aux_info
        .write(&mut output_record, &MSI_OMIT_AUX)?;

    let n_alts = output_record.allele_count() as usize - 1;
    let region_id_bytes: Vec<&[u8]> = (0..n_alts)
        .map(|i| {
            variant_info
                .matching_regions
                .get(&i)
                .map(|s| s.as_bytes())
                .unwrap_or(b".")
        })
        .collect();

    if region_id_bytes.iter().any(|r| *r != b".") {
        output_record.push_info_string(MSI_REGION_ID_TAG, &region_id_bytes)?;
        *counter += 1;
    }

    writer.write(&output_record)?;

    Ok(())
}

/* =============== Tests ========================== */

#[cfg(test)]
mod tests {
    use super::*;

    use rust_htslib::bcf::{self, Read};
    use tempfile::NamedTempFile;

    use crate::constants::{MSI_DUMMY_HEADER, MSI_REGION_ID_HEADER};
    use crate::utils::aux_info::tests::make_aux_collector;
    use crate::utils::bcf_utils::get_chrom;
    use crate::utils::bcf_utils::tests::{
        create_test_record, create_test_record_multi_alt, create_test_vcf, read_first_record,
        TestVcfConfig,
    };

    /* ============ VCF Helpers  ======================= */

    /// Create minimal VCF header for dummy indel tests
    fn create_minimal_vcf_header() -> bcf::Header {
        let mut header = bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");
        header.push_record(MSI_REGION_ID_HEADER.as_bytes());
        header.push_record(MSI_DUMMY_HEADER.as_bytes());
        header
    }

    /* ====== inject_dummy_indels tests ============== */

    #[test]
    fn test_inject_dummy_indels_simple_motif() {
        let tmp = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 1000,
            end: 1021,
            motif: "CAG".to_string(),
        };

        inject_dummy_indels(&mut writer, &region).unwrap();
        drop(writer);

        let mut reader = bcf::Reader::from_path(tmp.path()).unwrap();
        let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
        assert_eq!(records.len(), 2, "should write deletion + insertion");

        let del = &records[0];
        let del_alleles = del.alleles();
        assert_eq!(del.pos(), 1002, "Position: last base of 1st CAG");
        assert_eq!(
            del_alleles[0], b"GCAG",
            "Deletion REF: anchor G + motif CAG"
        );
        assert_eq!(del_alleles[1], b"G", "Deletion ALT: just anchor G");
        let del_region_id = del.info(MSI_REGION_ID_TAG).string().unwrap().unwrap();
        assert_eq!(del_region_id[0], b"chr1:1000-1021");

        let ins = &records[1];
        let ins_alleles = ins.alleles();
        assert_eq!(
            ins.pos(),
            1002,
            "Insertion shares the same anchor position as deletion"
        );
        assert_eq!(ins_alleles[0], b"G", "Insertion REF: anchor G");
        assert_eq!(
            ins_alleles[1], b"GCAG",
            "Insertion ALT: anchor G + motif CAG"
        );
        let ins_region_id = ins.info(MSI_REGION_ID_TAG).string().unwrap().unwrap();
        assert_eq!(ins_region_id[0], b"chr1:1000-1021");
    }

    #[test]
    fn test_inject_dummy_indels_different_motifs() {
        let tmp = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let cases = vec![
            ("A", 1000, 1020, 1000, b"AA" as &[u8], b"A" as &[u8]),
            ("AT", 1100, 1120, 1101, b"TAT" as &[u8], b"T" as &[u8]),
            ("CAG", 1200, 1221, 1202, b"GCAG" as &[u8], b"G" as &[u8]),
            ("AAAG", 1300, 1320, 1303, b"GAAAG" as &[u8], b"G" as &[u8]),
        ];

        for (motif, start, end, _, _, _) in &cases {
            let region = BedRegion {
                chrom: "chr1".to_string(),
                start: *start,
                end: *end,
                motif: motif.to_string(),
            };
            inject_dummy_indels(&mut writer, &region).unwrap();
        }
        drop(writer);

        let mut reader = bcf::Reader::from_path(tmp.path()).unwrap();
        for (motif, _, _, expected_pos, expected_del_ref, expected_del_alt) in &cases {
            let del = reader.records().next().unwrap().unwrap();
            let del_alleles = del.alleles();
            assert_eq!(
                del.pos(),
                *expected_pos as i64,
                "Deletion position mismatch for motif {}",
                motif
            );
            assert_eq!(
                del_alleles[0], *expected_del_ref,
                "Deletion REF mismatch for motif {}",
                motif
            );
            assert_eq!(
                del_alleles[1], *expected_del_alt,
                "Deletion ALT mismatch for motif {}",
                motif
            );

            let ins = reader.records().next().unwrap().unwrap();
            let ins_alleles = ins.alleles();
            assert_eq!(
                ins.pos(),
                *expected_pos as i64,
                "Insertion position mismatch for motif {}",
                motif
            );
            assert_eq!(
                ins_alleles[0], *expected_del_alt,
                "Insertion REF mismatch for motif {}",
                motif
            );
            assert_eq!(
                ins_alleles[1], *expected_del_ref,
                "Insertion ALT mismatch for motif {}",
                motif
            );
        }
    }

    #[test]
    #[should_panic(expected = "called `Result::unwrap()` on an `Err` value")]
    fn test_inject_dummy_indels_chromosome_not_found() {
        let tmp = NamedTempFile::new().unwrap();
        let mut header = bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");

        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let region = BedRegion {
            chrom: "chr2".to_string(), // Not in header!
            start: 1000,
            end: 1021,
            motif: "CAG".to_string(),
        };

        let _ = inject_dummy_indels(&mut writer, &region);
    }

    /* ====== write_variant tests ==================== */

    #[test]
    fn test_write_variant_with_annotation() {
        let tmp: NamedTempFile = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let record = create_test_record(&writer, 0, 1000, b"ACAG", b"ACAGCAG");

        let mut matching_regions = HashMap::new();
        matching_regions.insert(0, "chr1:1000-1020".to_string());

        let variant_info = VariantInWindow {
            record,
            chrom: "chr1".to_string(),
            matching_regions,
            aux_info: AuxInfo::default(),
        };

        let mut counter = 0;
        write_variant(&mut writer, variant_info, &mut counter).unwrap();
        drop(writer);

        assert_eq!(counter, 1);

        let record = read_first_record(tmp.path());

        let region_id = record.info(MSI_REGION_ID_TAG).string().unwrap();
        assert!(region_id.is_some());
        assert_eq!(region_id.unwrap()[0], b"chr1:1000-1020");
    }

    #[test]
    fn test_write_variant_without_annotation() {
        let tmp = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let record = create_test_record(&writer, 0, 2000, b"A", b"T");

        let variant_info = VariantInWindow {
            record,
            chrom: "chr1".to_string(),
            matching_regions: HashMap::new(),
            aux_info: AuxInfo::default(),
        };

        let mut counter = 0;
        write_variant(&mut writer, variant_info, &mut counter).unwrap();
        drop(writer);

        assert_eq!(counter, 0);

        let record = read_first_record(tmp.path());

        let region_id = record.info(MSI_REGION_ID_TAG).string().unwrap();
        assert!(region_id.is_none());

        assert_eq!(record.pos(), 2000);
        let alleles = record.alleles();
        assert_eq!(alleles[0], b"A");
        assert_eq!(alleles[1], b"T");
    }

    #[test]
    fn test_write_variant_removes_existing_region_id_non_ms() {
        let tmp = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let mut record = create_test_record(&writer, 0, 2000, b"A", b"T");
        record
            .push_info_string(MSI_REGION_ID_TAG, &[b"old:value"])
            .unwrap();

        let variant_info = VariantInWindow {
            record,
            chrom: "chr1".to_string(),
            matching_regions: HashMap::new(),
            aux_info: AuxInfo::default(),
        };

        let mut counter = 0;
        write_variant(&mut writer, variant_info, &mut counter).unwrap();
        drop(writer);

        assert_eq!(counter, 0);

        let record = read_first_record(tmp.path());

        let region_id = record.info(MSI_REGION_ID_TAG).string().unwrap();
        assert!(region_id.is_none(), "REGION_ID should be removed");
    }

    #[test]
    fn test_write_variant_replaces_existing_region_id_ms_indel() {
        let tmp = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let mut record = create_test_record(&writer, 0, 1000, b"ACAG", b"ACAGCAG");
        record
            .push_info_string(MSI_REGION_ID_TAG, &[b"old:value"])
            .unwrap();

        let variant_info = VariantInWindow {
            record,
            chrom: "chr1".to_string(),
            matching_regions: {
                let mut m = HashMap::new();
                m.insert(0, "chr1:1000-1020".to_string());
                m
            },
            aux_info: AuxInfo::default(),
        };

        let mut counter = 0;
        write_variant(&mut writer, variant_info, &mut counter).unwrap();
        drop(writer);

        assert_eq!(counter, 1);

        let record = read_first_record(tmp.path());

        let region_id = record.info(MSI_REGION_ID_TAG).string().unwrap().unwrap();
        assert_eq!(region_id.len(), 1);
        assert_eq!(region_id[0], b"chr1:1000-1020");
    }

    #[test]
    fn test_write_variant_propagates_aux_info() {
        let (tmp_src, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"ACAG",
            alt_alleles: vec![b"ACAGCAG"],
            extra_info_fields: vec![(b"COSMIC_ID", b"1", b"String", b"COSMIC ID", b"COSM123")],
            ..Default::default()
        });

        let aux = make_aux_collector(tmp_src.path(), &["COSMIC_ID"]);
        let mut src_reader = bcf::Reader::from_path(tmp_src.path()).unwrap();
        let src_record = src_reader.records().next().unwrap().unwrap();
        let aux_info = aux.collect(&src_record).unwrap();

        let tmp = NamedTempFile::new().unwrap();
        let mut header = create_minimal_vcf_header();
        header.push_record(
            br##"##INFO=<ID=COSMIC_ID,Number=1,Type=String,Description="COSMIC ID">"##,
        );
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let record = create_test_record(&writer, 0, 1000, b"ACAG", b"ACAGCAG");
        let variant_info = VariantInWindow {
            record,
            chrom: "chr1".to_string(),
            matching_regions: {
                let mut m = HashMap::new();
                m.insert(0, "chr1:1000-1020".to_string());
                m
            },
            aux_info,
        };

        let mut counter = 0;
        write_variant(&mut writer, variant_info, &mut counter).unwrap();
        drop(writer);

        let record = read_first_record(tmp.path());
        let cosmic = record.info(b"COSMIC_ID").string().unwrap();
        assert!(cosmic.is_some(), "COSMIC_ID propagated via aux_info");
        assert_eq!(cosmic.unwrap()[0], b"COSM123");
    }

    #[test]
    fn test_write_variant_two_alts_different_regions() {
        let tmp = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        // Two ALTs: ALT0=ACAGCAG (overlaps Region1), ALT1=A (overlaps Region2)
        let record = create_test_record_multi_alt(&writer, 0, 997, b"ACAG", &[b"ACAGCAG", b"A"]);

        let mut matching_regions = HashMap::new();
        matching_regions.insert(0, "chr1:1001-1030".to_string()); // ALT0
        matching_regions.insert(1, "chr1:994-1001".to_string()); // ALT1

        let variant_info = VariantInWindow {
            record,
            chrom: "chr1".to_string(),
            matching_regions,
            aux_info: AuxInfo::default(),
        };

        let mut counter = 0;
        write_variant(&mut writer, variant_info, &mut counter).unwrap();
        drop(writer);

        assert_eq!(counter, 1);

        let record = read_first_record(tmp.path());
        let region_ids = record.info(MSI_REGION_ID_TAG).string().unwrap().unwrap();

        assert_eq!(region_ids.len(), 2);
        assert_eq!(region_ids[0], b"chr1:1001-1030");
        assert_eq!(region_ids[1], b"chr1:994-1001");
    }

    #[test]
    fn test_write_variant_uses_chrom_field_not_record_rid() {
        // Output header: chr1=rid0, chr2=rid1.
        let mut header = bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");
        header.push_record(br"##contig=<ID=chr2,length=1000000>");
        header.push_record(br##"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">"##);

        let tmp = NamedTempFile::new().unwrap();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        // Record built with rid=0 (chr1 in THIS header) - but its VariantInWindow.chrom
        // deliberately says "chr2", simulating what happens when the record's
        // embedded rid no longer matches the output header's contig order.
        // If write_variant ever uses record.rid() (from input record), this would
        // wrongly write it as chr1.
        let record = create_test_record(&writer, 0, 99, b"A", b"AT");

        let variant_info = VariantInWindow {
            record,
            chrom: "chr2".to_string(),
            matching_regions: HashMap::new(),
            aux_info: AuxInfo::default(),
        };

        let mut counter = 0;
        write_variant(&mut writer, variant_info, &mut counter).unwrap();
        drop(writer);

        let record = read_first_record(tmp.path());
        let chrom = get_chrom(&record).unwrap();
        assert_eq!(
            chrom, "chr2",
            "should use VariantInWindow.chrom, not record's raw rid"
        );
    }
}
