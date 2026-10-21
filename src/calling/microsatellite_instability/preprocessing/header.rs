//! header.rs
//!
//! Output VCF header preparation for MSI preprocessing.
//!
//! This module provides:
//! 1. `prepare_header` function to create a modified VCF header for MSI preprocessing output

use std::path::Path;

use anyhow::Result;
use rust_htslib::bcf;

use crate::constants::{MSI_COPY_FIELDS, MSI_DUMMY_HEADER, MSI_REGION_ID_HEADER};
use crate::errors::Error;
use crate::utils::aux_info::AuxInfoCollector;
use crate::utils::ms_bed::collect_bed_chromosomes;

/// Prepare output VCF header for MSI preprocessing.
///
/// # Algorithm
/// 1. Collect contig names and the standard INFO field declarations (`MSI_COPY_FIELDS`) from the input header
/// 2. Collect the BED chromosomes and add those missing from the input header
/// 3. Write all contigs in lexicographically sorted order
/// 4. Write the INFO/REGION_ID and INFO/MSI_DUMMY declarations
/// 5. Write the collected standard INFO declarations, then the aux info fields
///
/// # Arguments
/// * `input_header` - HeaderView from input VCF reader
/// * `bed_path` - Path to BED file (for collecting chromosome names)
/// * `aux_info_collector` - Collector for auxiliary INFO fields to propagate (used for header preparation)
///
/// # Returns
/// * `Ok(Header)` - Prepared header ready for writer creation
/// * `Err` if:
///     - BED file cannot be opened or read, or a record in it is invalid
///     - BED file has no valid regions (`Error::BedFileNoValidRegions`)
///     - A copied INFO field lacks Number, Type or Description (`Error::VcfHeaderFieldMissing`)
///
/// # Note
/// Only the contig ID is kept for each contig (length and other attributes are dropped),
/// and header records other than the contig & INFO fields listed in the Algorithm are not carried over.
pub(super) fn prepare_header(
    input_header: &bcf::header::HeaderView,
    bed_path: &Path,
    aux_info_collector: &AuxInfoCollector,
) -> Result<bcf::Header> {
    let mut header = bcf::Header::new();

    // Copy Contigs & standard INFO field declarations from input header into vectors for later processing
    let mut all_contigs: Vec<String> = Vec::new();
    let mut standard_info_fields: Vec<(String, String, String, String)> = Vec::new();

    for rec in input_header.header_records() {
        match rec {
            bcf::header::HeaderRecord::Contig { values, .. } => {
                if let Some(id) = values.get("ID") {
                    all_contigs.push(id.clone());
                }
            }
            bcf::header::HeaderRecord::Info { values, .. } => {
                if let Some(id) = values.get("ID") {
                    if MSI_COPY_FIELDS.contains(&id.as_str()) {
                        match (
                            values.get("Number"),
                            values.get("Type"),
                            values.get("Description"),
                        ) {
                            (Some(number), Some(type_), Some(desc)) => {
                                standard_info_fields.push((
                                    id.clone(),
                                    number.clone(),
                                    type_.clone(),
                                    desc.clone(),
                                ));
                            }
                            _ => {
                                // Number, Type and Description are required by the VCF spec
                                return Err(Error::VcfHeaderFieldMissing {
                                    field: id.clone(),
                                    location: "INFO".to_string(),
                                }
                                .into());
                            }
                        }
                    }
                }
            }
            _ => {}
        }
    }

    // Collect BED chromosomes missing from input
    let bed_chroms = collect_bed_chromosomes(bed_path)?;

    if bed_chroms.is_empty() {
        return Err(Error::BedFileNoValidRegions.into());
    }

    // Only the BED chromosomes absent from the input header are
    // added - merged with the input's own contigs BEFORE sorting,
    // so the final list of contigs in header is in sorted order.
    let missing_bed_chroms: Vec<String> = bed_chroms
        .iter()
        .filter(|chrom| input_header.name2rid(chrom.as_bytes()).is_err())
        .cloned()
        .collect();

    all_contigs.extend(missing_bed_chroms);
    all_contigs.sort();

    for id in &all_contigs {
        header.push_record(format!("##contig=<ID={}>", id).as_bytes());
    }

    // Add MSI-specific INFO fields
    header.push_record(MSI_REGION_ID_HEADER.as_bytes());
    header.push_record(MSI_DUMMY_HEADER.as_bytes());

    // Standard INFO field declarations collected above.
    for (id, number, type_, desc) in &standard_info_fields {
        header.push_record(
            format!(
                "##INFO=<ID={},Number={},Type={},Description={}>",
                id, number, type_, desc
            )
            .as_bytes(),
        );
    }

    // User-specified propagated fields via --propagate-info-fields
    aux_info_collector.write_header_info(&mut header);

    Ok(header)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    use rust_htslib::bcf::{Format, Read};
    use tempfile::NamedTempFile;

    use crate::constants::MSI_REGION_ID_TAG;
    use crate::utils::aux_info::tests::make_aux_collector;

    /// Helper: Create minimal test VCF file with given contigs and INFO header lines
    fn create_test_vcf(contigs: &[&str], info_records: &[&str]) -> NamedTempFile {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "##fileformat=VCFv4.2").unwrap();
        for contig in contigs {
            writeln!(tmp, "##contig=<ID={}>", contig).unwrap();
        }
        for info in info_records {
            writeln!(tmp, "{}", info).unwrap();
        }
        writeln!(tmp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        tmp.flush().unwrap();
        tmp
    }

    /// Helper: Create test BED file with given chromosomes
    /// NOTE: Each region is 100-121 with a 7xCAG motif
    /// (inappropriate for MSI but good for testing header prep)
    fn create_test_bed(chroms: &[&str]) -> NamedTempFile {
        let mut tmp = NamedTempFile::new().unwrap();
        for chrom in chroms {
            writeln!(tmp, "{}\t100\t121\t7xCAG", chrom).unwrap();
        }
        tmp.flush().unwrap();
        tmp
    }

    /// Helper: Write header to file and read back as string for assertions
    fn header_to_string(header: &bcf::Header) -> String {
        let tmp = NamedTempFile::new().unwrap();
        let writer = bcf::Writer::from_path(tmp.path(), header, true, Format::Vcf).unwrap();
        drop(writer);
        std::fs::read_to_string(tmp.path()).unwrap()
    }

    // Helper: Field declaration fragment for REGION_ID
    fn region_id_field() -> String {
        format!("ID={}", String::from_utf8_lossy(MSI_REGION_ID_TAG))
    }

    /* ====== prepare_header tests =================== */

    #[test]
    fn test_prepare_header_adds_region_id_field() {
        let vcf_file = create_test_vcf(&["chr1"], &[]);
        let bed_file = create_test_bed(&["chr1"]);
        let aux = make_aux_collector(vcf_file.path(), &[]);

        let reader = bcf::Reader::from_path(vcf_file.path()).unwrap();
        let result = prepare_header(reader.header(), bed_file.path(), &aux).unwrap();

        let content = header_to_string(&result);
        assert!(content.contains(&region_id_field()));
    }

    #[test]
    fn test_prepare_header_adds_missing_contigs() {
        let vcf_file = create_test_vcf(&["chr1"], &[]);
        let bed_file = create_test_bed(&["chr1", "chr2"]);
        let aux = make_aux_collector(vcf_file.path(), &[]);

        let reader = bcf::Reader::from_path(vcf_file.path()).unwrap();
        let result = prepare_header(reader.header(), bed_file.path(), &aux).unwrap();

        let content = header_to_string(&result);
        assert!(content.contains(&region_id_field()));
        assert!(content.contains("contig=<ID=chr2>"));
    }

    #[test]
    fn test_prepare_header_errors_on_empty_bed() {
        let vcf_file = create_test_vcf(&["chr1"], &[]);
        let empty_bed = NamedTempFile::new().unwrap();
        let aux = make_aux_collector(vcf_file.path(), &[]);

        let reader = bcf::Reader::from_path(vcf_file.path()).unwrap();
        let err = prepare_header(reader.header(), empty_bed.path(), &aux).unwrap_err();

        assert_eq!(
            err.downcast_ref::<Error>(),
            Some(&Error::BedFileNoValidRegions)
        );
    }

    #[test]
    fn test_prepare_header_preserves_existing_contigs() {
        let vcf_file = create_test_vcf(&["chr1", "chr2"], &[]);
        let bed_file = create_test_bed(&["chr1"]);
        let aux = make_aux_collector(vcf_file.path(), &[]);

        let reader = bcf::Reader::from_path(vcf_file.path()).unwrap();
        let result = prepare_header(reader.header(), bed_file.path(), &aux).unwrap();

        let content = header_to_string(&result);
        // chr2 should still be in header (not removed)
        assert!(content.contains("contig=<ID=chr1>"));
        assert!(content.contains("contig=<ID=chr2>"));
    }

    #[test]
    fn test_prepare_header_no_duplicate_contigs() {
        let vcf_file = create_test_vcf(&["chr1"], &[]);
        let bed_file = create_test_bed(&["chr1"]);
        let aux = make_aux_collector(vcf_file.path(), &[]);

        let reader = bcf::Reader::from_path(vcf_file.path()).unwrap();
        let result = prepare_header(reader.header(), bed_file.path(), &aux).unwrap();

        let content = header_to_string(&result);

        let count = content.matches("contig=<ID=chr1>").count();
        assert_eq!(count, 1, "chr1 should appear exactly once, not duplicated");
    }

    #[test]
    fn test_prepare_header_sorts_contigs_globally() {
        let vcf_file = create_test_vcf(&["chr3", "chr5"], &[]);
        let bed_file = create_test_bed(&["chr1", "chr12", "chr13", "chr4"]);
        let aux = make_aux_collector(vcf_file.path(), &[]);

        let reader = bcf::Reader::from_path(vcf_file.path()).unwrap();
        let result = prepare_header(reader.header(), bed_file.path(), &aux).unwrap();
        let content = header_to_string(&result);

        let expected_order = ["chr1", "chr12", "chr13", "chr3", "chr4", "chr5"];
        let positions: Vec<usize> = expected_order
            .iter()
            .map(|chrom| {
                content
                    .find(&format!("contig=<ID={}>", chrom))
                    .unwrap_or_else(|| panic!("{} not found in header", chrom))
            })
            .collect();

        assert!(
            positions.windows(2).all(|w| w[0] < w[1]),
            "contigs not in lexicographic order: expected {:?}, got positions {:?}",
            expected_order,
            positions
        );
    }

    #[test]
    fn test_prepare_header_propagates_aux_info_fields() {
        let vcf_file = create_test_vcf(
            &["chr1"],
            &[r##"##INFO=<ID=COSMIC_ID,Number=1,Type=String,Description="COSMIC ID">"##],
        );
        let bed_file = create_test_bed(&["chr1"]);
        let aux = make_aux_collector(vcf_file.path(), &["COSMIC_ID"]);

        let reader = bcf::Reader::from_path(vcf_file.path()).unwrap();
        let result = prepare_header(reader.header(), bed_file.path(), &aux).unwrap();

        let content = header_to_string(&result);
        assert!(content.contains("ID=COSMIC_ID"));
        assert!(content.contains(&region_id_field()));
    }

    #[test]
    fn test_prepare_header_copies_standard_info_field() {
        let vcf_file = create_test_vcf(
            &["chr1"],
            &[
                r##"##INFO=<ID=HETEROZYGOSITY,Number=A,Type=Float,Description="Prior heterozygosity">"##,
            ],
        );
        let bed_file = create_test_bed(&["chr1"]);
        let aux = make_aux_collector(vcf_file.path(), &[]);

        let reader = bcf::Reader::from_path(vcf_file.path()).unwrap();
        let result = prepare_header(reader.header(), bed_file.path(), &aux).unwrap();

        let content = header_to_string(&result);
        assert!(content.contains(
            r#"ID=HETEROZYGOSITY,Number=A,Type=Float,Description="Prior heterozygosity""#
        ));
    }

    #[test]
    fn test_prepare_header_errors_on_info_field_missing_description() {
        let vcf_file = create_test_vcf(
            &["chr1"],
            &["##INFO=<ID=HETEROZYGOSITY,Number=A,Type=Float>"],
        );
        let bed_file = create_test_bed(&["chr1"]);
        let aux = make_aux_collector(vcf_file.path(), &[]);

        let reader = bcf::Reader::from_path(vcf_file.path()).unwrap();
        let err = prepare_header(reader.header(), bed_file.path(), &aux).unwrap_err();

        assert_eq!(
            err.downcast_ref::<Error>(),
            Some(&Error::VcfHeaderFieldMissing {
                field: "HETEROZYGOSITY".to_string(),
                location: "INFO".to_string(),
            })
        );
    }
}
