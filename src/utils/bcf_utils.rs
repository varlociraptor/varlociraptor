//! bcf_utils.rs
//!
//! Utilities for VCF/BCF file handling and parsing.
//!
//! This module provides:
//! 1. Sample information extraction
//! 2. Record field extraction (chromosome, SVLEN, probabilities, allele frequencies)
//! 3. Allele type classification (indel, symbolic, breakend, reference, spanning deletion)
//! 4. VCF/BCF fields validation
//! 5. VCF file validation
//!
//! NOTE: Helpful General Data Structs or functions (that can be copied) can be found in following PR:
//! 1. SampleInfo Struct for samples and their indexes map: https://github.com/rohan-ibn-tariq/varlociraptor/blob/56278ba1f36f8a89046c3bc9481d502ab7e0b377/src/utils/bcf_utils.rs#L34
//! 2. Function to extract sample names from VCF header: https://github.com/rohan-ibn-tariq/varlociraptor/blob/56278ba1f36f8a89046c3bc9481d502ab7e0b377/src/utils/bcf_utils.rs#L55
//! 3. Function to get combined probability that variant is absent (artifact): https://github.com/rohan-ibn-tariq/varlociraptor/blob/56278ba1f36f8a89046c3bc9481d502ab7e0b377/src/utils/bcf_utils.rs#L161
//! 4. Function to extract per-sample allele frequencies for a specific ALT allele: https://github.com/rohan-ibn-tariq/varlociraptor/blob/56278ba1f36f8a89046c3bc9481d502ab7e0b377/src/utils/bcf_utils.rs#L258
//! 5. Function to get SVLEN from INFO field or calculate dynamically: https://github.com/rohan-ibn-tariq/varlociraptor/blob/56278ba1f36f8a89046c3bc9481d502ab7e0b377/src/utils/bcf_utils.rs#L98
//!

use std::path::Path;

use anyhow::Result;
use bio::stats::{LogProb, Prob};
use log::info;
use rust_htslib::bcf::header::{HeaderView, TagLength, TagType};
use rust_htslib::bcf::{self, record::Numeric, Read};

use crate::constants::{
    MSI_FORMAT_AF_FIELD_LENGTH, MSI_FORMAT_AF_FIELD_TYPE, MSI_INFO_PROB_EVENT_FIELD_LENGTH,
    MSI_INFO_PROB_EVENT_FIELD_TYPE,
};
use crate::errors::Error;
use crate::utils::is_phred_scaled;
use crate::utils::stats::phred_to_prob;

/* ========= BCF Extraction Functions ============= */

/// Get chromosome name from a VCF record
///
/// # Arguments
/// * `record` - VCF record
///
/// # Returns
/// Chromosome name as String
///
/// # Errors
/// Returns error if:
/// - RID is missing in record
/// - Chromosome name resolution fails
///
/// # Example
/// assert_eq!(get_chrom(&record).unwrap(), "chr1");
pub(crate) fn get_chrom(record: &bcf::Record) -> Result<String> {
    let rid = record.rid().ok_or_else(|| Error::VcfRecordChromMissing {
        pos: record.pos() + 1,
    })?;

    let chrom_bytes =
        record
            .header()
            .rid2name(rid)
            .map_err(|_| Error::VcfRecordChromResolveFailed {
                pos: record.pos() + 1,
                rid,
                details: "Failed to resolve chromosome name".to_string(),
            })?;

    let chrom = String::from_utf8_lossy(chrom_bytes).to_string();

    Ok(chrom)
}

/// Get a sample's index in the VCF/BCF header.
///
/// # Returns
/// * `Ok(idx)` if the sample exists
/// * `Err(Error::VcfSamplesNotFound)` if it does not
pub(crate) fn get_sample_index(header: &HeaderView, sample: &str) -> Result<usize> {
    header.sample_id(sample.as_bytes()).ok_or_else(|| {
        Error::VcfSamplesNotFound {
            sample: sample.to_string(),
        }
        .into()
    })
}

/// Combine probabilities for user-specified events into
/// P(at least one specified event) for one ALT allele.
///
/// Note: As per current variant calling implementation in
/// Varlociraptor, these events are assigned to each ALT
/// allele separately (Number=A). See method, fn header(&self)
/// of Caller in src/variants/calling.rs for details. Therefore,
/// if in future implemenation changes this function may need to
/// be updated to handle different Number types (e.g. Number=1 or Number=G).
/// Also it considers float type probabilities, based on the current
/// implementation.
///
/// For each event name, reads `INFO/PROB_{EVENT}` (uppercased),
/// converts to log-space, and sums via `LogProb::ln_sum_exp`.
/// Returns the result as a linear probability.
///
/// Returns `None` if any event field is absent or missing at `alt_idx`.
/// Partial data is treated as unusable rather than silently summing
/// incomplete evidence.
///
/// # Arguments
/// * `record`   - VCF record
/// * `alt_idx`  - ALT allele index (0-based into the ALT array, not alleles array)
/// * `events`   - Event names (e.g. `["somatic", "germline_het"]`)
/// * `is_phred` - Whether `INFO/PROB_*` values are PHRED-scaled
///
/// # Returns
/// * `Ok(Some(p))` - P(at least one event) in [0.0, 1.0]
/// * `Ok(None)`    - Any `PROB_{EVENT}` field absent or missing at `alt_idx`
/// * `Err`         - NaN value in a probability field
pub(crate) fn get_events_probability(
    record: &bcf::Record,
    alt_idx: usize,
    events: &[String],
    is_phred: bool,
) -> Result<Option<f64>> {
    let mut log_probs: Vec<LogProb> = Vec::with_capacity(events.len());

    for event_name in events {
        let field_name = format!("PROB_{}", event_name.to_uppercase());
        let field_bytes = field_name.as_bytes();

        let raw = match record.info(field_bytes).float()? {
            Some(p) if alt_idx < p.len() => p[alt_idx],
            _ => return Ok(None),
        };

        if raw.is_missing() {
            return Ok(None);
        }

        if raw.is_nan() {
            return Err(Error::VcfProbabilityValueInvalid {
                field: field_name,
                value: raw,
                chrom: get_chrom(record)?,
                pos: record.pos() + 1,
            }
            .into());
        }

        let linear = if is_phred {
            phred_to_prob(raw as f64)
        } else {
            raw as f64
        };

        log_probs.push(LogProb::from(Prob(linear.clamp(0.0, 1.0))));
    }

    if log_probs.is_empty() {
        return Ok(None);
    }

    let prob_events = (*Prob::from(LogProb::ln_sum_exp(&log_probs))).clamp(0.0, 1.0);
    Ok(Some(prob_events))
}

/// Extract allele frequency for one sample and one ALT allele.
///
/// Reads `FORMAT/AF` at `sample_idx` and `alt_idx`.
///
/// # Arguments
/// * `record`     - VCF record
/// * `sample_idx` - Sample index in FORMAT columns
/// * `alt_idx`    - ALT allele index (0-based into the ALT array)
///
/// # Returns
/// * `Ok(Some(af))` - AF in [0.0, 1.0]
/// * `Ok(None)`     - Field absent or value missing
/// * `Err`          - AF outside [0.0, 1.0] or NaN
pub(crate) fn get_sample_af(
    record: &bcf::Record,
    sample_idx: usize,
    alt_idx: usize,
) -> Result<Option<f32>> {
    let afs = match record.format(b"AF").float() {
        Ok(a) => a,
        Err(_) => return Ok(None),
    };

    let af = match afs.get(sample_idx).and_then(|s| s.get(alt_idx)).copied() {
        Some(v) => v,
        None => return Ok(None),
    };

    if af.is_missing() {
        return Ok(None);
    }

    if af.is_nan() || !(0.0..=1.0).contains(&af) {
        return Err(Error::VcfAlleleFrequencyInvalid {
            sample: format!("sample_idx={}", sample_idx),
            af,
            chrom: get_chrom(record)?,
            pos: record.pos() + 1,
        }
        .into());
    }

    Ok(Some(af))
}

/* ================================================ */

/* ===== BCF Specification Check Functions ======== */

/// Check if VCF uses PHRED-scaled probabilities from file path.
///
/// Opens the VCF file and checks the description of PROB_{event}s
/// to determine if probabilities are PHRED-scaled or linear.
///
/// # Arguments
/// * `vcf_path` - Path to VCF/BCF file
///
/// # Returns
/// * `Ok(true)` if probabilities are PHRED-scaled
/// * `Ok(false)` if probabilities are linear
/// * `Err` if file cannot be opened
///
/// # Example
/// assert!(is_phred_scaled_from_path(Path::new("variants.vcf")).is_ok());
pub(crate) fn is_phred_scaled_from_path(vcf_path: &Path) -> Result<bool> {
    let vcf = bcf::Reader::from_path(vcf_path)?;
    Ok(is_phred_scaled(&vcf))
}

/* ================================================ */

/* ======== BCF Allele Type Check Functions ======= */

// /// Check if all alleles are SNVs (all same length as ref)
// /// Example: REF=A, ALT=T,G -> all length 1
// NOTE: Not required, can be toggled on if required in other utilities.
// pub fn all_alleles_snv(ref_allele: &[u8], alt_alleles: &[&[u8]]) -> bool {
//     let ref_allele_len = ref_allele.len();
//     alt_alleles.iter().all(|alt| alt.len() == ref_allele_len)
// }

/// Check if allele represents no variant (reference)
///
/// Matches alleles that indicate "no alternative allele":
/// - `.` - Missing/no ALT allele (VCF spec)
/// - `<REF>` - Explicit reference allele (rare)
///
/// # Arguments
/// * `allele` - Allele sequence as byte slice
///
/// # Returns
/// `true` if allele is reference, `false` otherwise
///
/// # Example
/// assert!(is_reference_allele(b"."));
pub(crate) fn is_reference_allele(allele: &[u8]) -> bool {
    allele == b"." || allele == b"<REF>"
}

/// Check if allele is symbolic (starts with <)
///
/// # Arguments
/// * `allele` - Allele sequence as byte slice
///
/// # Returns
/// `true` if allele is symbolic, `false` otherwise
///
/// Examples
///  assert!(is_symbolic(b"<DEL>"));
pub(crate) fn is_symbolic(allele: &[u8]) -> bool {
    allele.len() >= 3 && allele.starts_with(b"<") && allele.ends_with(b">")
}

/// Check if allele is a breakend (contains [ or ])
///
/// # Arguments
/// * `allele` - Allele sequence as byte slice
///
/// # Returns
/// `true` if allele is a breakend, `false` otherwise
///
/// # Examples
/// assert!(is_breakend(b"A[chr2:100["));
pub(crate) fn is_breakend(allele: &[u8]) -> bool {
    allele.iter().any(|&c| c == b'[' || c == b']')
}

/// Check if allele is a spanning deletion (*)
///
/// # Arguments
/// * `allele` - Allele sequence as byte slice
///
/// # Returns
/// `true` if allele is a spanning deletion, `false` otherwise
///
/// # Example
/// assert!(is_spanning_deletion(b"*"));
pub(crate) fn is_spanning_deletion(allele: &[u8]) -> bool {
    allele == b"*"
}

/* ================================================ */

/* ======= BCF Record Query Function tests ======== */

/// Check if a VCF record has a specific INFO string field present.
///
/// Useful for checking optional annotation fields without unwrapping.
///
/// # Arguments
/// * `record` - VCF record to check
/// * `field` - INFO field name as bytes (e.g., b"REGION_ID")
///
/// # Returns
/// `true` if field exists and has at least one value, `false` otherwise
///
/// # Example
/// assert!(record_has_info_string(&record, b"REGION_ID"));
/// assert!(!record_has_info_string(&record, b"MISSING_FIELD"));
pub fn record_has_info_string(record: &bcf::Record, field: &[u8]) -> bool {
    record.info(field).string().ok().flatten().is_some()
}

/// Get all values of an INFO string field from a VCF record.
///
/// Returns all values as a vector. The caller is responsible for
/// deciding how many values to use (e.g., first only, or all).
///
/// # Arguments
/// * `record` - VCF record to query
/// * `field` - INFO field name as bytes (e.g., b"REGION_ID")
///
/// # Returns
/// * `Some(Vec<String>)` - All values if field present
/// * `None` - If field is absent
/// * Panics if a present value is not valid UTF-8 - callers should only use
///   this on fields where UTF-8 validity is guaranteed by construction.
///
/// # Example
/// assert_eq!(get_info_strings(&record, b"REGION_ID").unwrap(), vec!["chr1:100-130".to_string()]);
pub fn get_info_strings(record: &bcf::Record, field: &[u8]) -> Option<Vec<String>> {
    record.info(field).string().ok().flatten().map(|v| {
        v.iter()
            .map(|s| String::from_utf8(s.to_vec()).unwrap())
            .collect()
    })
}

/// Check if a VCF record has a specific INFO flag set.
///
/// INFO flags are boolean - present means true, absent means false.
/// Returns false on any error (field missing, read error).
///
/// # Arguments
/// * `record` - VCF record to check
/// * `field` - INFO flag name as bytes (e.g., b"MSI_DUMMY")
///
/// # Returns
/// `true` if flag is present, `false` if absent or on error
///
/// # Example
/// assert!(record_has_info_flag(&record, b"MSI_DUMMY"));
/// assert!(!record_has_info_flag(&record, b"NONEXISTENT"));
pub fn record_has_info_flag(record: &bcf::Record, field: &[u8]) -> bool {
    record.info(field).flag().unwrap_or(false)
}

/// Read all records from a VCF/BCF file into a vector.
///
/// Convenience function for tests and small files. For large files,
/// use `bcf::Reader` directly.
///
/// # Arguments
/// * `path` - Path to VCF/BCF file
///
/// # Returns
/// * `Ok(Vec<bcf::Record>)` - All records in file order
/// * `Err` - If file cannot be opened or any record fails to parse
///
/// # Example
/// let records = read_bcf_records(Path::new("output.vcf")).unwrap();
/// assert_eq!(records.len(), 3);
pub fn read_bcf_records(path: &Path) -> Result<Vec<bcf::Record>> {
    let mut reader = bcf::Reader::from_path(path).map_err(|_| Error::VcfFileInvalid {
        path: path.to_path_buf(),
    })?;
    reader
        .records()
        .map(|r| {
            r.map_err(|e| {
                Error::VcfRecordReadFailed {
                    details: e.to_string(),
                }
                .into()
            })
        })
        .collect()
}

/* ================================================ */

/* ======= BCF Record INFO Copy Function ========= */

/// Copy specified INFO fields from source to destination BCF record.
///
/// Silently skips fields not declared in the source header. A field that
/// is declared but fails to read (e.g. a value not matching its declared
/// type or count) returns an `Err` rather than being silently dropped.
/// Handles all VCF INFO field types: Integer, Float, String, Flag.
///
/// Designed for use in preprocessing pipelines where a fresh output
/// record needs selected INFO fields carried over from the input.
///
/// # Arguments
/// * `source` - Source BCF record to copy fields from
/// * `dest` - Destination BCF record to write fields to
/// * `fields` - Field names to copy as string slices
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if reading source or writing destination fails
///
/// # Example
/// assert!(copy_info_fields(&source_record, &mut dest_record, &["SVLEN", "SVTYPE"]).is_ok());
pub(crate) fn copy_info_fields(
    source: &bcf::Record,
    dest: &mut bcf::Record,
    fields: &[&str],
) -> Result<()> {
    let header = source.header();

    for field in fields {
        let field_bytes = field.as_bytes();

        match header.info_type(field_bytes) {
            Ok((TagType::Integer, _)) => {
                if let Some(values) = source.info(field_bytes).integer()? {
                    dest.push_info_integer(field_bytes, &values)?;
                }
            }
            Ok((TagType::Float, _)) => {
                if let Some(values) = source.info(field_bytes).float()? {
                    dest.push_info_float(field_bytes, &values)?;
                }
            }
            Ok((TagType::String, _)) => {
                if let Some(values) = source.info(field_bytes).string()? {
                    dest.push_info_string(field_bytes, &values)?;
                }
            }
            Ok((TagType::Flag, _)) => {
                if source.info(field_bytes).flag()? {
                    dest.push_info_flag(field_bytes)?;
                }
            }
            Err(_) => {
                // Field absent in source - silently skip
            }
        }
    }
    Ok(())
}

/* ================================================ */

/* ========= BCF Validation Functions ============= */

/// Validate a single VCF header field has correct type and number.
///
/// # Arguments
/// * `field_result` - Result from `header.info_type()` or `header.format_type()`
/// * `location` - "INFO" or "FORMAT" (for error messages)
/// * `field_name` - Name of the field being validated
/// * `expected_type` - Expected TagType (e.g., `TagType::Float`)
/// * `expected_length` - Expected TagLength (e.g., `TagLength::AltAlleles`)
///
/// # Returns
/// * `Ok(())` if field exists with correct type and number
/// * `Err` if field missing or has wrong type/number
fn validate_vcf_header_field(
    field_result: std::result::Result<(TagType, TagLength), rust_htslib::errors::Error>,
    location: &str,
    field_name: &str,
    expected_type: TagType,
    expected_length: TagLength,
) -> Result<()> {
    match field_result {
        Ok((actual_type, actual_length)) => {
            if actual_type != expected_type || actual_length != expected_length {
                Err(Error::VcfHeaderFieldTypeInvalid {
                    location: location.to_string(),
                    field: field_name.to_string(),
                    expected: format!("Type={:?}, Number={:?}", expected_type, expected_length),
                    found: format!("Type={:?}, Number={:?}", actual_type, actual_length),
                }
                .into())
            } else {
                Ok(())
            }
        }
        Err(_) => Err(Error::VcfHeaderFieldMissing {
            field: field_name.to_string(),
            location: location.to_string(),
        }
        .into()),
    }
}

/// Validate that INFO fields exist in VCF header.
///
/// Checks VCF header for presence of each INFO field name. It does not
/// check field type, since fields may be any INFO type and intent here
/// is just to check the presence of the field.
///
/// # Arguments
/// * `header` - VCF header to check for INFO field declarations
/// * `fields` - Slice of INFO field names to validate (e.g., ["HETEROZYGOSITY"])
///
/// # Returns
/// * `Ok(())` if all fields exist
/// * `Err` if any field is missing
///
/// # Errors
/// Returns error if any requested INFO field is not declared in the VCF header.
/// Error message includes comma-separated list of missing fields.
///
/// # Example
/// assert!(validate_info_fields_exist(&header, &vec!["HETEROZYGOSITY".to_string()]).is_ok());
pub(crate) fn validate_info_fields_exist(header: &HeaderView, fields: &[String]) -> Result<()> {
    let missing: Vec<String> = fields
        .iter()
        .filter(|f| header.info_type(f.as_bytes()).is_err())
        .cloned()
        .collect();

    if !missing.is_empty() {
        return Err(Error::VcfInfoFieldsMissing {
            fields: missing.join(", "),
        }
        .into());
    }
    info!("  - INFO fields validated: {:?}", fields);
    Ok(())
}

/// Validate that required samples exist in VCF header.
///
/// Checks VCF header for presence of all requested sample names.
///
/// # Arguments
/// * `header` - VCF header to validate
/// * `required_samples` - Slice of sample names to validate
///
/// # Returns
/// * `Ok(())` if all samples exist
/// * `Err` if any sample is missing
///
/// # Errors
/// Returns error if any sample name is not found in VCF header.
/// Error message includes comma-separated list of missing samples.
///
/// # Example
/// assert!(validate_samples_exist(&header, &vec!["tumor".to_string()]).is_ok());
pub(crate) fn validate_samples_exist(
    header: &HeaderView,
    required_samples: &[String],
) -> Result<()> {
    let mut missing = Vec::new();

    for sample_name in required_samples {
        if header.sample_id(sample_name.as_bytes()).is_none() {
            missing.push(sample_name.clone());
        }
    }

    if !missing.is_empty() {
        return Err(Error::VcfSamplesNotFound {
            sample: missing.join(", "),
        }
        .into());
    }

    info!("  - Samples validated: {:?}", required_samples);

    Ok(())
}

/// Validate that required events exist in VCF header, with correct type.
///
/// Checks VCF header for presence of INFO/PROB_{EVENT} fields
/// for each requested event name, and that each present field has the expected
/// shape emitted by varlociraptor's calling model: `Type=Float, Number=A`
/// (one probability per ALT allele). Event names are automatically
/// converted to uppercase and prefixed with "PROB_".
///
/// # Arguments
/// * `header` - VCF header to check for event fields
/// * `event_names` - Slice of event names to validate (e.g., ["somatic_tumor"])
///
/// # Returns
/// * `Ok(())` if all event fields exist with the correct type/shape
/// * `Err` if any event field is missing, or exists with the wrong type/shape
///
/// # Errors
/// * `VcfEventsMissing` if any INFO/PROB_{EVENT} field is not found. Error
///   message includes comma-separated list of missing events.
/// * `VcfHeaderFieldTypeInvalid` if a PROB_{EVENT} field is found but isn't
///   `Type=Float, Number=A`.
///
/// # Event Field Mapping
/// Event name          -> INFO field checked
/// "somatic_tumor"     -> INFO/PROB_SOMATIC_TUMOR
/// "germline_normal"   -> INFO/PROB_GERMLINE_NORMAL
///
/// # Example
/// assert!(validate_events_exist(&header, &vec!["somatic_tumor".to_string()]).is_ok());
pub(crate) fn validate_events_exist(header: &HeaderView, event_names: &[String]) -> Result<()> {
    let mut missing_events = Vec::new();

    for event_name in event_names {
        let field_name = format!("PROB_{}", event_name.to_uppercase());

        match header.info_type(field_name.as_bytes()) {
            Err(_) => missing_events.push(event_name.clone()),
            Ok(_) => {
                validate_vcf_header_field(
                    header.info_type(field_name.as_bytes()),
                    "INFO",
                    &field_name,
                    MSI_INFO_PROB_EVENT_FIELD_TYPE,
                    MSI_INFO_PROB_EVENT_FIELD_LENGTH,
                )?;
            }
        }
    }

    if !missing_events.is_empty() {
        return Err(Error::VcfEventsMissing {
            events: missing_events.join(", "),
        }
        .into());
    }

    info!("  - Events validated: {:?}", event_names);
    Ok(())
}

/// Validate VCF header contains required fields for MSI analysis.
///
/// Validates that mandatory INFO or FORMAT fields exist with correct types:
/// - FORMAT:AF (Type=Float, Number=A)
///
/// # Arguments
/// * `header` - VCF header to validate
///
/// # Returns
/// * `Ok(())` if all required fields present with correct types
/// * `Err` if any field is missing or has incorrect type
pub(crate) fn validate_required_vcf_fields_msi(header: &HeaderView) -> Result<()> {
    validate_vcf_header_field(
        header.format_type(b"AF"),
        "FORMAT",
        "AF",
        MSI_FORMAT_AF_FIELD_TYPE,
        MSI_FORMAT_AF_FIELD_LENGTH,
    )?;

    Ok(())
}

/// Validate VCF file has at least one variant.
///
/// Simple sanity check - ensures file is not empty and logs first variant position.
///
/// # Arguments
/// * `vcf` - VCF reader (must be at start of file)
///
/// # Returns
/// * `Ok(())` if file has at least one variant
/// * `Err` if file is empty or read fails
///
/// # Note
/// Advances the reader past the first record - reopen the file before
/// further iteration.
pub(crate) fn validate_vcf_file(vcf: &mut bcf::Reader) -> Result<()> {
    match vcf.records().next() {
        None => Err(Error::VcfFileEmpty.into()),
        Some(Err(e)) => Err(Error::VcfRecordReadFailed {
            details: e.to_string(),
        }
        .into()),
        Some(Ok(record)) => {
            let chrom = get_chrom(&record)?;
            let pos = record.pos();
            info!("  - First variant: {}:{}", chrom, pos + 1);
            info!("  - VCF file validated successfully");
            Ok(())
        }
    }
}

/* ================================================ */

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    use crate::constants::test_constants::{TEST_EPSILON, TEST_EPSILON_F32, TEST_EPSILON_LOOSE};
    use crate::utils::genomics::is_indel;

    /// Configuration for test VCF creation
    pub(crate) struct TestVcfConfig<'a> {
        pub ref_allele: &'a [u8],
        pub alt_alleles: Vec<&'a [u8]>,
        pub af_values: Option<Vec<f32>>,
        pub prob_absent: Option<Vec<f32>>,
        pub prob_artifact: Option<Vec<f32>>,
        pub prob_somatic: Option<Vec<f32>>,
        pub prob_high_vaf: Option<Vec<f32>>,
        pub num_samples: usize,
        pub use_phred: bool,
        /// Extra INFO fields: (id, number, type, description, value)
        /// e.g. (b"COSMIC_ID", b"1", b"String", b"COSMIC ID", b"COSM123")
        pub extra_info_fields: Vec<(&'a [u8], &'a [u8], &'a [u8], &'a [u8], &'a [u8])>,
        pub write_prob_somatic: bool,
        pub write_prob_high_vaf: bool,
    }

    impl<'a> Default for TestVcfConfig<'a> {
        fn default() -> Self {
            Self {
                ref_allele: b"A",
                alt_alleles: vec![b"AT"],
                af_values: None,
                prob_absent: None,
                prob_artifact: None,
                prob_somatic: None,
                prob_high_vaf: None,
                num_samples: 2,
                use_phred: false,
                extra_info_fields: vec![],
                write_prob_somatic: true,
                write_prob_high_vaf: true,
            }
        }
    }

    pub(crate) fn create_test_vcf(config: TestVcfConfig) -> (NamedTempFile, Vec<String>) {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();

        let mut header = rust_htslib::bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");
        header.push_record(br##"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">"##);

        // Extra INFO field declarations
        for (id, number, type_, desc, _) in &config.extra_info_fields {
            header.push_record(
                format!(
                    "##INFO=<ID={},Number={},Type={},Description=\"{}\">",
                    String::from_utf8_lossy(id),
                    String::from_utf8_lossy(number),
                    String::from_utf8_lossy(type_),
                    String::from_utf8_lossy(desc),
                )
                .as_bytes(),
            );
        }

        // Conditionally set PROB_ABSENT, PROB_ARTIFACT and PROB_SOMATIC headers based on use_phred
        if config.use_phred {
            header.push_record(br##"##INFO=<ID=PROB_ABSENT,Number=A,Type=Float,Description="Probability absent (PHRED)">"##);
            header.push_record(br##"##INFO=<ID=PROB_ARTIFACT,Number=A,Type=Float,Description="Probability artifact (PHRED)">"##);
            header.push_record(br##"##INFO=<ID=PROB_SOMATIC,Number=A,Type=Float,Description="Probability somatic (PHRED)">"##);
            header.push_record(br##"##INFO=<ID=PROB_HIGH_VAF,Number=A,Type=Float,Description="Probability high VAF (PHRED)">"##);
        } else {
            header.push_record(br##"##INFO=<ID=PROB_ABSENT,Number=A,Type=Float,Description="Probability absent (linear)">"##);
            header.push_record(br##"##INFO=<ID=PROB_ARTIFACT,Number=A,Type=Float,Description="Probability artifact (linear)">"##);
            header.push_record(br##"##INFO=<ID=PROB_SOMATIC,Number=A,Type=Float,Description="Probability somatic (linear)">"##);
            header.push_record(br##"##INFO=<ID=PROB_HIGH_VAF,Number=A,Type=Float,Description="Probability high VAF (linear)">"##);
        }

        header.push_record(br##"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"##);
        header.push_record(
            br##"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">"##,
        );

        let mut sample_names = Vec::new();
        for i in 0..config.num_samples {
            let name = format!("sample{}", i + 1);
            header.push_sample(name.as_bytes());
            sample_names.push(name);
        }

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

        // Build alleles: [REF, ALT1, ALT2, ...]
        let mut alleles_vec = vec![config.ref_allele];
        alleles_vec.extend(config.alt_alleles.iter());
        rec.set_alleles(&alleles_vec).unwrap();

        // SVLEN for each ALT
        let svlen_values: Vec<i32> = config
            .alt_alleles
            .iter()
            .map(|alt| (alt.len() as i32) - (config.ref_allele.len() as i32))
            .collect();
        rec.push_info_integer(b"SVLEN", &svlen_values).unwrap();

        // Extra INFO fields with specified values
        for (id, _, type_, _, value) in &config.extra_info_fields {
            match *type_ {
                b"Integer" => {
                    let v: i32 = std::str::from_utf8(value).unwrap().parse().unwrap();
                    rec.push_info_integer(id, &[v]).unwrap();
                }
                b"Float" => {
                    let v: f32 = std::str::from_utf8(value).unwrap().parse().unwrap();
                    rec.push_info_float(id, &[v]).unwrap();
                }
                b"Flag" => {
                    rec.push_info_flag(id).unwrap();
                }
                b"String" => {
                    rec.push_info_string(id, &[*value]).unwrap();
                }
                _ => {
                    panic!(
                        "Unsupported INFO field type in test config: {}",
                        String::from_utf8_lossy(type_)
                    );
                }
            }
        }

        // PROB_ABSENT
        let prob_absent = match config.prob_absent {
            Some(pa) => pa,
            None => {
                if config.use_phred {
                    vec![20.0; config.alt_alleles.len()]
                } else {
                    vec![0.005; config.alt_alleles.len()]
                }
            }
        };
        rec.push_info_float(b"PROB_ABSENT", &prob_absent).unwrap();

        // PROB_ARTIFACT
        let prob_artifact = match config.prob_artifact {
            Some(pa) => pa,
            None => {
                if config.use_phred {
                    vec![10.0; config.alt_alleles.len()]
                } else {
                    vec![0.005; config.alt_alleles.len()]
                }
            }
        };
        rec.push_info_float(b"PROB_ARTIFACT", &prob_artifact)
            .unwrap();

        // PROB_SOMATIC
        if config.write_prob_somatic {
            let prob_somatic = match config.prob_somatic {
                Some(ps) => ps,
                None => {
                    if config.use_phred {
                        vec![30.0; config.alt_alleles.len()]
                    } else {
                        vec![0.001; config.alt_alleles.len()]
                    }
                }
            };
            rec.push_info_float(b"PROB_SOMATIC", &prob_somatic).unwrap();
        }

        // PROB_HIGH_VAF
        if config.write_prob_high_vaf {
            let prob_high_vaf = match config.prob_high_vaf {
                Some(phv) => phv,
                None => {
                    if config.use_phred {
                        vec![40.0; config.alt_alleles.len()]
                    } else {
                        vec![0.0001; config.alt_alleles.len()]
                    }
                }
            };
            rec.push_info_float(b"PROB_HIGH_VAF", &prob_high_vaf)
                .unwrap();
        }

        // GT
        let mut genotypes = Vec::new();
        for _ in 0..config.num_samples {
            genotypes.push(i32::from(bcf::record::GenotypeAllele::Unphased(0)));
            genotypes.push(i32::from(bcf::record::GenotypeAllele::Unphased(1)));
        }
        rec.push_format_integer(b"GT", &genotypes).unwrap();

        // AF values
        if let Some(afs) = config.af_values {
            rec.push_format_float(b"AF", &afs).unwrap();
        }

        wtr.write(&rec).unwrap();

        (tmp, sample_names)
    }

    fn create_header_view(header_lines: &[&[u8]]) -> HeaderView {
        let mut header = bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");

        for line in header_lines {
            header.push_record(line);
        }

        let tmp = NamedTempFile::new().unwrap();
        let writer = bcf::Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();
        drop(writer);

        let reader = bcf::Reader::from_path(tmp.path()).unwrap();
        reader.header().clone()
    }

    /// Read the first record from a VCF/BCF file (without header).
    ///
    /// # Returns
    /// bcf::Record for testing
    pub fn read_first_record(path: &Path) -> bcf::Record {
        let mut reader = bcf::Reader::from_path(path).unwrap();
        reader.records().next().unwrap().unwrap()
    }

    /// Create a test VCF record with one or more ALT alleles.
    ///
    /// # Arguments
    /// * `writer`      - BCF writer
    /// * `rid`         - Reference ID (chromosome index in header)
    /// * `pos`         - Position (0-based)
    /// * `ref_allele`  - Reference allele bytes (e.g. `b"ACAG"`)
    /// * `alt_alleles` - Slice of ALT allele byte slices (e.g. `&[b"ACAGCAG", b"A"]`)
    ///
    /// # Returns
    /// Configured BCF record ready for testing.
    pub(crate) fn create_test_record_multi_alt(
        writer: &bcf::Writer,
        rid: u32,
        pos: i64,
        ref_allele: &[u8],
        alt_alleles: &[&[u8]],
    ) -> bcf::Record {
        let mut record = writer.empty_record();
        record.set_rid(Some(rid));
        record.set_pos(pos);
        let mut alleles: Vec<&[u8]> = vec![ref_allele];
        alleles.extend_from_slice(alt_alleles);
        record.set_alleles(&alleles).unwrap();
        record
    }

    /// Create a test VCF record with a single ALT allele.
    ///
    /// # Arguments
    /// * `writer` - VCF writer
    /// * `rid` - Reference ID (chromosome index in header)
    /// * `pos` - Position (0-based)
    /// * `ref_allele` - Reference allele (e.g., b"A")
    /// * `alt_allele` - Alternate allele (e.g., b"T")
    ///
    /// # Returns
    /// Configured VCF record ready for testing
    pub fn create_test_record(
        writer: &bcf::Writer,
        rid: u32,
        pos: i64,
        ref_allele: &[u8],
        alt_allele: &[u8],
    ) -> bcf::Record {
        create_test_record_multi_alt(writer, rid, pos, ref_allele, &[alt_allele])
    }

    /// Minimal multi-record VCF builder for position/overlap-driven tests.
    ///
    /// # Note:
    /// Unlike [`create_test_vcf`], which builds exactly one feature-rich
    /// record (samples, AF, PROB_*, extra INFO) from a [`TestVcfConfig`],
    /// this builds a header + N bare records (chrom/pos/REF/ALT only, no
    /// samples/FORMAT). Use this when a test needs several records at
    /// once and doesn't care about sample/format.
    ///
    /// # Arguments
    /// * `header_lines` - Extra `##...` header lines (e.g. `##contig=<...>`)
    ///   to push after `##fileformat=VCFv4.2`. Caller is responsible for
    ///   declaring any contigs referenced by `records`.
    /// * `records` - `(rid, pos, ref_allele, alt_alleles)` tuples, written
    ///   in the given order. `rid` indexes into the contigs declared via
    ///   `header_lines`, in declaration order.
    ///
    /// # Returns
    /// A `NamedTempFile` containing the resulting VCF.
    pub(crate) fn create_minimal_vcf(
        header_lines: &[&[u8]],
        records: &[(u32, i64, &[u8], &[&[u8]])],
    ) -> NamedTempFile {
        let tmp = NamedTempFile::new().unwrap();
        let mut header = bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        for line in header_lines {
            header.push_record(line);
        }
        let mut wtr = bcf::Writer::from_path(tmp.path(), &header, true, bcf::Format::Vcf).unwrap();
        for (rid, pos, ref_allele, alts) in records {
            let rec = create_test_record_multi_alt(&wtr, *rid, *pos, ref_allele, alts);
            wtr.write(&rec).unwrap();
        }
        tmp
    }

    /* ==== BCF Extraction Function(s) tests ========= */

    #[test]
    fn test_get_chrom() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig::default());

        let record = read_first_record(tmp_vcf.path());
        let chrom = get_chrom(&record).unwrap();
        assert_eq!(chrom, "chr1");
    }

    #[test]
    fn test_get_sample_index_found() {
        let (tmp_vcf, sample_names) = create_test_vcf(TestVcfConfig {
            num_samples: 2,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let result = get_sample_index(header, &sample_names[1]);
        assert_eq!(result.unwrap(), 1);
        let result = get_sample_index(header, &sample_names[0]);
        assert_eq!(result.unwrap(), 0);
    }

    #[test]
    fn test_get_sample_index_missing() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let result = get_sample_index(header, "nonexistent_sample");

        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("nonexistent_sample"));
    }

    #[test]
    fn test_get_events_probability_single_event() {
        // PROB_SOMATIC = 0.9: P(events) = 0.9
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            prob_somatic: Some(vec![0.9]),
            num_samples: 1,
            ..Default::default()
        });
        let record = read_first_record(tmp_vcf.path());

        let p = get_events_probability(&record, 0, &["somatic".to_string()], false)
            .unwrap()
            .unwrap();
        assert!((p - 0.9).abs() < TEST_EPSILON_LOOSE);
    }

    #[test]
    fn test_get_events_probability_two_events_log_sum() {
        // PROB_SOMATIC=0.3, PROB_HIGH_VAF=0.4
        // ln_sum_exp(ln(0.3), ln(0.4)) = ln(0.7) : P(events) ≈ 0.7
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            prob_somatic: Some(vec![0.3]),
            prob_high_vaf: Some(vec![0.4]),
            num_samples: 1,
            ..Default::default()
        });
        let record = read_first_record(tmp_vcf.path());

        let p = get_events_probability(
            &record,
            0,
            &["somatic".to_string(), "high_vaf".to_string()],
            false,
        )
        .unwrap()
        .unwrap();
        assert!((p - 0.7).abs() < TEST_EPSILON_LOOSE);
    }

    #[test]
    fn test_get_events_probability_phred_scaled() {
        // PHRED 10: linear 0.1
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            use_phred: true,
            prob_somatic: Some(vec![10.0]),
            num_samples: 1,
            ..Default::default()
        });
        let record = read_first_record(tmp_vcf.path());

        let p = get_events_probability(&record, 0, &["somatic".to_string()], true)
            .unwrap()
            .unwrap();
        assert!((p - 0.1).abs() < TEST_EPSILON_LOOSE);
    }

    #[test]
    fn test_get_events_probability_at_boundary_one() {
        // P=1.0 should not error or exceed 1.0 after log-space
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            prob_somatic: Some(vec![1.0]),
            num_samples: 1,
            ..Default::default()
        });
        let record = read_first_record(tmp_vcf.path());

        let p = get_events_probability(&record, 0, &["somatic".to_string()], false)
            .unwrap()
            .unwrap();
        assert!((p - 1.0).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_get_events_probability_field_in_header_not_on_record() {
        // write_prob_somatic: false, PROB_SOMATIC in header but no value on record
        // function returns Ok(None)
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            write_prob_somatic: false,
            num_samples: 1,
            ..Default::default()
        });
        let record = read_first_record(tmp_vcf.path());

        assert!(
            get_events_probability(&record, 0, &["somatic".to_string()], false)
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn test_get_sample_af_valid_two_samples() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            af_values: Some(vec![0.75, 0.25]),
            num_samples: 2,
            ..Default::default()
        });
        let record = read_first_record(tmp_vcf.path());

        assert!((get_sample_af(&record, 0, 0).unwrap().unwrap() - 0.75).abs() < TEST_EPSILON_F32);
        assert!((get_sample_af(&record, 1, 0).unwrap().unwrap() - 0.25).abs() < TEST_EPSILON_F32);
    }

    #[test]
    fn test_get_sample_af_missing_returns_none() {
        use rust_htslib::bcf::record::Numeric;
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            af_values: Some(vec![f32::missing()]),
            num_samples: 1,
            ..Default::default()
        });
        let record = read_first_record(tmp_vcf.path());

        assert!(get_sample_af(&record, 0, 0).unwrap().is_none());
    }

    #[test]
    fn test_get_sample_af_out_of_range_errors() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            af_values: Some(vec![1.5]),
            num_samples: 1,
            ..Default::default()
        });
        let record = read_first_record(tmp_vcf.path());

        assert!(get_sample_af(&record, 0, 0).is_err());
    }

    #[test]
    fn test_get_sample_af_no_field_returns_none() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            af_values: None,
            num_samples: 1,
            ..Default::default()
        });
        let record = read_first_record(tmp_vcf.path());

        assert!(get_sample_af(&record, 0, 0).unwrap().is_none());
    }

    #[test]
    fn test_get_sample_af_multi_alt() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"A",
            alt_alleles: vec![b"AT", b"ATT"],
            af_values: Some(vec![0.6, 0.3]), // sample1: ALT1=0.6, ALT2=0.3
            num_samples: 1,
            ..Default::default()
        });
        let record = read_first_record(tmp_vcf.path());

        let af0 = get_sample_af(&record, 0, 0).unwrap().unwrap();
        let af1 = get_sample_af(&record, 0, 1).unwrap().unwrap();
        assert!((af0 - 0.6).abs() < TEST_EPSILON_F32);
        assert!((af1 - 0.3).abs() < TEST_EPSILON_F32);
    }

    #[test]
    fn test_get_sample_af_sample_idx_out_of_bounds() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            af_values: Some(vec![0.5]),
            num_samples: 1,
            ..Default::default()
        });
        let record = read_first_record(tmp_vcf.path());

        // sample_idx=5 doesn't exist — should return None not panic
        assert!(get_sample_af(&record, 5, 0).unwrap().is_none());
    }

    /* ====== BCF Specification check tests ===== ==== */

    #[test]
    fn test_is_phred_scaled_from_path_linear() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            use_phred: false,
            num_samples: 1,
            ..Default::default()
        });

        let is_phred = is_phred_scaled_from_path(tmp_vcf.path()).unwrap();
        assert!(!is_phred, "Should detect linear probabilities");
    }

    #[test]
    fn test_is_phred_scaled_from_path_phred() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            use_phred: true,
            num_samples: 1,
            ..Default::default()
        });

        let is_phred = is_phred_scaled_from_path(tmp_vcf.path()).unwrap();
        assert!(is_phred, "Should detect PHRED-scaled probabilities");
    }

    #[test]
    fn test_is_phred_scaled_from_path_invalid_file() {
        let result = is_phred_scaled_from_path(Path::new("/nonexistent/file.vcf"));
        assert!(result.is_err(), "Should fail for nonexistent file");
    }

    /* ====== BCF ALLELE Type tests ================== */

    #[test]
    fn test_allele_type_checks() {
        // Reference allele
        assert!(is_reference_allele(b"."));
        assert!(is_reference_allele(b"<REF>"));

        // Indel
        assert!(is_indel(b"AC", b"ACG")); // insertion
        assert!(is_indel(b"ACG", b"AC")); // deletion
        assert!(!is_indel(b"AC", b"AG"));

        // Symbolic
        assert!(is_symbolic(b"<DEL>"));
        assert!(is_symbolic(b"<DUP:TANDEM>"));
        assert!(!is_symbolic(b"ACG"));

        // Breakend
        assert!(is_breakend(b"A[chr2:100["));
        assert!(!is_breakend(b"ACG"));

        // Spanning deletion
        assert!(is_spanning_deletion(b"*"));
        assert!(!is_spanning_deletion(b"AC"));
    }

    /* ======= BCF Record Query Function tests ======== */

    /// Create a minimal VCF writer for testing INFO field operations.
    ///
    /// Creates a temporary VCF file with a minimal header containing
    /// only `fileformat` and `chr1` contig, plus any additional header
    /// records provided by the caller (e.g., INFO field definitions).
    ///
    /// The returned `NamedTempFile` must be kept alive for the duration
    /// of the test - dropping it deletes the file.
    ///
    /// # Arguments
    /// * `info_records` - Additional header lines to include
    ///   (e.g., `br##"##INFO=<ID=REGION_ID,Number=1,Type=String,Description="Region">"##`)
    ///
    /// # Returns
    /// Tuple of `(tmp_file, writer)` - caller writes records via writer,
    /// then drops writer.
    ///
    /// # Example
    /// let (tmp, mut writer) = create_info_test_vcf(&[
    ///     br##"##INFO=<ID=REGION_ID,Number=1,Type=String,Description="Region">"##
    /// ]);
    fn create_info_test_vcf(info_records: &[&[u8]]) -> (NamedTempFile, bcf::Writer) {
        let tmp = NamedTempFile::new().unwrap();
        let mut header = bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");
        for record in info_records {
            header.push_record(record);
        }
        let writer = bcf::Writer::from_path(tmp.path(), &header, true, bcf::Format::Vcf).unwrap();
        (tmp, writer)
    }

    #[test]
    fn test_record_has_info_string_present_and_absent() {
        // Present: REGION_ID set - true
        // Absent:  NONEXISTENT   - false
        let (tmp, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=REGION_ID,Number=1,Type=String,Description="Region">"##,
        ]);
        let mut record = create_test_record(&writer, 0, 100, b"A", b"AT");
        record
            .push_info_string(b"REGION_ID", &[b"chr1:100-130"])
            .unwrap();
        writer.write(&record).unwrap();
        drop(writer);

        let record = read_first_record(tmp.path());
        assert!(record_has_info_string(&record, b"REGION_ID"));
        assert!(!record_has_info_string(&record, b"NONEXISTENT"));
    }

    #[test]
    fn test_get_info_strings_single_value() {
        // Number=1 field - Vec with one element
        let (tmp, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=REGION_ID,Number=1,Type=String,Description="Region">"##,
        ]);
        let mut record = create_test_record(&writer, 0, 100, b"A", b"AT");
        record
            .push_info_string(b"REGION_ID", &[b"chr1:100-130"])
            .unwrap();
        writer.write(&record).unwrap();
        drop(writer);

        let record = read_first_record(tmp.path());
        assert_eq!(
            get_info_strings(&record, b"REGION_ID"),
            Some(vec!["chr1:100-130".to_string()])
        );
    }

    #[test]
    fn test_get_info_strings_multiple_values() {
        // Number=. field - Vec with all elements, caller decides how many to use
        let (tmp, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=TAGS,Number=.,Type=String,Description="Tags">"##,
        ]);
        let mut record = create_test_record(&writer, 0, 100, b"A", b"AT");
        record
            .push_info_string(b"TAGS", &[b"val1", b"val2"])
            .unwrap();
        writer.write(&record).unwrap();
        drop(writer);

        let record = read_first_record(tmp.path());
        assert_eq!(
            get_info_strings(&record, b"TAGS"),
            Some(vec!["val1".to_string(), "val2".to_string()])
        );
    }

    #[test]
    fn test_get_info_strings_absent_field() {
        // Absent field - None
        let (tmp, mut writer) = create_info_test_vcf(&[]);
        let record = create_test_record(&writer, 0, 100, b"A", b"AT");
        writer.write(&record).unwrap();
        drop(writer);

        let record = read_first_record(tmp.path());
        assert_eq!(get_info_strings(&record, b"NONEXISTENT"), None);
    }

    #[test]
    fn test_record_has_info_flag_present_and_absent() {
        // Present: MSI_DUMMY set - true
        // Absent:  NONEXISTENT   - false
        let (tmp, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=MSI_DUMMY,Number=0,Type=Flag,Description="Dummy">"##,
        ]);
        let mut record = create_test_record(&writer, 0, 100, b"A", b"AT");
        record.push_info_flag(b"MSI_DUMMY").unwrap();
        writer.write(&record).unwrap();
        drop(writer);

        let record = read_first_record(tmp.path());
        assert!(record_has_info_flag(&record, b"MSI_DUMMY"));
        assert!(!record_has_info_flag(&record, b"NONEXISTENT"));
    }

    #[test]
    fn test_record_has_info_flag_not_set() {
        // Flag declared in header but not set on record - false
        let (tmp, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=MSI_DUMMY,Number=0,Type=Flag,Description="Dummy">"##,
        ]);
        let record = create_test_record(&writer, 0, 100, b"A", b"AT");
        writer.write(&record).unwrap();
        drop(writer);

        let record = read_first_record(tmp.path());
        assert!(!record_has_info_flag(&record, b"MSI_DUMMY"));
    }

    #[test]
    fn test_read_bcf_records_multiple() {
        // Three records at different positions - all returned in order
        let (tmp, mut writer) = create_info_test_vcf(&[]);
        for pos in [100, 200, 300] {
            let record = create_test_record(&writer, 0, pos, b"A", b"AT");
            writer.write(&record).unwrap();
        }
        drop(writer);

        let records = read_bcf_records(tmp.path()).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].pos(), 100);
        assert_eq!(records[1].pos(), 200);
        assert_eq!(records[2].pos(), 300);
    }

    #[test]
    fn test_read_bcf_records_empty_file() {
        // No records written - empty Vec, not an error
        let (tmp, writer) = create_info_test_vcf(&[]);
        drop(writer);

        let records = read_bcf_records(tmp.path()).unwrap();
        assert_eq!(records.len(), 0);
    }

    #[test]
    fn test_read_bcf_records_nonexistent_file() {
        // Nonexistent path - VcfFileInvalid error
        let result = read_bcf_records(Path::new("/nonexistent/file.vcf"));
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("invalid"));
    }

    /* ======= BCF Record INFO Copy Functions ========= */

    #[test]
    fn test_copy_info_fields_integer() {
        // Source: record with SVLEN=3
        let (tmp_src, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">"##,
        ]);
        let mut record = create_test_record(&writer, 0, 100, b"A", b"AT");
        record.push_info_integer(b"SVLEN", &[3]).unwrap();
        writer.write(&record).unwrap();
        drop(writer);

        // Dest: fresh record, copy, assert
        let (tmp_dst, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">"##,
        ]);
        let src = read_first_record(tmp_src.path());
        let mut dst = create_test_record(&writer, 0, 100, b"A", b"AT");
        copy_info_fields(&src, &mut dst, &["SVLEN"]).unwrap();
        writer.write(&dst).unwrap();
        drop(writer);

        let result = read_first_record(tmp_dst.path());
        assert_eq!(result.info(b"SVLEN").integer().unwrap().unwrap()[0], 3);
    }

    #[test]
    fn test_copy_info_fields_float() {
        let (tmp_src, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=SCORE,Number=A,Type=Float,Description="Score">"##,
        ]);
        let mut record = create_test_record(&writer, 0, 100, b"A", b"AT");
        record.push_info_float(b"SCORE", &[0.42]).unwrap();
        writer.write(&record).unwrap();
        drop(writer);

        let (tmp_dst, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=SCORE,Number=A,Type=Float,Description="Score">"##,
        ]);
        let src = read_first_record(tmp_src.path());
        let mut dst = create_test_record(&writer, 0, 100, b"A", b"AT");
        copy_info_fields(&src, &mut dst, &["SCORE"]).unwrap();
        writer.write(&dst).unwrap();
        drop(writer);

        let result = read_first_record(tmp_dst.path());
        assert!(
            (result.info(b"SCORE").float().unwrap().unwrap()[0] - 0.42).abs() < TEST_EPSILON as f32
        );
    }

    #[test]
    fn test_copy_info_fields_string() {
        let (tmp_src, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">"##,
        ]);
        let mut record = create_test_record(&writer, 0, 100, b"A", b"AT");
        record.push_info_string(b"SVTYPE", &[b"INS"]).unwrap();
        writer.write(&record).unwrap();
        drop(writer);

        let (tmp_dst, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">"##,
        ]);
        let src = read_first_record(tmp_src.path());
        let mut dst = create_test_record(&writer, 0, 100, b"A", b"AT");
        copy_info_fields(&src, &mut dst, &["SVTYPE"]).unwrap();
        writer.write(&dst).unwrap();
        drop(writer);

        let result = read_first_record(tmp_dst.path());
        assert_eq!(result.info(b"SVTYPE").string().unwrap().unwrap()[0], b"INS");
    }

    #[test]
    fn test_copy_info_fields_flag() {
        let (tmp_src, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise">"##,
        ]);
        let mut record = create_test_record(&writer, 0, 100, b"A", b"AT");
        record.push_info_flag(b"IMPRECISE").unwrap();
        writer.write(&record).unwrap();
        drop(writer);

        let (tmp_dst, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise">"##,
        ]);
        let src = read_first_record(tmp_src.path());
        let mut dst = create_test_record(&writer, 0, 100, b"A", b"AT");
        copy_info_fields(&src, &mut dst, &["IMPRECISE"]).unwrap();
        writer.write(&dst).unwrap();
        drop(writer);

        let result = read_first_record(tmp_dst.path());
        assert!(result.info(b"IMPRECISE").flag().unwrap());
    }

    #[test]
    fn test_copy_info_fields_absent_silently_skipped() {
        // Source has no SVLEN - copy should succeed, dest field absent
        let (tmp_src, mut writer) = create_info_test_vcf(&[]);
        let record = create_test_record(&writer, 0, 100, b"A", b"AT");
        writer.write(&record).unwrap();
        drop(writer);

        let (tmp_dst, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">"##,
        ]);
        let src = read_first_record(tmp_src.path());
        let mut dst = create_test_record(&writer, 0, 100, b"A", b"AT");

        assert!(copy_info_fields(&src, &mut dst, &["SVLEN"]).is_ok());
        writer.write(&dst).unwrap();
        drop(writer);

        let result = read_first_record(tmp_dst.path());
        assert!(result.info(b"SVLEN").integer().unwrap().is_none());
    }

    #[test]
    fn test_copy_info_fields_number_mismatch_errors() {
        // Header declares Number=1 (exactly one value expected)
        let (tmp_src, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=BADCOUNT,Number=1,Type=Integer,Description="Test">"##,
        ]);
        let mut record = create_test_record(&writer, 0, 100, b"A", b"AT");
        // Write 3 values where the header says exactly 1 is expected
        record.push_info_integer(b"BADCOUNT", &[1, 2, 3]).unwrap();
        writer.write(&record).unwrap();
        drop(writer);

        let (_tmp_dst, writer) = create_info_test_vcf(&[]);
        let src = read_first_record(tmp_src.path());
        let mut dst = create_test_record(&writer, 0, 100, b"A", b"AT");

        let result = copy_info_fields(&src, &mut dst, &["BADCOUNT"]);
        assert!(
            result.is_err(),
            "Number=1 field with 3 stored values should error, got {:?}",
            result
        );
    }

    #[test]
    fn test_copy_info_fields_type_mismatch_errors() {
        // Header declares Type=Integer
        let (tmp_src, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=BADTYPE,Number=1,Type=Integer,Description="Test">"##,
        ]);
        let mut record = create_test_record(&writer, 0, 100, b"A", b"AT");
        // But actually write it as a Float
        record.push_info_float(b"BADTYPE", &[1.5]).unwrap();
        writer.write(&record).unwrap();
        drop(writer);

        let (_tmp_dst, writer) = create_info_test_vcf(&[]);
        let src = read_first_record(tmp_src.path());
        let mut dst = create_test_record(&writer, 0, 100, b"A", b"AT");

        let result = copy_info_fields(&src, &mut dst, &["BADTYPE"]);
        assert!(
            result.is_err(),
            "Integer-declared field actually stored as Float should error, got {:?}",
            result
        );
    }

    #[test]
    fn test_copy_info_fields_multiple() {
        // Integer + String copied in one call
        let (tmp_src, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">"##,
            br##"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">"##,
        ]);
        let mut record = create_test_record(&writer, 0, 100, b"A", b"AT");
        record.push_info_integer(b"SVLEN", &[3]).unwrap();
        record.push_info_string(b"SVTYPE", &[b"INS"]).unwrap();
        writer.write(&record).unwrap();
        drop(writer);

        let (tmp_dst, mut writer) = create_info_test_vcf(&[
            br##"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">"##,
            br##"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">"##,
        ]);
        let src = read_first_record(tmp_src.path());
        let mut dst = create_test_record(&writer, 0, 100, b"A", b"AT");
        copy_info_fields(&src, &mut dst, &["SVLEN", "SVTYPE"]).unwrap();
        writer.write(&dst).unwrap();
        drop(writer);

        let result = read_first_record(tmp_dst.path());
        assert_eq!(result.info(b"SVLEN").integer().unwrap().unwrap()[0], 3);
        assert_eq!(result.info(b"SVTYPE").string().unwrap().unwrap()[0], b"INS");
    }

    /* ======== validate_vcf_file tests ============== */

    #[test]
    fn test_validate_vcf_header_field_correct() {
        let header =
            create_header_view(&[br##"##INFO=<ID=TEST,Number=A,Type=Float,Description="Test">"##]);
        assert!(validate_vcf_header_field(
            header.info_type(b"TEST"),
            "INFO",
            "TEST",
            TagType::Float,
            TagLength::AltAlleles,
        )
        .is_ok());
    }

    #[test]
    fn test_validate_vcf_header_field_wrong_type() {
        let header = create_header_view(&[
            br##"##INFO=<ID=TEST,Number=A,Type=Integer,Description="Test">"##,
        ]);
        let result = validate_vcf_header_field(
            header.info_type(b"TEST"),
            "INFO",
            "TEST",
            TagType::Float,
            TagLength::AltAlleles,
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_vcf_header_field_missing() {
        let header = create_header_view(&[]);
        let result = validate_vcf_header_field(
            header.info_type(b"TEST"),
            "INFO",
            "TEST",
            TagType::Float,
            TagLength::AltAlleles,
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_info_fields_exist_all_present() {
        let header = create_header_view(&[
            br##"##INFO=<ID=CUSTOM1,Number=1,Type=String,Description="Test">"##,
        ]);
        assert!(validate_info_fields_exist(&header, &["CUSTOM1".to_string()]).is_ok());
    }

    #[test]
    fn test_validate_info_fields_exist_missing() {
        let header = create_header_view(&[]);
        let result = validate_info_fields_exist(&header, &["TEST".to_string()]);
        assert!(result.is_err());
        assert!(format!("{}", result.unwrap_err()).contains("TEST"));
    }

    #[test]
    fn test_validate_info_fields_exist_empty_list() {
        let header = create_header_view(&[]);
        assert!(validate_info_fields_exist(&header, &[]).is_ok());
    }

    #[test]
    fn test_validate_info_fields_exist_partial_missing() {
        let header =
            create_header_view(&[br##"##INFO=<ID=OK,Number=1,Type=String,Description="Test">"##]);
        let result =
            validate_info_fields_exist(&header, &["OK".to_string(), "MISSING".to_string()]);
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("MISSING") && !msg.contains("OK"));
    }

    #[test]
    fn test_validate_samples_exist_single_sample() {
        let (tmp_vcf, sample_names) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let result = validate_samples_exist(header, &[sample_names[0].clone()]);

        assert!(result.is_ok());
    }

    #[test]
    fn test_validate_samples_exist_multiple_samples() {
        let (tmp_vcf, sample_names) = create_test_vcf(TestVcfConfig {
            num_samples: 3,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let to_validate = vec![sample_names[0].clone(), sample_names[2].clone()];

        let result = validate_samples_exist(header, &to_validate);

        assert!(result.is_ok());
    }

    #[test]
    fn test_validate_samples_exist_missing_sample() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 2,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let result = validate_samples_exist(header, &["nonexistent_sample".to_string()]);

        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("nonexistent_sample"));
    }

    #[test]
    fn test_validate_samples_exist_multiple_missing() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let result =
            validate_samples_exist(header, &["missing1".to_string(), "missing2".to_string()]);

        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("missing1"));
        assert!(err_msg.contains("missing2"));
    }

    #[test]
    fn test_validate_samples_exist_case_sensitive() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        // Sample is "sample1" (lowercase)
        let result = validate_samples_exist(
            header,
            &["Sample1".to_string()], // Uppercase S
        );

        assert!(result.is_err(), "Sample names are case-sensitive");
    }

    #[test]
    fn test_validate_events_exist_single_event() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let result = validate_events_exist(header, &["somatic".to_string()]);

        assert!(result.is_ok());
    }

    #[test]
    fn test_validate_events_exist_multiple_events() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let result =
            validate_events_exist(header, &["somatic".to_string(), "high_vaf".to_string()]);

        assert!(result.is_ok());
    }

    #[test]
    fn test_validate_events_exist_missing_event() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let result = validate_events_exist(header, &["nonexistent_event".to_string()]);

        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("nonexistent_event"));
        assert!(err_msg.contains("PROB_"));
    }

    #[test]
    fn test_validate_events_exist_case_insensitive() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        // Event field is PROB_SOMATIC (uppercase)
        // User provides lowercase
        let result = validate_events_exist(header, &["somatic".to_string()]);

        assert!(result.is_ok(), "Should handle case conversion");
    }

    #[test]
    fn test_validate_events_exist_multiple_missing() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let result =
            validate_events_exist(header, &["missing1".to_string(), "missing2".to_string()]);

        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("missing1"));
        assert!(err_msg.contains("missing2"));
    }

    #[test]
    fn test_validate_events_exist_partial_match() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        // One valid, one missing
        let result = validate_events_exist(
            header,
            &["somatic".to_string(), "missing_event".to_string()],
        );

        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("missing_event"));
        assert!(!err_msg.contains("somatic")); // Valid one not in error
    }

    #[test]
    fn test_validate_events_exist_wrong_type() {
        let header = create_header_view(&[
            br##"##INFO=<ID=PROB_FOO,Number=A,Type=Integer,Description="Test">"##,
        ]);
        let result = validate_events_exist(&header, &["foo".to_string()]);
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_events_exist_wrong_length() {
        let header = create_header_view(&[
            br##"##INFO=<ID=PROB_FOO,Number=1,Type=Float,Description="Test">"##,
        ]);
        let result = validate_events_exist(&header, &["foo".to_string()]);
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_required_vcf_fields_all_valid() {
        let header =
            create_header_view(&[br##"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Test">"##]);

        assert!(validate_required_vcf_fields_msi(&header).is_ok());
    }

    #[test]
    fn test_validate_required_vcf_fields_missing_field() {
        let header = create_header_view(&[
            br##"##INFO=<ID=PROB_ARTIFACT,Number=A,Type=Float,Description="Test">"##,
        ]);

        let result = validate_required_vcf_fields_msi(&header);
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("AF") && err_msg.contains("missing"));
    }

    #[test]
    fn test_validate_required_vcf_fields_wrong_type() {
        let header = create_header_view(&[
            br##"##FORMAT=<ID=AF,Number=A,Type=Integer,Description="Test">"##,
        ]);

        let result = validate_required_vcf_fields_msi(&header);
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("incorrect type"));
    }

    #[test]
    fn test_validate_vcf_file_with_variant() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"A",
            alt_alleles: vec![b"T"],
            af_values: Some(vec![0.5]),
            num_samples: 1,
            ..Default::default()
        });

        let mut vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();

        let result = validate_vcf_file(&mut vcf);
        assert!(result.is_ok());
    }

    #[test]
    fn test_validate_vcf_file_empty() {
        let tmp = NamedTempFile::new().unwrap();
        let mut header = bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");
        header.push_record(br##"##FORMAT=<ID=AF,Number=A,Type=Float,Description="AF">"##);
        header.push_sample(b"sample1");

        let writer = bcf::Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();
        drop(writer); // No records written

        let mut vcf = bcf::Reader::from_path(tmp.path()).unwrap();

        let result = validate_vcf_file(&mut vcf);
        assert!(result.is_err());

        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("empty") || err_msg.contains("no variant"));
    }

    #[test]
    fn test_validate_vcf_file_corrupted() {
        let tmp = NamedTempFile::new().unwrap();

        // Write valid header but invalid record line
        std::fs::write(
            tmp.path(),
            b"##fileformat=VCFv4.2\n\
            ##contig=<ID=chr1,length=1000>\n\
            ##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n\
            chr1\tINVALID_POS\t.\tA\tT\t.\t.\t.\tAF\t0.5\n",
        )
        .unwrap();

        let mut vcf = bcf::Reader::from_path(tmp.path()).unwrap();

        let result = validate_vcf_file(&mut vcf);
        assert!(result.is_err());

        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("read") || err_msg.contains("parse"));
    }
}
