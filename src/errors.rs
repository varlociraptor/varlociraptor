//! src/errors.rs
//!
//!  Error types for varlociraptor
//!
//! # Naming Convention
//!
//! All error variants should follow a consistent naming pattern:
//! **`<Component><Action><ErrorType>`**
//!
//! Note: Please register new components, actions, and error types when adding errors,
//! if they do not already exist, so that the naming convention remains consistent and
//! unrepetitive.
//!
//! ## Components:
//! - `Bed` - BED file operations (generic)
//! - `Vcf` - VCF/BCF file operations (generic)
//! - `Thread` - Threading and concurrency
//! - `MsiBed` - MSI-specific BED operations
//! - `MsiConfig` - MSI configuration
//! - `MsiVcf` - MSI-specific VCF operations

//!
//! ## Actions(when applicable):
//! - `File` - File-level operations
//! - `Record` - Record-level operations
//! - `Chrom` - Chromosome operations
//! - `Samples` - Sample operations
//! - `Sample` - Single sample operations
//! - `Allele` - Allele operations
//! - `Frequency` - Frequency calculations
//! - `Probability` - Probability calculations
//! - `Exclusion` - Sample exclusion operations
//! - `Count` - Counting operations
//! - `Threshold` - Threshold validation
//! - `Motif` - Motif pattern operations (MSI)
//! - `Output` - Output configuration (MSI)
//!
//! ## Error Types:
//! - `Invalid` - Validation failure
//! - `Missing` - Required data absent
//! - `Empty` - No data found
//! - `Failed` - Operation failure
//! - `Mismatch` - Data inconsistency
//!
//! ## Examples:
//! - `BedFileInvalid` - Invalid BED file path
//! - `VcfRecordReadFailed` - Failed to read VCF record
//! - `MsiConfigThresholdInvalid` - Invalid MSI threshold value

use bio_types::genome::Locus;
use std::path::PathBuf;

use thiserror::Error;

use crate::cli::MIN_THREAD_COUNT;
use crate::estimation::microsatellite_instability as msi;

#[derive(Error, Debug, PartialEq)]
pub(crate) enum Error {
    /* ======================= Generic Errors ======================== */
    /* -------------------- File Validation -------------------------- */
    /* 1. Bed File Errors */
    #[error("invalid BED file path (expected .bed extension): {path}")]
    BedFileInvalid { path: PathBuf },
    #[error("BED file is empty (no microsatellite loci found)")]
    BedFileEmpty,
    #[error("invalid BED record at {chrom}:{pos}: {msg}")]
    BedRecordInvalid {
        chrom: String,
        pos: i64,
        msg: String,
    },
    #[error("failed to read BED record at line {line}: {details}")]
    BedRecordReadFailed { line: usize, details: String },
    /* 2. VCF/BCF File Errors */
    #[error(
        "invalid VCF/BCF file path (expected .vcf, .vcf.gz, .bcf, or .bcf.gz extension): {path}"
    )]
    VcfFileInvalid { path: PathBuf },
    #[error("VCF/BCF file contains no samples")]
    VcfSamplesMissing,
    #[error("VCF/BCF file is empty (no variant records)")]
    VcfFileEmpty,
    #[error("VCF/BCF record at position {pos} is missing chromosome reference (RID)")]
    VcfRecordChromMissing { pos: i64 },
    #[error("VCF/BCF record at position {pos} failed to resolve chromosome name for rid {rid}: {details}")]
    VcfRecordChromResolveFailed { pos: i64, rid: u32, details: String },
    #[error("failed to read VCF/BCF record: {details}")]
    VcfRecordReadFailed { details: String },
    #[error(
        "Invalid allele frequency {af} for sample '{sample}' at {chrom}:{pos} (must be 0.0-1.0)"
    )]
    VcfAlleleFrequencyInvalid {
        sample: String,
        af: f32,
        chrom: String,
        pos: i64,
    },
    #[error("Invalid probability value in field '{field}' at {chrom}:{pos} (value={value})")]
    VcfProbabilityValueInvalid {
        field: String,
        value: f32,
        chrom: String,
        pos: i64,
    },
    #[error("VCF/BCF file does not contain the following excluded sample(s): {samples}")]
    VcfSampleExclusionInvalid { samples: String },
    #[error("VCF/BCF file contains no samples after excluding specified samples")]
    VcfSamplesEmptyAfterExclusion,
    #[error("VCF/BCF header missing or malformed required {location} field: {field}")]
    VcfHeaderFieldMissing {
        field: String,
        location: String, // e.g "INFO" or "FORMAT"
    },
    #[error("VCF/BCF header field {location}:{field} has incorrect type (expected {expected}, found {found})")]
    VcfHeaderFieldTypeInvalid {
        location: String, // e.g "INFO" or "FORMAT"
        field: String,
        expected: String,
        found: String,
    },
    /* -------------------- Concurrency ------------------------------ */
    #[error(
        "invalid thread count: must be at least {}, got {count}",
        MIN_THREAD_COUNT
    )]
    ThreadCountInvalid { count: usize },
    /* --------------------------------------------------------------- */
    #[error("formula refers to unknown sample {name}")]
    InvalidSampleName { name: String },
    #[error("contamination refers to unknown sample {name}; it is not defined in the scenario")]
    InvalidContaminationSampleName { name: String },
    #[error(
        "inheritance definition refers to unknown sample {name}; it is not defined in the scenario"
    )]
    InvalidInheritanceSampleName { name: String },
    #[error("observation files must be provided as samplename=path")]
    InvalidObservationsSpec,
    #[error("invalid locus for --testcase-locus. Use CHROM:POS syntax")]
    InvalidLocus,
    #[error("no candidate variant at the given locus")]
    NoCandidateFound,
    #[error("testcase prefix must be given with --testcase-prefix")]
    MissingPrefix,
    #[error("expected tag {name} missing from BCF record")]
    MissingBCFTag { name: String },
    #[error("invalid BCF record at {chrom}:{pos}: {msg}")]
    InvalidBCFRecord {
        chrom: String,
        pos: i64,
        msg: String,
    },
    #[error("unable to estimate TMB because no valid records were found in the given BCF/VCF")]
    NoRecordsFound,
    #[error("contig {contig} not found in universe definition and no 'all' defined")]
    UniverseContigNotFound { contig: String },
    #[error("contig {contig} not found in ploidy definition and no 'all' defined")]
    PloidyContigNotFound { contig: String },
    #[error("record {i} in candidate BCF/VCF does not define a chromosome")]
    RecordMissingChrom { i: usize },
    #[error("inconsistent observations: input observation BCF files do not contain exactly the same records")]
    InconsistentObservations,
    #[error("sample {name} (given by --obs) cannot be found in the scenario")]
    InvalidObservationSampleName { name: String },
    #[error("invalid observations: varlociraptor cannot be parsed from given observations ({path}); either the file has not been preprocessed with varlociraptor or with a too old version")]
    InvalidObservations { path: PathBuf },
    #[error("invalid observations: varlociraptor cannot read given observations; either the file has not been preprocessed with varlociraptor or with a too old version")]
    InvalidObservationFormat,
    #[error("invalid BND record: ALT {spec} does not follow BND spec")]
    InvalidBNDRecordAlt { spec: String },
    #[error("at least one BCF with observations must be provided")]
    EmptyObservations,
    #[error(
        "undefined expression {identifier}; please define under 'expressions:' in your scenario"
    )]
    UndefinedExpression { identifier: String },
    #[error("invalid prior configuration: {msg}")]
    InvalidPriorConfiguration { msg: String },
    #[error("read position determined from cigar string exceeds record length")]
    ReadPosOutOfBounds,
    #[error("invalid strand information '{value}', must be '+', '-', '*' or '.'")]
    InvalidStrandInfo { value: char },
    #[error("invalid read orientation information '{value}', must be 'F1R2', 'F2R1', etc.")]
    InvalidReadOrientationInfo { value: String },
    #[error("the following events are not disjunct: {expressions}")]
    OverlappingEvents { expressions: String },
    #[error("the input VCF/BCF is not sorted: {previous_locus:?} > {current_locus:?}")]
    UnsortedVariantFile {
        previous_locus: Locus,
        current_locus: Locus,
    },
    // #[error("invalid phase set, PS tag only supported for single sample VCF/BCF, may only contain a single value, and records may only contain a single ALT allele")]
    // InvalidPhaseSet,
    #[error("haplotype block consisting of normal variants in combination with breakends: this is currently unsupported")]
    HaplotypeBlockWithBreakend,
    #[error("invalid prior contamination estimate. Both --prior-estiate and --prior-considered-cells have to be specified. The latter has to be >0.")]
    InvalidPriorContaminationEstimate,
    #[error("breakend with MATEID found that does not have its own ID set: this is currently unsupported, as there is no way to uniquely identify the pair")]
    BreakendMateidWithoutRecid,
    #[error("invalid FDR control events, no events provided or none of the given events matches the events found in the callset")]
    InvalidFDRControlEvents,
    #[error("unrealistic insert size distribution: the standard deviation is 0.0, consider sampling more reads for estimating alignment properties")]
    UnrealisticIsizeSd,
    #[error("given field for variant heterozygosity or variant somatic effective mutation rate has to have as many entries as ALT alleles in the record")]
    InvalidVariantPrior,

    /* ======================= MSI: Estimation Errors ================ */
    /* -------------------- Configuration ---------------------------- */
    #[error(
        "invalid MSI threshold: must be > {} (default: {}), got {threshold}",
        msi::MIN_MSI_THRESHOLD,
        msi::DEFAULT_MSI_THRESHOLD
    )]
    MsiConfigThresholdInvalid { threshold: f64 },
    #[error("at least one output must be specified: use --plot-pseudotime, --plot-distribution, --data-pseudotime, or --data-distribution")]
    MsiConfigOutputMissing,
    /* -------------------- BED File Errors -------------------------- */
    #[error(
        "invalid motif format in BED record: expected 'NxMOTIF' (e.g., '15xCAG'), got '{motif}'"
    )]
    MsiBedMotifInvalid { motif: String },
    #[error("BED record missing required name field (4th column) containing motif information")]
    MsiBedMotifNameMissing,
    /* -------------------- Processing Errors ------------------------ */
    #[error("No chromosome match between BED and VCF files. Verify chromosome naming is consistent (e.g., 'chr1' vs '1')")]
    MsiVcfChromMismatch,
}

pub(crate) fn invalid_bcf_record(chrom: &str, pos: i64, msg: &str) -> Error {
    Error::InvalidBCFRecord {
        chrom: chrom.to_owned(),
        pos,
        msg: msg.to_owned(),
    }
}
