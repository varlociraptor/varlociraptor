use bio_types::genome::Locus;
use std::path::PathBuf;

use thiserror::Error;

use crate::cli::MIN_THREAD_COUNT;
use crate::estimation::microsatellite_instability as msi;

#[derive(Error, Debug, PartialEq)]
pub(crate) enum Error {
    /* Generic Errors */
    /* ======================== Filing Errors ======================== */
    #[error("invalid BED file path (expected .bed extension): {path}")]
    InvalidBedFile { path: PathBuf },
    #[error(
        "invalid VCF/BCF file path (expected .vcf, .vcf.gz, .bcf, or .bcf.gz extension): {path}"
    )]
    InvalidVcfFile { path: PathBuf },
    /* =============================================================== */
    /* =================== Concurrency Errors ======================== */
    #[error(
        "invalid thread count: must be at least {}, got {count}",
        MIN_THREAD_COUNT
    )]
    InvalidThreadCount { count: usize },
    /* =============================================================== */
    /* ************** */
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
    /* Estimation Errors */
    /* =============== Estimation: MSI Errors ======================== */
    #[error(
        "invalid MSI threshold: must be > {} (default: {}), got {threshold}",
        msi::MIN_MSI_THRESHOLD,
        msi::DEFAULT_MSI_THRESHOLD
    )]
    InvalidMsiThreshold { threshold: f64 },
    #[error("at least one output must be specified: use --plot-pseudotime, --plot-distribution, --data-pseudotime, or --data-distribution")]
    NoMsiOutputSpecified,
    #[error("BED file is empty (no microsatellite loci found)")]
    BedFileEmpty,
    #[error(
        "invalid motif format in BED record: expected 'NxMOTIF' (e.g., '15xCAG'), got '{motif}'"
    )]
    InvalidMsiBedMotif { motif: String },
    #[error("BED record missing required name field (4th column) containing motif information")]
    BedRecordMissingMotifName,
    #[error("VCF/BCF file contains no samples")]
    VcfNoSamples,
    #[error("VCF/BCF file contains invalid sample exclusion(s): sample '{samples}' not found in VCF/BCF")]
    InvalidSampleExclusion { samples: String },
    #[error("VCF/BCF file contains no samples after excluding specified samples")]
    NoSamplesAfterExclusion,
    #[error("VCF/BCF file is empty (no variant records)")]
    VcfEmpty,
    #[error("VCF/BCF record at {chrom}:{pos} is missing chromosome information")]
    VcfRecordMissingChromosome { chrom: String, pos: i64 },
    #[error("VCF/BCF record at position {pos} failed to resolve chromosome name for rid {rid}: {details}")]
    VcfChromResolutionFailed { pos: i64, rid: u32, details: String },
    #[error("failed to read VCF/BCF record: {details}")]
    VcfRecordRead { details: String },
    #[error(
        "Invalid allele frequency {af} for sample '{sample}' at {chrom}:{pos} (must be 0.0-1.0)"
    )]
    InvalidAlleleFrequency {
        sample: String,
        af: f32,
        chrom: String,
        pos: i64,
    },
    #[error("Invalid probability value in field '{field}' at {chrom}:{pos} (value={value})")]
    InvalidProbabilityValue {
        field: String,
        value: f32,
        chrom: String,
        pos: i64,
    },
    #[error("No chromosome match between BED and VCF files. Ensure chromosome naming is consistent (both with or without 'chr' prefix).")]
    NoChromosomeMatch,
    #[error("Failed to read BED record at line {line}: {details}")]
    BedRecordRead { line: usize, details: String },
}

pub(crate) fn invalid_bcf_record(chrom: &str, pos: i64, msg: &str) -> Error {
    Error::InvalidBCFRecord {
        chrom: chrom.to_owned(),
        pos,
        msg: msg.to_owned(),
    }
}
