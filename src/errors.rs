use std::path::PathBuf;

use thiserror::Error;

#[derive(Error, Debug, PartialEq)]
pub(crate) enum Error {
    #[error("formula refers to unknown sample {name}")]
    InvalidSampleName { name: String },
    #[error("contamination refers to unknown sample {name}; it is not defined in the scenario")]
    InvalidContaminationSampleName { name: String },
    #[error("observation files must be provided as samplename=path")]
    InvalidObservationsSpec,
    #[error(
        "invalid variant index given, must be not higher than the number of variants at the locus"
    )]
    InvalidIndex,
    #[error("invalid locus for --testcase-locus. Use CHROM:POS syntax")]
    InvalidLocus,
    #[error("no candidate variant at the given locus")]
    NoCandidateFound,
    #[error("testcase prefix must be given with --testcase-prefix")]
    MissingPrefix,
    #[error("expected tag {name} missing from BCF record")]
    MissingBCFTag { name: String },
    #[error("invalid BCF record: {msg}")]
    InvalidBCFRecord { msg: String },
    #[error("unable to estimate TMB because no valid records were found in the given BCF/VCF")]
    NoRecordsFound,
    #[error("contig {contig} not found in universe definition and no 'all' defined")]
    UniverseContigNotFound { contig: String },
    #[error("record {i} in candidate BCF/VCF does not define a chromosome")]
    RecordMissingChrom { i: usize },
    #[error("inconsistent observations: input observation BCF files do not contain exactly the same records")]
    InconsistentObservations,
    #[error("no observations given for sample {name}")]
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
        "undefined expression ${identifier}; please define under 'expressions:' in your scenario"
    )]
    UndefinedExpression { identifier: String },
    #[error("read position determined from cigar string exceeds record length")]
    ReadPosOutOfBounds,
    #[error("invalid strand information '{value}', must be '+', '-', or '*'")]
    InvalidStrandInfo { value: char },
    #[error("invalid read orientation information '{value}', must be 'F1R2', 'F2R1', etc.")]
    InvalidReadOrientationInfo { value: String }
}
