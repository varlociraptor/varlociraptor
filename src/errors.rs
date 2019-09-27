pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Snafu, Debug, PartialEq)]
#[snafu(visibility = "pub")]
pub enum Error {
    #[snafu(display("formula refers to unknown sample {}", name))]
    InvalidSampleName { name: String },
    #[snafu(display("event {} does not define VAF range for all samples", event_name))]
    MissingSampleEvent { event_name: String },
    #[snafu(display("no BAM file given for sample {}", name))]
    InvalidBAMSampleName { name: String },
    #[snafu(display(
        "contamination refers to unknown sample {}; it is not defined in the scenario",
        name
    ))]
    InvalidContaminationSampleName { name: String },
    #[snafu(display("alignment property files must be provided as name=path"))]
    InvalidAlignmentPropertiesSpec,
    #[snafu(display("BAM files must be provided as name=path"))]
    InvalidBAMSpec,
    #[snafu(display(
        "invalid variant index given, must be not higher than the number of variants at the locus"
    ))]
    InvalidIndex,
    #[snafu(display("invalid locus for --testcase-locus. Use CHROM:POS syntax"))]
    InvalidLocus,
    #[snafu(display("no candidate variant at the given locus"))]
    NoCandidateFound,
    #[snafu(display("testcase prefix must be given with --testcase-prefix"))]
    MissingPrefix,
    #[snafu(display("candidate variants must be provided via --candidates"))]
    MissingCandidates,
    #[snafu(display("--min-bayes-factor must be between 0.0 and 1.0"))]
    InvalidMinBayesFactor,
    #[snafu(display("expected tag {} missing from BCF record", name))]
    MissingBCFTag { name: String },
    #[snafu(display("invalid BCF record: {}", msg))]
    InvalidBCFRecord { msg: String },
    #[snafu(display(
        "unable to estimate TMB because no valid records were found in the given BCF/VCF"
    ))]
    NoRecordsFound,
}
