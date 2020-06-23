use std::path::PathBuf;

#[derive(Snafu, Debug, PartialEq)]
#[snafu(visibility = "pub(crate)")]
pub(crate) enum Error {
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
    #[snafu(display("observation files must be provided as samplename=path"))]
    InvalidObservationsSpec,
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
    #[snafu(display(
        "contig {} not found in universe definition and no 'all' defined",
        contig
    ))]
    UniverseContigNotFound { contig: String },
    #[snafu(display("contig {} not found in the reference", contig))]
    ReferenceContigNotFound { contig: String },
    #[snafu(display("record {} in candidate BCF/VCF does not define a chromosome", i))]
    RecordMissingChrom { i: usize },
    #[snafu(display("inconsistent observations: input observation BCF files do not contain exactly the same records"))]
    InconsistentObservations,
    #[snafu(display("No observations given for sample {}.", name))]
    InvalidObservationSampleName { name: String },
    #[snafu(display("invalid observations: varlociraptor cannot be parsed from given observations ({}); either the file has not been preprocessed with varlociraptor or with a too old version", path.display()))]
    InvalidObservations { path: PathBuf },
    #[snafu(display("invalid observations: varlociraptor cannot read given observations; either the file has not been preprocessed with varlociraptor or with a too old version"))]
    InvalidObservationFormat,
}
