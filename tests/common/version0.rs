use std::path::PathBuf;
use std::str;

use serde_json;
use yaml_rust::Yaml;

use crate::common::Testcase;
use varlociraptor::cli::{PreprocessKind, Varlociraptor};
use varlociraptor::testcase::Mode;

#[derive(Debug)]
pub(crate) struct TestcaseVersion0 {
    pub(crate) inner: Vec<Yaml>,
    pub(crate) path: PathBuf,
}

impl TestcaseVersion0 {
    fn options(&self) -> cli::Varlociraptor {
        serde_json::from_str(self.yaml()["options"].as_str().unwrap()).unwrap()
    }
}

impl Testcase for TestcaseVersion0 {
    fn inner(&self) -> &[Yaml] {
        &self.inner
    }

    fn path(&self) -> &PathBuf {
        &self.path
    }

    fn mode(&self) -> Mode {
        let options = self.options();
        match options {
            cli::Varlociraptor::Call {
                kind: cli::CallKind::Variants { mode, .. },
            } => match mode {
                cli::VariantCallMode::Generic { .. } => Mode::Generic,
                cli::VariantCallMode::TumorNormal { .. } => Mode::TumorNormal,
            },
            _ => panic!("unexpected options"),
        }
    }

    fn sample_alignment_properties(&self, sample_name: &str) -> String {
        let mut props: serde_json::Value =
            serde_json::from_str(self.sample(sample_name)["properties"].as_str().unwrap()).unwrap();
        props.as_object_mut().unwrap().insert(
            "max_read_len".to_owned(),
            serde_json::Value::Number(serde_json::Number::from(100)),
        );

        props.to_string()
    }

    fn purity(&self) -> Option<f64> {
        match self.options() {
            cli::Varlociraptor::Call {
                kind:
                    cli::CallKind::Variants {
                        mode: cli::VariantCallMode::TumorNormal { purity, .. },
                        ..
                    },
            } => Some(purity),
            _ => panic!("Invalid options"),
        }
    }

    fn preprocess_options(&self, _: &str) -> String {
        let options = self.options();

        match options {
            cli::Varlociraptor::Call {
                kind:
                    cli::CallKind::Variants {
                        reference,
                        spurious_ins_rate,
                        spurious_del_rate,
                        spurious_insext_rate,
                        spurious_delext_rate,
                        protocol_strandedness,
                        indel_window,
                        max_depth,
                        ..
                    },
            } => {
                let options = Varlociraptor::Preprocess {
                    kind: PreprocessKind::Variants {
                        reference,
                        spurious_ins_rate,
                        spurious_del_rate,
                        spurious_insext_rate,
                        spurious_delext_rate,
                        protocol_strandedness,
                        realignment_window: indel_window as u64,
                        max_depth,
                        // The rest will be overwritten.
                        alignment_properties: None,
                        bam: PathBuf::from("dummy"),
                        candidates: self.candidates(),
                        output: None,
                        omit_insert_size: false,
                        threads: 1,
                        buffer_capacity: 1,
                    },
                };

                serde_json::to_string(&options).unwrap()
            }
            _ => panic!("invalid options for testcase"),
        }
    }
}

// old cli
pub(crate) mod cli {
    use std::path::PathBuf;

    use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
    use itertools::Itertools;
    use serde::{Deserialize, Serialize};
    use structopt::StructOpt;

    use varlociraptor::variants::model::VariantType;
    use varlociraptor::variants::sample::ProtocolStrandedness;

    #[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
    #[structopt(
        name = "varlociraptor",
        about = "A caller for SNVs and indels in tumor-normal pairs.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    pub(crate) enum Varlociraptor {
        #[structopt(
            name = "call",
            about = "Call variants.",
            setting = structopt::clap::AppSettings::ColoredHelp,
        )]
        Call {
            #[structopt(subcommand)]
            kind: CallKind,
        },
        #[structopt(
            name = "filter-calls",
            about = "Filter calls by either controlling the false discovery rate (FDR) at given level, or by posterior odds against the given events.",
            setting = structopt::clap::AppSettings::ColoredHelp,
        )]
        FilterCalls {
            #[structopt(subcommand)]
            method: FilterMethod,
        },
        #[structopt(
            name = "decode-phred",
            about = "Decode PHRED-scaled values to human readable probabilities.",
            setting = structopt::clap::AppSettings::ColoredHelp,
        )]
        DecodePHRED,
        #[structopt(
            name = "estimate",
            about = "Perform estimations.",
            setting = structopt::clap::AppSettings::ColoredHelp,
        )]
        #[structopt(
            name = "estimate",
            about = "Perform estimations.",
            setting = structopt::clap::AppSettings::ColoredHelp,
        )]
        Estimate {
            #[structopt(subcommand)]
            kind: EstimateKind,
        },
    }

    #[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
    pub(crate) enum EstimateKind {
        #[structopt(
            name = "tmb",
            about = "Estimate tumor mutational burden. Takes Varlociraptor calls (must be annotated \
                    with e.g. snpEFF) from STDIN, prints TMB estimate in Vega-lite JSON format to STDOUT. \
                    It can be converted to an image via vega-lite-cli (see conda package).",
            usage = "varlociraptor estimate tmb --coding-genome-size 3e7 --somatic-tumor-events SOMATIC_TUMOR \
                    --tumor-sample tumor < calls.bcf | vg2svg > tmb.svg",
            setting = structopt::clap::AppSettings::ColoredHelp,
        )]
        TMB {
            #[structopt(
                long = "somatic-tumor-events",
                default_value = "SOMATIC_TUMOR",
                help = "Events to consider (e.g. SOMATIC_TUMOR)."
            )]
            somatic_tumor_events: Vec<String>,
            #[structopt(
                long = "tumor-sample",
                default_value = "tumor",
                help = "Name of the tumor sample in the given VCF/BCF."
            )]
            tumor_sample: String,
            #[structopt(
                long = "coding-genome-size",
                default_value = "3e7",
                help = "Size (in bases) of the covered coding genome."
            )]
            coding_genome_size: f64,
        },
    }

    #[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
    pub(crate) enum CallKind {
        #[structopt(
            name = "variants",
            about = "Call variants.",
            setting = structopt::clap::AppSettings::ColoredHelp,
        )]
        Variants {
            #[structopt(subcommand)]
            mode: VariantCallMode,
            #[structopt(
                parse(from_os_str),
                help = "FASTA file with reference genome. Has to be indexed with samtools faidx."
            )]
            reference: PathBuf,
            #[structopt(
                parse(from_os_str),
                long,
                help = "VCF/BCF file to process (if omitted, read from STDIN)."
            )]
            candidates: Option<PathBuf>,
            #[structopt(
                parse(from_os_str),
                long,
                help = "BCF file that shall contain the results (if omitted, write to STDOUT)."
            )]
            output: Option<PathBuf>,
            #[structopt(
                long = "spurious-ins-rate",
                default_value = "2.8e-6",
                help = "Rate of spuriously inserted bases by the sequencer (Illumina: 2.8e-6, see Schirmer et al. BMC Bioinformatics 2016)."
            )]
            spurious_ins_rate: f64,
            #[structopt(
                long = "spurious-del-rate",
                default_value = "5.1e-6",
                help = "Rate of spuriosly deleted bases by the sequencer (Illumina: 5.1e-6, see Schirmer et al. BMC Bioinformatics 2016)."
            )]
            spurious_del_rate: f64,
            #[structopt(
                long = "spurious-insext-rate",
                default_value = "0.0",
                help = "Extension rate of spurious insertions by the sequencer (Illumina: 0.0, see Schirmer et al. BMC Bioinformatics 2016)"
            )]
            spurious_insext_rate: f64,
            #[structopt(
                long = "spurious-delext-rate",
                default_value = "0.0",
                help = "Extension rate of spurious deletions by the sequencer (Illumina: 0.0, see Schirmer et al. BMC Bioinformatics 2016)"
            )]
            spurious_delext_rate: f64,
            #[structopt(
                long = "strandedness",
                default_value = "opposite",
                possible_values = { use strum::IntoEnumIterator; &ProtocolStrandedness::iter().map(|v| v.into()).collect_vec() },
                help = "Strandedness of sequencing protocol in case of paired-end (opposite strand as usual or same strand as with mate-pair sequencing.)"
            )]
            protocol_strandedness: ProtocolStrandedness,
            #[structopt(long = "omit-snvs", help = "Don't call SNVs.")]
            omit_snvs: bool,
            #[structopt(long = "omit-indels", help = "Don't call Indels.")]
            omit_indels: bool,
            #[structopt(
                long = "max-indel-len",
                default_value = "1000",
                help = "Omit longer indels when calling."
            )]
            max_indel_len: u32,
            #[structopt(
                long = "indel-window",
                default_value = "64",
                help = "Number of bases to consider left and right of indel breakpoint when \
                        calculating read support. This number should not be too large in order to \
                        avoid biases caused by other close variants. Currently implemented maximum \
                        value is 64."
            )]
            indel_window: u32,
            #[structopt(
                long = "max-depth",
                default_value = "200",
                help = "Maximum number of observations to use for calling. If locus is exceeding this \
                        number, downsampling is performed."
            )]
            max_depth: usize,
            #[structopt(
                long = "testcase-locus",
                help = "Create a test case for the given locus. Locus must be given in the form \
                        CHROM:POS[:IDX]. IDX is thereby an optional value to select a particular \
                        variant at the locus, counting from 1. If IDX is not specified, the first \
                        variant will be chosen. Alternatively, for single variant VCFs, you can \
                        specify 'all'."
            )]
            testcase_locus: Option<String>,
            #[structopt(
                long = "testcase-prefix",
                help = "Create test case files in the given directory."
            )]
            testcase_prefix: Option<String>,
        },
        #[structopt(
            name = "cnvs",
            about = "Call CNVs in tumor-normal sample pairs. This is experimental (do not use it yet).",
            setting = structopt::clap::AppSettings::ColoredHelp,
        )]
        CNVs {
            #[structopt(
                parse(from_os_str),
                long,
                help = "VCF/BCF file (generated by varlociraptor call-tumor-normal) to process \
                        (if omitted, read from STDIN)."
            )]
            calls: Option<PathBuf>,
            #[structopt(
                parse(from_os_str),
                long,
                help = "BCF file that shall contain the results (if omitted, write to STDOUT)."
            )]
            output: Option<PathBuf>,
            #[structopt(long, short = "p", help = "Tumor purity.")]
            purity: f64,
            #[structopt(
                long = "min-bayes-factor",
                default_value = "1.01",
                help = "Minimum bayes factor (> 1.0) between likelihoods of CNV and no CNV to consider. \
                        The higher this value, the fewer candidate CNVs will be investigated. \
                        Note that this can be usually left unchanged, because every CNV is provided \
                        with a posterior probability that can be used for filtering, e.g., via \
                        'varlociraptor filter-calls control-fdr'."
            )]
            min_bayes_factor: f64,
            #[structopt(
                long,
                default_value = "1000",
                help = "Maximum distance between supporting loci in a CNV."
            )]
            max_dist: u32,
            #[structopt(long, short = "t", help = "Number of threads to use.")]
            threads: usize,
        },
    }

    #[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
    pub(crate) enum VariantCallMode {
        #[structopt(
            name = "tumor-normal",
            about = "Call somatic and germline variants from a tumor-normal sample pair and a VCF/BCF with candidate variants.",
            usage = "varlociraptor call variants reference.fa tumor-normal --purity 0.75 tumor.bam normal.bam < candidates.bcf > calls.bcf",
            setting = structopt::clap::AppSettings::ColoredHelp,
        )]
        TumorNormal {
            #[structopt(parse(from_os_str), help = "BAM file with reads from tumor sample.")]
            tumor: PathBuf,
            #[structopt(parse(from_os_str), help = "BAM file with reads from normal sample.")]
            normal: PathBuf,
            #[structopt(short, long, default_value = "1.0", help = "Purity of tumor sample.")]
            purity: f64,
            #[structopt(
                parse(from_os_str),
                long = "tumor-alignment-properties",
                help = "Alignment properties JSON file for tumor sample. If not provided, properties \
                        will be estimated from the given BAM file."
            )]
            tumor_alignment_properties: Option<PathBuf>,
            #[structopt(
                parse(from_os_str),
                long = "normal-alignment-properties",
                help = "Alignment properties JSON file for normal sample. If not provided, properties \
                        will be estimated from the given BAM file."
            )]
            normal_alignment_properties: Option<PathBuf>,
        },
        #[structopt(
            name = "generic",
            about = "Call variants for a given scenario specified with the varlociraptor calling \
                    grammar and a VCF/BCF with candidate variants.",
            usage = "varlociraptor call variants reference.fa generic --bams relapse=relapse.bam \
                    tumor=tumor.bam normal=normal.bam < candidates.bcf > calls.bcf",
            setting = structopt::clap::AppSettings::ColoredHelp,
        )]
        Generic {
            #[structopt(
                parse(from_os_str),
                long,
                help = "Scenario defined in the varlociraptor calling grammar."
            )]
            scenario: PathBuf,
            #[structopt(long, help = "BAM files with aligned reads for each sample.")]
            bams: Vec<String>,
            #[structopt(
                long = "alignment-properties",
                help = "Alignment properties JSON file for normal sample. If not provided, properties \
                        will be estimated from the given BAM file."
            )]
            alignment_properties: Vec<String>,
        },
    }

    #[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
    pub(crate) enum FilterMethod {
        #[structopt(
            name = "control-fdr",
            about = "Filter variant calls by controlling FDR. Filtered calls are printed to STDOUT.",
            usage = "varlociraptor filter-calls control-fdr calls.bcf --events SOMATIC_TUMOR --fdr 0.05 \
                    --var SNV > calls.filtered.bcf",
            setting = structopt::clap::AppSettings::ColoredHelp,
        )]
        ControlFDR {
            #[structopt(parse(from_os_str), help = "BCF file with varlociraptor calls.")]
            calls: PathBuf,
            #[structopt(
                long = "var",
                possible_values = { use strum::IntoEnumIterator; &VariantType::iter().map(|v| v.into()).collect_vec() },
                help = "Variant type to consider."
            )]
            vartype: VariantType,
            #[structopt(long, help = "FDR to control for.")]
            fdr: f64,
            #[structopt(long, help = "Events to consider.")]
            events: Vec<String>,
            #[structopt(long, help = "Minimum indel length to consider.")]
            minlen: Option<u32>,
            #[structopt(long, help = "Maximum indel length to consider (exclusive).")]
            maxlen: Option<u32>,
        },
        #[structopt(
            name = "posterior-odds",
            about = "Filter variant calls by posterior odds of given events against the rest of events. \
                    Calls are taken from STDIN, filtered calls are printed to STDOUT.",
            usage = "varlociraptor filter-calls posterior-odds --events SOMATIC_TUMOR --odds strong < calls.bcf",
            setting = structopt::clap::AppSettings::ColoredHelp,
        )]
        PosteriorOdds {
            #[structopt(
                possible_values = { use strum::IntoEnumIterator; &KassRaftery::iter().map(|v| v.into()).collect_vec() },
                help = "Kass-Raftery score to filter against."
            )]
            odds: KassRaftery,
            #[structopt(long, help = "Events to consider.")]
            events: Vec<String>,
        },
    }
}
