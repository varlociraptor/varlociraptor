// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cell::RefCell;
use std::collections::HashMap;
use std::convert::{From, TryFrom};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{bail, Context, Result};
use bio::io::fasta;
use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::{LogProb, Prob};
use itertools::Itertools;
use structopt::StructOpt;
use strum::IntoEnumIterator;

use crate::calling;
use crate::calling::variants::calling::{
    call_generic, CallWriter, DefaultCandidateFilter, SampleInfos,
};
use crate::calling::variants::preprocessing::haplotype_feature_index::HaplotypeFeatureIndex;
use crate::conversion;
use crate::errors;
use crate::estimation;
use crate::estimation::alignment_properties::AlignmentProperties;
//use crate::estimation::sample_variants;
//use crate::estimation::tumor_mutational_burden;
use crate::filtration;
use crate::grammar;
use crate::reference;
use crate::testcase;
use crate::utils::PathMap;
use crate::variants::evidence::realignment;

use crate::variants::model::prior::CheckablePrior;
use crate::variants::model::prior::{Inheritance, Prior};
use crate::variants::model::{AlleleFreq, Contamination, VariantType};
use crate::variants::sample::{estimate_alignment_properties, ProtocolStrandedness};
use crate::SimpleEvent;

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
#[structopt(
    name = "varlociraptor",
    about = "Flexible, arbitrary-scenario, uncertainty-aware variant calling \
    with parameter free filtration via FDR control for small and structural variants.",
    setting = structopt::clap::AppSettings::ColoredHelp,
)]
pub enum Varlociraptor {
    #[structopt(
        name = "preprocess",
        about = "Preprocess variants",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Preprocess {
        #[structopt(subcommand)]
        kind: PreprocessKind,
    },
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
        usage = "varlociraptor decode-phred < in.bcf > out.bcf",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    DecodePHRED,
    #[structopt(
        name = "estimate",
        about = "Perform estimations.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Estimate {
        #[structopt(subcommand)]
        kind: EstimateKind,
    },
    #[structopt(
        name = "plot",
        about = "Create plots",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Plot {
        #[structopt(subcommand)]
        kind: PlotKind,
    },
    #[structopt(
        name = "genotype",
        about = "Infer classical genotypes from Varlociraptor's AF field (1.0: 1/1, 0.5: 0/1, 0.0: 0/0, otherwise: 0/1). This assumes diploid samples.",
        usage = "varlociraptor genotype < in.bcf > out.bcf",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Genotype,
}

pub struct PreprocessInput {
    reference: PathBuf,
    bam: PathBuf,
}

impl Varlociraptor {
    fn preprocess_input(&self) -> PreprocessInput {
        if let Varlociraptor::Preprocess {
            kind:
                PreprocessKind::Variants {
                    ref reference,
                    ref bam,
                    ..
                },
        } = &self
        {
            PreprocessInput {
                reference: reference.to_owned(),
                bam: bam.to_owned(),
            }
        } else {
            panic!("bug: these are not preprocess options.");
        }
    }
}

fn default_reference_buffer_size() -> usize {
    10
}

fn default_pairhmm_mode() -> String {
    "exact".to_owned()
}

fn default_log_mode() -> String {
    "default".to_owned()
}

fn default_min_bam_refetch_distance() -> u64 {
    1
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum PreprocessKind {
    #[structopt(
        name = "variants",
        about = "Preprocess given variants by obtaining internal observations (allele likelihoods, strand information, ...)\
                 for each fragment. \
                 The obtained observations are printed to STDOUT in BCF format. Note that the resulting BCFs \
                 will be very large and are only intended for internal use (e.g. for piping into 'varlociraptor \
                 call variants generic').",
        usage = "varlociraptor preprocess variants reference.fasta --alignment-properties alignment-properties.json \
                 --candidates candidates.bcf --bam sample.bam > sample.observations.bcf",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Variants {
        #[structopt(
            parse(from_os_str),
            help = "FASTA file with reference genome. Has to be indexed with samtools faidx."
        )]
        reference: PathBuf,
        #[structopt(
            parse(from_os_str),
            long,
            required = true,
            help = "sorted VCF/BCF file to process (if omitted, read from STDIN)."
        )]
        candidates: PathBuf,
        #[structopt(
            long,
            required = true,
            help = "BAM file with aligned reads from a single sample."
        )]
        bam: PathBuf,
        #[structopt(
            long,
            help = "Report fragment IDs in output BCF. This information can be used for phasing."
        )]
        #[serde(default)]
        report_fragment_ids: bool,
        #[structopt(
            long = "reference-buffer-size",
            short = "b",
            default_value = "10",
            help = "Number of reference sequences to keep in buffer. Use a smaller value \
                    to save memory at the expense of sometimes reduced parallelization."
        )]
        #[serde(default = "default_reference_buffer_size")]
        reference_buffer_size: usize,
        #[structopt(
            long = "min-bam-refetch-distance",
            default_value = "1",
            help = "Base pair distance to last fetched BAM interval such that a refetching is performed \
                  instead of reading through until the next interval is reached. Making this too small \
                  can cause unnecessary random access. Making this too large can lead to unneccessary \
                  iteration over irrelevant records. Benchmarking has shown that at least for short reads, \
                  a value of 1 (e.g. always refetch) does not incur additional costs and is a reasonable \
                  default."
        )]
        #[serde(default = "default_min_bam_refetch_distance")]
        min_bam_refetch_distance: u64,
        #[structopt(
            long = "alignment-properties",
            help = "Alignment properties JSON file for sample. If not provided, properties \
                    will be estimated from the given BAM file. It is recommended to estimate alignment \
                    properties separately, see 'varlociraptor estimate alignment-properties --help'."
        )]
        alignment_properties: Option<PathBuf>,
        #[structopt(
            parse(from_os_str),
            long,
            help = "BCF file that shall contain the results (if omitted, write to STDOUT)."
        )]
        output: Option<PathBuf>,
        #[structopt(
            long = "propagate-info-fields",
            help = "Additional INFO fields in the input BCF that shall be propagated to the output BCF \
            (not supported for variants connected via EVENT tags except breakends)."
        )]
        propagate_info_fields: Vec<String>,
        #[structopt(
            long = "strandedness",
            default_value = "opposite",
            possible_values = &ProtocolStrandedness::iter().map(|v| v.into()).collect_vec(),
            help = "Strandedness of sequencing protocol in case of paired-end (opposite strand as usual or same strand as with mate-pair sequencing.)"
        )]
        protocol_strandedness: ProtocolStrandedness,
        #[structopt(
            long = "indel-window",
            default_value = "64",
            help = "Number of bases to consider left and right of breakpoint when \
                    calculating read support. Currently implemented maximum \
                    value is 64."
        )]
        realignment_window: u64,
        #[structopt(
            long = "max-depth",
            default_value = "200",
            help = "Maximum number of observations to use for calling. If locus is exceeding this \
                    number, downsampling is performed."
        )]
        max_depth: usize,
        #[structopt(
            long = "omit-insert-size",
            help = "Do not consider insert size when calculating support for a variant. Use this flag when \
                    processing amplicon data, where indels do not impact the observed insert size"
        )]
        #[serde(default)]
        omit_insert_size: bool,
        #[structopt(
            long = "pairhmm-mode",
            possible_values = &["fast", "exact", "homopolymer"],
            default_value = "exact",
            help = "PairHMM computation mode (either fast, exact or homopolymer). Fast mode means that only the best \
                    alignment path is considered for probability calculation. In rare cases, this can lead \
                    to wrong results for single reads. Hence, we advice to not use it when \
                    discrete allele frequences are of interest (0.5, 1.0). For continuous \
                    allele frequencies, fast mode should cause almost no deviations from the \
                    exact results. Also, if per sample allele frequencies are irrelevant (e.g. \
                    in large cohorts), fast mode can be safely used. \
                    Note that fast and exact mode are not suitable for ONT data.\
                    The homopolymer mode should be used for ONT data; it is similar to the exact mode \
                    but considers homopolymer errors to be different from gaps."
        )]
        #[serde(default = "default_pairhmm_mode")]
        pairhmm_mode: String,
        #[structopt(
            long = "log-mode",
            possible_values = &["default", "each-record"],
            default_value = "default",
            help = "Specify how progress should be logged. By default, a record count will be printed. With 'each-record', \
            Varlociraptor will additionally print contig and position of each processed candidate variant record. This is \
            useful for debugging."
        )]
        #[serde(default = "default_log_mode")]
        log_mode: String,
        #[structopt(
            parse(from_os_str),
            long = "output-raw-observations",
            help = "Output raw observations to the given TSV file path. Attention, only use this for debugging when processing \
            a single variant. Otherwise it will cause a huge file and significant performance hits."
        )]
        #[serde(default)]
        output_raw_observations: Option<PathBuf>,
    },
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum PlotKind {
    #[structopt(
        name = "variant-calling-prior",
        about = "Plot variant calling prior given a scenario. Plot is printed to STDOUT in Vega-lite format.",
        usage = "varlociraptor plot variant-calling-prior --scenario scenario.yaml > plot.vl.json",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    VariantCallingPrior {
        #[structopt(
            parse(from_os_str),
            long = "scenario",
            required = true,
            help = "Variant calling scenario that configures the prior."
        )]
        scenario: PathBuf,
        #[structopt(
            long = "contig",
            required = true,
            help = "Contig to consider for ploidy information."
        )]
        contig: String,
        #[structopt(long = "sample", required = true, help = "Sample to plot.")]
        sample: String,
    },
    #[structopt(
        name = "scatter",
        about = "Plot variant allelic fraction scatter plot overlayed with a contour plot between two sample groups",
        usage = "varlociraptor plot scatter --sample-x sample1 \
        --sample-y sample2 sample3 < calls.bcf | vg2svg > scatter.svg",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Scatter {
        #[structopt(
            long = "sample-x",
            help = "Name of the first sample in the given VCF/BCF."
        )]
        sample_x: String,
        #[structopt(
            long = "sample-y",
            help = "Name(s) of the alternative sample(s) in the given VCF/BCF. Multiple samples can be given."
        )]
        sample_y: Vec<String>,
    },
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum EstimateKind {
    #[structopt(
        name = "alignment-properties",
        about = "Estimate properties like insert size, maximum softclip length, and the PCR homopolymer error model.",
        usage = "varlociraptor estimate alignment-properties reference.fasta --bam sample.bam > sample.alignment-properties.json",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    AlignmentProperties {
        #[structopt(
            parse(from_os_str),
            help = "FASTA file with reference genome. Has to be indexed with samtools faidx."
        )]
        reference: PathBuf,

        #[structopt(
            long,
            required = true,
            help = "BAM file with aligned reads from a single sample."
        )]
        bam: PathBuf,

        #[structopt(long, help = "Number of records to sample from the BAM file")]
        num_records: Option<usize>,

        #[structopt(
            long,
            help = "Estimated gap extension probabilities below this threshold will be set to zero.",
            default_value = "1e-10"
        )]
        epsilon_gap: f64,
    },
    #[structopt(
        name = "contamination",
        about = "Estimate contamination between samples (still experimental). Takes two samples, the one that \
        is presumably contaminated and the contaminating one, and calculates the posterior distribution of \
        the contamination fraction, which is printed into a TSV table. The most likely fraction comes at the top. \
        Uses a uniform prior distribution by default. Optionally, e.g. a pathologist's estimate can be provided, \
        such that the calculated posterior distribution can be sharpened by the model with the given prior \
        knowledge. For now, it is always advisable to check the visual output generated by --output-plot.",
        usage = "varlociraptor estimate contamination --sample sample-a.bcf --contaminant sample-b.bcf --output-plot contamination-qc-plot.vl.json > contamination.tsv"
    )]
    Contamination {
        #[structopt(long = "sample", help = "Presumably contaminated sample.")]
        sample: PathBuf,
        #[structopt(long = "contaminant", help = "Presumably contaminating sample.")]
        contaminant: PathBuf,
        #[structopt(
            long = "prior-estimate",
            help = "Prior estimate of the contamination (1-purity), e.g. \
            obtained by counting tumor cells in a histology sample. Has to be used \
            together with the number of cells considered for the estimate (--prior-considered-cells)."
        )]
        prior_estimate: Option<f64>,
        #[structopt(
            long = "prior-considered-cells",
            help = "Number of cells considered for the prior estimate (see --prior-estimate)."
        )]
        prior_considered_cells: Option<u32>,
        #[structopt(
            long = "output",
            help = "Output file; if not specified, output is printed to STDOUT."
        )]
        output: Option<PathBuf>,
        #[structopt(
            long = "output-plot",
            help = "Path to store vega-lite plot of contamination/purity."
        )]
        output_plot: Option<PathBuf>,
    },
    #[structopt(
        name = "mutational-burden",
        about = "Estimate mutational burden. Takes Varlociraptor calls (must be annotated \
                 with e.g. VEP but using ANN instead of CSQ) from STDIN, prints mutational burden estimate in Vega-lite JSON format to STDOUT. \
                 It can be converted to an image via vega-lite-cli (see conda package).",
        usage = "varlociraptor estimate mutational-burden --coding-genome-size 3e7 --events SOMATIC_TUMOR \
                 --sample tumor < calls.bcf | vg2svg > tmb.svg",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    MutationalBurden {
        #[structopt(long = "events", help = "Events to consider (e.g. SOMATIC_TUMOR).")]
        events: Vec<String>,
        #[structopt(
            long = "sample",
            help = "Name(s) of the sample(s) in the given VCF/BCF. Multiple samples can be given when using the multibar plot mode."
        )]
        sample: Vec<String>,
        #[structopt(
            long = "coding-genome-size",
            default_value = "3e7",
            help = "Size (in bases) of the covered coding genome."
        )]
        coding_genome_size: f64,
        #[structopt(
            long = "plot-mode",
            possible_values = &estimation::mutational_burden::PlotMode::iter().map(|v| v.into()).collect_vec(),
            help = "How to plot (as stratified curve, histogram or multi-sample barplot)."
        )]
        mode: estimation::mutational_burden::PlotMode,
        #[structopt(
            long = "vaf-cutoff",
            default_value = "0.2",
            help = "Minimal variant allelic fraction to consider for mutli-sample barplot"
        )]
        cutoff: f64,
    },
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum CallKind {
    #[structopt(
        name = "variants",
        about = "Call variants.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Variants {
        #[structopt(subcommand)]
        mode: VariantCallMode,
        #[structopt(
            long = "propagate-info-fields",
            help = "Additional INFO fields in the input BCF that shall be propagated to the output BCF."
        )]
        propagate_info_fields: Vec<String>,
        #[structopt(
            long = "omit-strand-bias",
            help = "Do not consider strand bias when calculating the probability of an artifact.\
                    Use this flag when processing (panel) sequencing data, where the wet-lab \
                    methodology leads to strand bias in the coverage of genuine variants."
        )]
        #[serde(default)]
        omit_strand_bias: bool,
        #[structopt(
            long = "omit-read-orientation-bias",
            help = "Do not consider read orientation bias when calculating the probability of an \
                    artifact."
        )]
        #[serde(default)]
        omit_read_orientation_bias: bool,
        #[structopt(
            long = "omit-read-position-bias",
            help = "Do not consider read position bias when calculating the probability of an \
                    artifact. Use this flag when processing (panel) sequencing data, where the \
                    wet-lab methodology leads to stacks of reads starting at the same position."
        )]
        #[serde(default)]
        omit_read_position_bias: bool,
        #[structopt(
            long = "omit-softclip-bias",
            help = "Do not consider softclip bias when calculating the probability of an \
                    artifact. Use this flag when processing (panel) sequencing data, where the \
                    wet-lab methodology leads to stacks of reads starting at the same position."
        )]
        #[serde(default)]
        omit_softclip_bias: bool,
        #[structopt(
            long = "omit-homopolymer-artifact-detection",
            help = "Do not perform PCR homopolymer artifact detection when calculating the probability of an \
                    artifact. If you are sure that your protocol did not use any PCR you should use this flag."
        )]
        #[serde(default)]
        omit_homopolymer_artifact_detection: bool,
        #[structopt(
            long = "omit-alt-locus-bias",
            help = "Do not consider alt locus bias when calculating the probability of an artifact. \
                   Use this flag when you have e.g. prior knowledge about all candidate variants being not \
                   caused by mapping artifacts."
        )]
        #[serde(default)]
        omit_alt_locus_bias: bool,
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
        #[structopt(
            long = "testcase-anonymous",
            help = "Anonymize any identifiers (via uuid4) and the sequences (by randomly permuting the alphabet) in the test case."
        )]
        testcase_anonymous: bool,
        #[structopt(
            long,
            short,
            help = "Output variant calls to given path (in BCF format). If omitted, prints calls to STDOUT."
        )]
        output: Option<PathBuf>,
        #[structopt(
            long = "log-mode",
            possible_values = &["default", "each-record"],
            default_value = "default",
            help = "Specify how progress should be logged. By default, a record count will be printed. With 'each-record', \
            Varlociraptor will additionally print contig and position of each processed candidate variant record. This is \
            useful for debugging."
        )]
        #[serde(default = "default_log_mode")]
        log_mode: String,
    },
    // #[structopt(
    //     name = "cnvs",
    //     about = "Call CNVs in tumor-normal sample pairs. This is experimental (do not use it yet).",
    //     setting = structopt::clap::AppSettings::ColoredHelp,
    // )]
    // CNVs {
    //     #[structopt(
    //         parse(from_os_str),
    //         long,
    //         help = "VCF/BCF file (generated by varlociraptor call-tumor-normal) to process \
    //                 (if omitted, read from STDIN)."
    //     )]
    //     calls: Option<PathBuf>,
    //     #[structopt(
    //         parse(from_os_str),
    //         long,
    //         help = "BCF file that shall contain the results (if omitted, write to STDOUT)."
    //     )]
    //     output: Option<PathBuf>,
    //     #[structopt(long, short = "p", help = "Tumor purity.")]
    //     purity: f64,
    //     #[structopt(
    //         long = "min-bayes-factor",
    //         default_value = "1.01",
    //         help = "Minimum bayes factor (> 1.0) between likelihoods of CNV and no CNV to consider. \
    //                 The higher this value, the fewer candidate CNVs will be investigated. \
    //                 Note that this can be usually left unchanged, because every CNV is provided \
    //                 with a posterior probability that can be used for filtering, e.g., via \
    //                 'varlociraptor filter-calls control-fdr'."
    //     )]
    //     min_bayes_factor: f64,
    //     #[structopt(
    //         long,
    //         default_value = "1000",
    //         help = "Maximum distance between supporting loci in a CNV."
    //     )]
    //     max_dist: u32,
    //     #[structopt(long, short = "t", help = "Number of threads to use.")]
    //     threads: usize,
    // },
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum VariantCallMode {
    #[structopt(
        name = "tumor-normal",
        about = "Call somatic and germline variants from a tumor-normal sample pair and a VCF/BCF with candidate variants.",
        usage = "varlociraptor call variants tumor-normal --purity 0.75 --tumor tumor.bcf --normal normal.bcf > calls.bcf",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    TumorNormal {
        #[structopt(
            parse(from_os_str),
            long = "tumor",
            required = true,
            help = "BCF file with varlociraptor preprocess results for the tumor sample."
        )]
        tumor_observations: PathBuf,
        #[structopt(
            parse(from_os_str),
            long = "normal",
            required = true,
            help = "BCF file with varlociraptor preprocess results for the normal sample."
        )]
        normal_observations: PathBuf,
        #[structopt(short, long, default_value = "1.0", help = "Purity of tumor sample.")]
        purity: f64,
    },
    #[structopt(
        name = "generic",
        about = "Call variants for a given scenario specified with the varlociraptor calling \
                 grammar and a VCF/BCF with candidate variants.",
        usage = "varlociraptor call variants --output calls.bcf generic --scenario scenario.yaml \
                 --obs relapse=relapse.bcf tumor=tumor.bcf normal=normal.bcf",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Generic {
        #[structopt(
            parse(from_os_str),
            long,
            required = true,
            help = "Scenario defined in the varlociraptor calling grammar."
        )]
        scenario: PathBuf,
        #[structopt(
            long = "obs",
            required = true,
            help = "BCF file with varlociraptor preprocess results for samples defined in the given scenario (given as samplename=path/to/calls.bcf). It is possible to omit a sample here (e.g. model tumor/normal in the scenario, but only call on the tumor sample when there is no normal sequenced). In that case, the resulting probabilities will be accordingly uncertain, because observations of the omitted sample are missing (which is equivalent to having no coverage in the sample)."
        )]
        sample_observations: Vec<String>,
    },
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum FilterMethod {
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
            possible_values = &VariantType::iter().map(|v| v.into()).collect_vec(),
            help = "Variant type to consider. When controlling global FDR (not using --local) this should \
            be used to control FDR for each type separately. Otherwise, less certain variant types will be \
            underrepresented."
        )]
        vartype: Option<VariantType>,
        #[structopt(long, help = "FDR to control for.")]
        fdr: f64,
        #[structopt(
            long = "local",
            help = "Control local FDR instead of global FDR. This means that for each record, the posterior \
            of the selected events has to be at least 1-fdr."
        )]
        local: bool,
        #[structopt(long, help = "Events to consider.")]
        events: Vec<String>,
        #[structopt(long, help = "Minimum indel length to consider.")]
        minlen: Option<u64>,
        #[structopt(long, help = "Maximum indel length to consider (exclusive).")]
        maxlen: Option<u64>,
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
            long,
            possible_values = &KassRaftery::iter().map(|v| v.into()).collect_vec(),
            help = "Kass-Raftery score to filter against."
        )]
        odds: KassRaftery,
        #[structopt(long, help = "Events to consider.")]
        events: Vec<String>,
    },
}

fn parse_key_values(values: &[String]) -> Option<PathMap> {
    let mut map = HashMap::new();
    for value in values {
        let kv = value.split_terminator('=').collect_vec();
        if kv.len() == 2 {
            map.insert(kv[0].to_owned(), PathBuf::from(kv[1]));
        } else {
            return None;
        }
    }
    Some(map)
}

impl Default for Varlociraptor {
    fn default() -> Self {
        Varlociraptor::from_iter(vec!["--help"])
    }
}

pub fn run(opt: Varlociraptor) -> Result<()> {
    let opt_clone = opt.clone();
    match opt {
        Varlociraptor::Preprocess { kind } => {
            match kind {
                PreprocessKind::Variants {
                    reference,
                    candidates,
                    bam,
                    report_fragment_ids,
                    alignment_properties,
                    output,
                    propagate_info_fields,
                    protocol_strandedness,
                    realignment_window,
                    max_depth,
                    omit_insert_size,
                    pairhmm_mode,
                    reference_buffer_size,
                    min_bam_refetch_distance,
                    log_mode,
                    output_raw_observations,
                } => {
                    // TODO: handle testcases
                    if realignment_window > (128 / 2) {
                        return Err(
                            structopt::clap::Error::with_description(
                                "Command-line option --indel-window requires a value <= 64 with the current implementation.", 
                                structopt::clap::ErrorKind::ValueValidation
                            ).into()
                        );
                    };

                    let mut reference_buffer = Arc::new(reference::Buffer::new(
                        fasta::IndexedReader::from_file(&reference)
                            .context("Unable to read genome reference.")?,
                        reference_buffer_size,
                    ));

                    let alignment_properties = est_or_load_alignment_properties(
                        &alignment_properties,
                        &bam,
                        omit_insert_size,
                        Arc::get_mut(&mut reference_buffer).unwrap(),
                        Some(crate::estimation::alignment_properties::NUM_FRAGMENTS),
                        0.,
                    )?;

                    let gap_params = alignment_properties.gap_params.clone();

                    let log_each_record = log_mode == "each-record";

                    let propagate_info_fields = propagate_info_fields
                        .iter()
                        .map(|s| s.as_bytes().to_owned())
                        .collect();

                    match pairhmm_mode.as_ref() {
                        "homopolymer" => {
                            let hop_params = alignment_properties.hop_params.clone();
                            let mut processor =
                                calling::variants::preprocessing::ObservationProcessor::builder()
                                    .report_fragment_ids(report_fragment_ids)
                                    .alignment_properties(alignment_properties)
                                    .protocol_strandedness(protocol_strandedness)
                                    .max_depth(max_depth)
                                    .inbam(bam)
                                    .min_bam_refetch_distance(min_bam_refetch_distance)
                                    .reference_buffer(Arc::clone(&reference_buffer))
                                    .haplotype_feature_index(HaplotypeFeatureIndex::new(
                                        &candidates,
                                    )?)
                                    .inbcf(candidates)
                                    .aux_info_fields(propagate_info_fields)
                                    .options(opt_clone)
                                    .outbcf(output)
                                    .raw_observation_output(output_raw_observations)
                                    .log_each_record(log_each_record)
                                    .realigner(realignment::HomopolyPairHMMRealigner::new(
                                        reference_buffer,
                                        gap_params,
                                        hop_params,
                                        realignment_window,
                                    ))
                                    .build();
                            processor.process()?;
                        }
                        "fast" => {
                            let mut processor =
                                calling::variants::preprocessing::ObservationProcessor::builder()
                                    .report_fragment_ids(report_fragment_ids)
                                    .alignment_properties(alignment_properties)
                                    .protocol_strandedness(protocol_strandedness)
                                    .max_depth(max_depth)
                                    .inbam(bam)
                                    .min_bam_refetch_distance(min_bam_refetch_distance)
                                    .reference_buffer(Arc::clone(&reference_buffer))
                                    .haplotype_feature_index(HaplotypeFeatureIndex::new(
                                        &candidates,
                                    )?)
                                    .inbcf(candidates)
                                    .aux_info_fields(propagate_info_fields)
                                    .options(opt_clone)
                                    .outbcf(output)
                                    .raw_observation_output(output_raw_observations)
                                    .log_each_record(log_each_record)
                                    .realigner(realignment::PathHMMRealigner::new(
                                        gap_params,
                                        realignment_window,
                                        reference_buffer,
                                    ))
                                    .build();
                            processor.process()?;
                        }
                        "exact" => {
                            let mut processor =
                                calling::variants::preprocessing::ObservationProcessor::builder()
                                    .report_fragment_ids(report_fragment_ids)
                                    .alignment_properties(alignment_properties)
                                    .protocol_strandedness(protocol_strandedness)
                                    .max_depth(max_depth)
                                    .inbam(bam)
                                    .min_bam_refetch_distance(min_bam_refetch_distance)
                                    .reference_buffer(Arc::clone(&reference_buffer))
                                    .haplotype_feature_index(HaplotypeFeatureIndex::new(
                                        &candidates,
                                    )?)
                                    .inbcf(candidates)
                                    .aux_info_fields(propagate_info_fields)
                                    .options(opt_clone)
                                    .outbcf(output)
                                    .raw_observation_output(output_raw_observations)
                                    .log_each_record(log_each_record)
                                    .realigner(realignment::PairHMMRealigner::new(
                                        reference_buffer,
                                        gap_params,
                                        realignment_window,
                                    ))
                                    .build();
                            processor.process()?;
                        }
                        _ => panic!("Unknown pairhmm mode '{}'", pairhmm_mode),
                    };
                }
            }
        }
        Varlociraptor::Call { kind } => {
            match kind {
                CallKind::Variants {
                    mode,
                    omit_strand_bias,
                    omit_read_orientation_bias,
                    omit_read_position_bias,
                    omit_softclip_bias,
                    omit_homopolymer_artifact_detection,
                    omit_alt_locus_bias,
                    testcase_locus,
                    testcase_prefix,
                    testcase_anonymous,
                    output,
                    log_mode,
                    propagate_info_fields,
                } => {
                    let testcase_builder = if let Some(testcase_locus) = testcase_locus {
                        if let Some(testcase_prefix) = testcase_prefix {
                            // TODO obtain sample information from input bcfs?
                            Some(
                                testcase::builder::TestcaseBuilder::default()
                                    .prefix(PathBuf::from(testcase_prefix))
                                    .anonymize(testcase_anonymous)
                                    .locus(&testcase_locus)?,
                            )
                        } else {
                            return Err(errors::Error::MissingPrefix.into());
                        }
                    } else {
                        None
                    };

                    let log_each_record = log_mode == "each-record";

                    match mode {
                        VariantCallMode::Generic {
                            scenario,
                            sample_observations,
                        } => {
                            if let Some(sample_observations) =
                                parse_key_values(&sample_observations)
                            {
                                if let Some(mut testcase_builder) = testcase_builder {
                                    for (i, (sample_name, obspath)) in
                                        sample_observations.iter().enumerate()
                                    {
                                        let options = calling::variants::preprocessing::read_preprocess_options(obspath)?;
                                        let preprocess_input = options.preprocess_input();
                                        testcase_builder = testcase_builder.register_sample(
                                            sample_name,
                                            preprocess_input.bam,
                                            options,
                                        )?;
                                        if i == 0 {
                                            testcase_builder =
                                                testcase_builder.candidates(obspath.to_owned());
                                            testcase_builder = testcase_builder
                                                .reference(preprocess_input.reference)?;
                                        }
                                    }

                                    let mut testcase = testcase_builder
                                        .scenario(Some(scenario))
                                        .mode(testcase::builder::Mode::Generic)
                                        .build()
                                        .unwrap();
                                    info!("Writing testcase.");
                                    testcase.write()?;
                                    return Ok(());
                                }

                                let scenario = grammar::Scenario::from_path(scenario)?;

                                call_generic(
                                    scenario,
                                    sample_observations,
                                    omit_strand_bias,
                                    omit_read_orientation_bias,
                                    omit_read_position_bias,
                                    omit_softclip_bias,
                                    omit_homopolymer_artifact_detection,
                                    omit_alt_locus_bias,
                                    output,
                                    log_each_record,
                                    CallWriter::new(),
                                    DefaultCandidateFilter::new(),
                                    propagate_info_fields,
                                )?;
                            } else {
                                return Err(errors::Error::InvalidObservationsSpec.into());
                            }
                        }
                        VariantCallMode::TumorNormal {
                            tumor_observations,
                            normal_observations,
                            purity,
                        } => {
                            if let Some(testcase_builder) = testcase_builder {
                                let tumor_options =
                                    calling::variants::preprocessing::read_preprocess_options(
                                        &tumor_observations,
                                    )?;
                                let normal_options =
                                    calling::variants::preprocessing::read_preprocess_options(
                                        &normal_observations,
                                    )?;
                                let mut testcase = testcase_builder
                                    .candidates(tumor_observations)
                                    .reference(tumor_options.preprocess_input().reference)?
                                    .register_sample(
                                        "tumor",
                                        tumor_options.preprocess_input().bam,
                                        tumor_options,
                                    )?
                                    .register_sample(
                                        "normal",
                                        normal_options.preprocess_input().bam,
                                        normal_options,
                                    )?
                                    .scenario(None)
                                    .mode(testcase::builder::Mode::TumorNormal)
                                    .purity(Some(purity))
                                    .build()
                                    .unwrap();

                                testcase.write()?;
                                return Ok(());
                            }

                            let scenario = grammar::Scenario::try_from(
                                format!(
                                    r#"
                            samples:
                              tumor:
                                resolution: 0.01
                                contamination:
                                  by: normal
                                  fraction: {impurity}
                                universe: "[0.0,1.0]"
                              normal:
                                resolution: 0.1
                                universe: "[0.0,0.5[ | 0.5 | 1.0"
                            events:
                              somatic_tumor:  "tumor:]0.0,1.0] & normal:0.0"
                              somatic_normal: "tumor:]0.0,1.0] & normal:]0.0,0.5["
                              germline_het:   "tumor:]0.0,1.0] & normal:0.5"
                              germline_hom:   "tumor:]0.0,1.0] & normal:1.0"
                            "#,
                                    impurity = 1.0 - purity
                                )
                                .as_str(),
                            )?;

                            let mut observations = PathMap::default();
                            observations.insert("tumor".to_owned(), tumor_observations);
                            observations.insert("normal".to_owned(), normal_observations);

                            call_generic(
                                scenario,
                                observations,
                                omit_strand_bias,
                                omit_read_orientation_bias,
                                omit_read_position_bias,
                                omit_softclip_bias,
                                omit_homopolymer_artifact_detection,
                                omit_alt_locus_bias,
                                output,
                                log_each_record,
                                CallWriter::new(),
                                DefaultCandidateFilter::new(),
                                propagate_info_fields,
                            )?;
                        }
                    }
                } // CallKind::CNVs {
                  //     calls,
                  //     output,
                  //     min_bayes_factor,
                  //     threads,
                  //     purity,
                  //     max_dist,
                  // } => {
                  //     rayon::ThreadPoolBuilder::new()
                  //         .num_threads(threads)
                  //         .build_global()?;

                  //     if min_bayes_factor <= 1.0 {
                  //         Err(errors::Error::InvalidMinBayesFactor)?
                  //     }

                  //     let mut caller = calling::cnvs::CallerBuilder::default()
                  //         .bcfs(calls.as_ref(), output.as_ref())?
                  //         .min_bayes_factor(min_bayes_factor)
                  //         .purity(purity)
                  //         .max_dist(max_dist)
                  //         .build()
                  //         .unwrap();
                  //     caller.call()?;
                  // }
            }
        }
        Varlociraptor::FilterCalls { method } => match method {
            FilterMethod::ControlFDR {
                calls,
                events,
                fdr,
                local,
                vartype,
                minlen,
                maxlen,
            } => {
                let events = events
                    .into_iter()
                    .map(|event| SimpleEvent { name: event })
                    .collect_vec();
                let vartype = match (vartype, minlen, maxlen) {
                    (Some(VariantType::Insertion(None)), Some(minlen), Some(maxlen)) => {
                        Some(VariantType::Insertion(Some(minlen..maxlen)))
                    }
                    (Some(VariantType::Deletion(None)), Some(minlen), Some(maxlen)) => {
                        Some(VariantType::Deletion(Some(minlen..maxlen)))
                    }
                    (vartype, _, _) => vartype,
                };

                filtration::fdr::control_fdr::<_, &PathBuf, &str>(
                    &calls,
                    None,
                    &events,
                    vartype.as_ref(),
                    LogProb::from(Prob::checked(fdr)?),
                    local,
                )?;
            }
            FilterMethod::PosteriorOdds { ref events, odds } => {
                let events = events
                    .iter()
                    .map(|event| SimpleEvent {
                        name: event.to_owned(),
                    })
                    .collect_vec();

                filtration::posterior_odds::filter_by_odds::<_, &PathBuf, &PathBuf>(
                    None, None, &events, odds,
                )?;
            }
        },
        Varlociraptor::DecodePHRED => {
            conversion::decode_phred::decode_phred()?;
        }
        Varlociraptor::Genotype => {
            conversion::genotype::genotype()?;
        }
        Varlociraptor::Estimate { kind } => match kind {
            EstimateKind::Contamination {
                sample,
                contaminant,
                prior_estimate,
                prior_considered_cells,
                output,
                output_plot,
            } => {
                let prior_estimate = match (prior_estimate, prior_considered_cells) {
                    (Some(p), Some(n)) if n > 0 => Some(
                        estimation::contamination::PriorEstimate::new(AlleleFreq(p), n),
                    ),
                    (None, None) => None,
                    _ => bail!(errors::Error::InvalidPriorContaminationEstimate),
                };
                estimation::contamination::estimate_contamination(
                    sample,
                    contaminant,
                    output,
                    output_plot,
                    prior_estimate,
                )?;
            }
            EstimateKind::MutationalBurden {
                events,
                sample,
                coding_genome_size,
                mode,
                cutoff,
            } => estimation::mutational_burden::collect_estimates(
                &events,
                &sample,
                coding_genome_size as u64,
                mode,
                cutoff as f64,
            )?,
            EstimateKind::AlignmentProperties {
                reference,
                bam,
                num_records,
                epsilon_gap,
            } => {
                let mut reference_buffer =
                    reference::Buffer::new(fasta::IndexedReader::from_file(&reference)?, 1);
                let alignment_properties = estimate_alignment_properties(
                    bam,
                    false,
                    &mut reference_buffer,
                    num_records,
                    epsilon_gap,
                )?;
                println!("{}", serde_json::to_string_pretty(&alignment_properties)?);
            }
        },
        Varlociraptor::Plot { kind } => match kind {
            PlotKind::VariantCallingPrior {
                scenario,
                contig,
                sample,
            } => {
                let scenario = grammar::Scenario::from_path(scenario)?;
                let sample_infos = SampleInfos::try_from(&scenario)?;

                let mut universes = scenario.sample_info();
                let mut ploidies = scenario.sample_info();
                for (sample_name, sample) in scenario.samples().iter() {
                    universes = universes.push(
                        sample_name,
                        sample.contig_universe(&contig, scenario.species())?,
                    );
                    ploidies = ploidies.push(
                        sample_name,
                        sample.contig_ploidy(&contig, scenario.species())?,
                    );
                }
                let universes = universes.build();
                let ploidies = ploidies.build();

                let prior = Prior::builder()
                    .variant_type_fractions(scenario.variant_type_fractions())
                    .ploidies(Some(ploidies))
                    .universe(Some(universes))
                    .uniform(sample_infos.uniform_prior().clone())
                    .germline_mutation_rate(sample_infos.germline_mutation_rates().clone())
                    .somatic_effective_mutation_rate(
                        sample_infos.somatic_effective_mutation_rates().clone(),
                    )
                    .inheritance(sample_infos.inheritance().clone())
                    .heterozygosity(scenario.species().as_ref().and_then(|species| {
                        species.heterozygosity().map(|het| LogProb::from(Prob(het)))
                    }))
                    .variant_type(Some(VariantType::Snv))
                    .build();
                prior.check()?;

                prior.plot(&sample, sample_infos.names())?;
            }
            PlotKind::Scatter { sample_x, sample_y } => {
                estimation::sample_variants::vaf_scatter(&sample_x, &sample_y)?
            }
        },
    }
    Ok(())
}

pub(crate) fn est_or_load_alignment_properties(
    alignment_properties_file: &Option<impl AsRef<Path>>,
    bam_file: impl AsRef<Path>,
    omit_insert_size: bool,
    reference_buffer: &mut reference::Buffer,
    num_records: Option<usize>,
    epsilon_gap: f64,
) -> Result<AlignmentProperties> {
    if let Some(alignment_properties_file) = alignment_properties_file {
        Ok(serde_json::from_reader(File::open(
            alignment_properties_file,
        )?)?)
    } else {
        estimate_alignment_properties(
            bam_file,
            omit_insert_size,
            reference_buffer,
            num_records,
            epsilon_gap,
        )
    }
}
