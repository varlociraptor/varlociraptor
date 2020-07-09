// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::HashMap;
use std::convert::{From, TryFrom};
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use bio::io::fasta;
use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::{LogProb, Prob};
use itertools::Itertools;
use rust_htslib::{bam, bcf};
use structopt::StructOpt;
use strum::IntoEnumIterator;

use crate::calling;
use crate::conversion;
use crate::errors;
use crate::estimation;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::filtration;
use crate::grammar;
use crate::testcase;
use crate::variants::evidence::realignment::pairhmm::GapParams;
use crate::variants::model::modes::generic::{FlatPrior, GenericModelBuilder};
use crate::variants::model::{Contamination, VariantType};
use crate::variants::sample::{estimate_alignment_properties, ProtocolStrandedness, SampleBuilder};
use crate::variants::types::breakends::BreakendIndex;
use crate::SimpleEvent;

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
#[structopt(
    name = "varlociraptor",
    about = "A caller for SNVs and indels in tumor-normal pairs.",
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

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum PreprocessKind {
    #[structopt(
        name = "variants",
        about = "Preprocess given variants by obtaining internal observations (allele likelihoods, strand information, ...)\
                 for each fragment. \
                 The obtained observations are printed to STDOUT in BCF format. Note that the resulting BCFs \
                 will be very large and are only intended for internal use (e.g. for piping into 'varlociraptor \
                 call variants generic').",
        usage = "varlociraptor preprocess variants reference.fasta --candidates candidates.bcf --bam sample.bam --output sample.observations.bcf",
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
            help = "VCF/BCF file to process (if omitted, read from STDIN)."
        )]
        candidates: PathBuf,
        #[structopt(
            long,
            required = true,
            help = "BAM file with aligned reads from a single sample."
        )]
        bam: PathBuf,
        #[structopt(
            long = "alignment-properties",
            help = "Alignment properties JSON file for sample. If not provided, properties \
                    will be estimated from the given BAM file."
        )]
        alignment_properties: Option<PathBuf>,
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
    },
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum EstimateKind {
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
        #[structopt(
            long = "plot-mode",
            possible_values = &estimation::tumor_mutational_burden::PlotMode::iter().map(|v| v.into()).collect_vec(),
            help = "How to plot (as stratified curve or as histogram)."
        )]
        mode: estimation::tumor_mutational_burden::PlotMode,
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
            long,
            short,
            help = "Output variant calls to given path (in BCF format). If omitted, prints calls to STDOUT."
        )]
        output: Option<PathBuf>,
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
        usage = "varlociraptor call variants --output calls.bcf tumor-normal --purity 0.75 --tumor tumor.bcf --normal normal.bcf",
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
            help = "BCF file with varlociraptor preprocess results for each sample defined in the given scenario (given as samplename=path/to/calls.bcf)."
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
            help = "Variant type to consider."
        )]
        vartype: VariantType,
        #[structopt(long, help = "FDR to control for.")]
        fdr: f64,
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

pub(crate) type PathMap = HashMap<String, PathBuf>;

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
                    alignment_properties,
                    output,
                    spurious_ins_rate,
                    spurious_del_rate,
                    spurious_insext_rate,
                    spurious_delext_rate,
                    protocol_strandedness,
                    realignment_window,
                    max_depth,
                    omit_insert_size,
                } => {
                    // TODO: handle testcases

                    let spurious_ins_rate = Prob::checked(spurious_ins_rate)?;
                    let spurious_del_rate = Prob::checked(spurious_del_rate)?;
                    let spurious_insext_rate = Prob::checked(spurious_insext_rate)?;
                    let spurious_delext_rate = Prob::checked(spurious_delext_rate)?;
                    if realignment_window > (128 / 2) {
                        return Err(
                            structopt::clap::Error::with_description(
                                "Command-line option --indel-window requires a value <= 64 with the current implementation.", 
                                structopt::clap::ErrorKind::ValueValidation
                            ).into()
                        );
                    };

                    // If we omit the insert size information for calculating the evidence, we can savely allow hardclips here.
                    let allow_hardclips = omit_insert_size;
                    let alignment_properties = est_or_load_alignment_properites(
                        &alignment_properties,
                        &bam,
                        allow_hardclips,
                    )?;

                    let bam_reader = bam::IndexedReader::from_path(bam)
                        .context("Unable to read BAM/CRAM file.")?;

                    let gap_params = GapParams {
                        prob_insertion_artifact: LogProb::from(spurious_ins_rate),
                        prob_deletion_artifact: LogProb::from(spurious_del_rate),
                        prob_insertion_extend_artifact: LogProb::from(spurious_insext_rate),
                        prob_deletion_extend_artifact: LogProb::from(spurious_delext_rate),
                    };

                    let sample = SampleBuilder::default()
                        .max_depth(max_depth)
                        .protocol_strandedness(protocol_strandedness)
                        .alignments(bam_reader, alignment_properties)
                        .use_fragment_evidence(!omit_insert_size)
                        .build()
                        .unwrap();

                    let mut processor =
                        calling::variants::preprocessing::ObservationProcessorBuilder::default()
                            .sample(sample)
                            .reference(
                                fasta::IndexedReader::from_file(&reference)
                                    .context("Unable to read genome reference.")?,
                            )
                            .realignment(gap_params, realignment_window)
                            .breakend_index(BreakendIndex::new(&candidates)?)
                            .inbcf(candidates)
                            .outbcf(output, &opt_clone)?
                            .build()
                            .unwrap();

                    processor.process()?
                }
            }
        }
        Varlociraptor::Call { kind } => {
            match kind {
                CallKind::Variants {
                    mode,
                    testcase_locus,
                    testcase_prefix,
                    output,
                } => {
                    let testcase_builder = if let Some(testcase_locus) = testcase_locus {
                        if let Some(testcase_prefix) = testcase_prefix {
                            // TODO obtain sample information from input bcfs?
                            Some(
                                testcase::TestcaseBuilder::default()
                                    .prefix(PathBuf::from(testcase_prefix))
                                    .locus(&testcase_locus)?,
                            )
                        } else {
                            return Err(errors::Error::MissingPrefix.into());
                        }
                    } else {
                        None
                    };

                    let call_generic = |scenario: grammar::Scenario,
                                        observations: PathMap|
                     -> Result<()> {
                        let mut contaminations = scenario.sample_info();
                        let mut resolutions = scenario.sample_info();
                        let mut sample_names = scenario.sample_info();
                        let mut sample_observations = scenario.sample_info();

                        // parse samples
                        for (sample_name, sample) in scenario.samples().iter() {
                            let contamination = if let Some(contamination) = sample.contamination()
                            {
                                let contaminant = scenario.idx(contamination.by()).ok_or(
                                    errors::Error::InvalidContaminationSampleName {
                                        name: sample_name.to_owned(),
                                    },
                                )?;
                                Some(Contamination {
                                    by: contaminant,
                                    fraction: *contamination.fraction(),
                                })
                            } else {
                                None
                            };
                            contaminations = contaminations.push(sample_name, contamination);
                            resolutions = resolutions.push(sample_name, *sample.resolution());

                            if let Some(obs) = observations.get(sample_name) {
                                sample_observations = sample_observations
                                    .push(sample_name, Some(bcf::Reader::from_path(obs)?));
                            } else {
                                sample_observations = sample_observations.push(sample_name, None);
                            }

                            sample_names = sample_names.push(sample_name, sample_name.to_owned());
                        }

                        let sample_names = sample_names.build();

                        for obs_sample_name in observations.keys() {
                            if !sample_names.as_slice().contains(obs_sample_name) {
                                return Err(errors::Error::InvalidObservationSampleName {
                                    name: obs_sample_name.to_owned(),
                                }
                                .into());
                            }
                        }

                        let model = GenericModelBuilder::default()
                            // TODO allow to define prior in the grammar
                            .prior(FlatPrior::new())
                            .contaminations(contaminations.build())
                            .resolutions(resolutions.build())
                            .build()
                            .unwrap();

                        // setup caller
                        let mut caller = calling::variants::CallerBuilder::default()
                            .samplenames(sample_names)
                            .observations(sample_observations.build())
                            .scenario(scenario)
                            .model(model)
                            .outbcf(output.as_ref())?
                            .build()
                            .unwrap();

                        // call
                        caller.call()?;

                        Ok(())
                    };

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
                                            &sample_name,
                                            preprocess_input.bam,
                                            &options,
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
                                        .mode(testcase::Mode::Generic)
                                        .build()
                                        .unwrap();
                                    info!("Writing testcase.");
                                    testcase.write()?;
                                    return Ok(());
                                }

                                let mut scenario_content = String::new();
                                File::open(scenario)?.read_to_string(&mut scenario_content)?;

                                let scenario: grammar::Scenario =
                                    serde_yaml::from_str(&scenario_content)?;

                                call_generic(scenario, sample_observations)?;
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
                                        &tumor_options,
                                    )?
                                    .register_sample(
                                        "normal",
                                        normal_options.preprocess_input().bam,
                                        &normal_options,
                                    )?
                                    .scenario(None)
                                    .mode(testcase::Mode::TumorNormal)
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
                                resolution: 100
                                contamination:
                                  by: normal
                                  fraction: {impurity}
                                universe: "[0.0,1.0]"
                              normal:
                                resolution: 5
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

                            call_generic(scenario, observations)?;
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
                vartype,
                minlen,
                maxlen,
            } => {
                let events = events
                    .into_iter()
                    .map(|event| SimpleEvent { name: event })
                    .collect_vec();
                let vartype = match (vartype, minlen, maxlen) {
                    (VariantType::Insertion(None), Some(minlen), Some(maxlen)) => {
                        VariantType::Insertion(Some(minlen..maxlen))
                    }
                    (VariantType::Deletion(None), Some(minlen), Some(maxlen)) => {
                        VariantType::Deletion(Some(minlen..maxlen))
                    }
                    (vartype, _, _) => vartype,
                };

                filtration::fdr::control_fdr::<_, &PathBuf, &str>(
                    &calls,
                    None,
                    &events,
                    &vartype,
                    LogProb::from(Prob::checked(fdr)?),
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
        Varlociraptor::Estimate { kind } => match kind {
            EstimateKind::TMB {
                somatic_tumor_events,
                tumor_sample,
                coding_genome_size,
                mode,
            } => estimation::tumor_mutational_burden::estimate(
                &somatic_tumor_events,
                &tumor_sample,
                coding_genome_size as u64,
                mode,
            )?,
        },
    }
    Ok(())
}

pub(crate) fn est_or_load_alignment_properites(
    alignment_properties_file: &Option<impl AsRef<Path>>,
    bam_file: impl AsRef<Path>,
    omit_insert_size: bool,
) -> Result<AlignmentProperties> {
    if let Some(alignment_properties_file) = alignment_properties_file {
        Ok(serde_json::from_reader(File::open(
            alignment_properties_file,
        )?)?)
    } else {
        estimate_alignment_properties(bam_file, omit_insert_size)
    }
}
