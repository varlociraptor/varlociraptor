// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::error::Error;
use std::fs::File;
use std::path::{Path, PathBuf};

use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::bayesian::model::Model;
use bio::stats::{LogProb, Prob};
use itertools::Itertools;
use rayon;
use rust_htslib::bam;
use structopt::StructOpt;

use crate::call::CallerBuilder;
use crate::call_cnvs;
use crate::conversion;
use crate::errors;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::filtration;
use crate::model::modes::common::FlatPrior;
use crate::model::modes::tumor::{TumorNormalLikelihood, TumorNormalPair, TumorNormalPosterior};
use crate::model::sample::{estimate_alignment_properties, SampleBuilder};
use crate::model::{ContinuousAlleleFreqs, VariantType};
use crate::testcase::TestcaseBuilder;
use crate::SimpleEvent;

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
#[structopt(
    name = "varlociraptor",
    about = "A caller for SNVs and indels in tumor-normal pairs."
)]
#[structopt(raw(setting = "structopt::clap::AppSettings::ColoredHelp"))]
pub enum Varlociraptor {
    #[structopt(
        name = "call-tumor-normal",
        about = "Call somatic and germline variants from a tumor-normal sample pair and a VCF/BCF with candidate variants."
    )]
    #[structopt(raw(setting = "structopt::clap::AppSettings::ColoredHelp"))]
    CallTumorNormal {
        #[structopt(parse(from_os_str), help = "BAM file with reads from tumor sample.")]
        tumor: PathBuf,
        #[structopt(parse(from_os_str), help = "BAM file with reads from normal sample.")]
        normal: PathBuf,
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
        #[structopt(short, long, default_value = "1.0", help = "Purity of tumor sample.")]
        purity: f64,
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
            default_value = "100",
            help = "Number of bases to consider left and right of indel breakpoint when \
                    calculating read support. This number should not be too large in order to \
                    avoid biases caused by other close variants."
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
        #[structopt(
            parse(from_os_str),
            long,
            help = "Alignment properties JSON file for tumor sample. If not provided, properties \
                    will be estimated from the given BAM file."
        )]
        tumor_alignment_properties: Option<PathBuf>,
        #[structopt(
            parse(from_os_str),
            long,
            help = "Alignment properties JSON file for normal sample. If not provided, properties \
                    will be estimated from the given BAM file."
        )]
        normal_alignment_properties: Option<PathBuf>,
    },
    #[structopt(name = "call-cnvs", about = "Call CNVs in tumor-normal sample pairs.")]
    #[structopt(raw(setting = "structopt::clap::AppSettings::ColoredHelp"))]
    CallCNVs {
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
            long,
            default_value = "0.00001",
            help = "Prior probability for CNV. This can be any small enough value."
        )]
        prior: f64,
        #[structopt(long, short = "t", help = "Number of threads to use.")]
        threads: usize,
    },
    #[structopt(
        name = "filter-calls",
        about = "Filter calls by either controlling the false discovery rate (FDR) at given level, or by posterior odds against the given events."
    )]
    #[structopt(raw(setting = "structopt::clap::AppSettings::ColoredHelp"))]
    FilterCalls {
        #[structopt(subcommand)]
        method: FilterMethod,
    },
    #[structopt(
        name = "decode-phred",
        about = "Decode PHRED-scaled values to human readable probabilities."
    )]
    DecodePHRED,
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum FilterMethod {
    #[structopt(name = "control-fdr")]
    #[structopt(raw(setting = "structopt::clap::AppSettings::ColoredHelp"))]
    ControlFDR {
        #[structopt(parse(from_os_str), help = "BCF file with varlociraptor calls.")]
        calls: PathBuf,
        #[structopt(
            long = "var",
            raw(
                possible_values = "{use strum::IntoEnumIterator; &VariantType::iter().map(|v| v.into()).collect_vec()}"
            ),
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
    #[structopt(name = "posterior-odds")]
    #[structopt(raw(setting = "structopt::clap::AppSettings::ColoredHelp"))]
    PosteriorOdds {
        #[structopt(
            raw(
                possible_values = "{use strum::IntoEnumIterator; &KassRaftery::iter().map(|v| v.into()).collect_vec()}"
            ),
            help = "Kass-Raftery score to filter against."
        )]
        odds: KassRaftery,
        #[structopt(long, help = "Events to consider.")]
        events: Vec<String>,
    },
}

impl Default for Varlociraptor {
    fn default() -> Self {
        Varlociraptor::from_iter(vec!["--help"])
    }
}

pub fn run(opt: Varlociraptor) -> Result<(), Box<Error>> {
    match opt {
        Varlociraptor::CallTumorNormal {
            ref tumor,
            ref normal,
            spurious_ins_rate,
            spurious_del_rate,
            spurious_insext_rate,
            spurious_delext_rate,
            indel_window,
            omit_snvs,
            omit_indels,
            max_indel_len,
            max_depth,
            ref reference,
            ref candidates,
            purity,
            ref output,
            ref testcase_locus,
            ref testcase_prefix,
            ref tumor_alignment_properties,
            ref normal_alignment_properties,
        } => {
            if let Some(testcase_locus) = testcase_locus {
                if let Some(testcase_prefix) = testcase_prefix {
                    if let Some(candidates) = candidates {
                        // just write a testcase and quit
                        let mut testcase = TestcaseBuilder::default()
                            .prefix(PathBuf::from(testcase_prefix))
                            .options(opt.clone())
                            .locus(testcase_locus)?
                            .reference(reference)?
                            .candidates(candidates)?
                            .register_bam("tumor", tumor)
                            .register_bam("normal", normal)
                            .build()?;
                        testcase.write()?;
                        return Ok(());
                    } else {
                        Err(errors::TestcaseError::MissingCandidates)?;
                    }
                } else {
                    Err(errors::TestcaseError::MissingPrefix)?;
                }
            }

            let tumor_alignment_properties =
                est_or_load_alignment_properites(tumor_alignment_properties, tumor)?;
            let normal_alignment_properties =
                est_or_load_alignment_properites(normal_alignment_properties, normal)?;
            info!("Estimated alignment properties:");
            info!("{:?}", tumor_alignment_properties);
            info!("{:?}", normal_alignment_properties);

            let tumor_bam = bam::IndexedReader::from_path(tumor)?;
            let normal_bam = bam::IndexedReader::from_path(normal)?;

            let spurious_ins_rate = Prob::checked(spurious_ins_rate)?;
            let spurious_del_rate = Prob::checked(spurious_del_rate)?;
            let spurious_insext_rate = Prob::checked(spurious_insext_rate)?;
            let spurious_delext_rate = Prob::checked(spurious_delext_rate)?;

            let sample_builder = || {
                SampleBuilder::default()
                    .error_probs(
                        spurious_ins_rate,
                        spurious_del_rate,
                        spurious_insext_rate,
                        spurious_delext_rate,
                        indel_window as u32,
                    )
                    .max_depth(max_depth)
            };

            let tumor_sample = sample_builder()
                .name("tumor".to_owned())
                .alignments(tumor_bam, tumor_alignment_properties)
                .build()?;
            let normal_sample = sample_builder()
                .name("normal".to_owned())
                .alignments(normal_bam, normal_alignment_properties)
                .build()?;

            let model = Model::new(
                TumorNormalLikelihood::new(purity),
                FlatPrior::new(),
                TumorNormalPosterior::new(),
            );

            let mut caller = CallerBuilder::default()
                .samples(vec![tumor_sample, normal_sample])
                .reference(reference)?
                .inbcf(candidates.as_ref())?
                .model(model)
                .omit_snvs(omit_snvs)
                .omit_indels(omit_indels)
                .max_indel_len(max_indel_len)
                .event(
                    "germline_het",
                    TumorNormalPair {
                        tumor: ContinuousAlleleFreqs::inclusive(0.0..1.0),
                        normal: ContinuousAlleleFreqs::singleton(0.5),
                    },
                )
                .event(
                    "germline_hom",
                    TumorNormalPair {
                        tumor: ContinuousAlleleFreqs::inclusive(0.0..1.0),
                        normal: ContinuousAlleleFreqs::singleton(1.0),
                    },
                )
                .event(
                    "somatic_tumor",
                    TumorNormalPair {
                        tumor: ContinuousAlleleFreqs::left_exclusive(0.0..1.0).min_observations(2),
                        normal: ContinuousAlleleFreqs::absent(),
                    },
                )
                .event(
                    "somatic_normal",
                    TumorNormalPair {
                        tumor: ContinuousAlleleFreqs::left_exclusive(0.0..1.0),
                        normal: ContinuousAlleleFreqs::exclusive(0.0..0.5).min_observations(2),
                    },
                )
                .event(
                    "absent",
                    TumorNormalPair {
                        tumor: ContinuousAlleleFreqs::absent(),
                        normal: ContinuousAlleleFreqs::absent(),
                    },
                )
                .outbcf(output.as_ref())?
                .build()?;

            caller.call()?;
        }
        Varlociraptor::CallCNVs {
            ref calls,
            ref output,
            prior,
            threads,
            purity,
        } => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()?;
            let prior = Prob::checked(prior)?;
            let mut caller = call_cnvs::CallerBuilder::default()
                .bcfs(calls.as_ref(), output.as_ref())?
                .prior(LogProb::from(prior))
                .purity(purity)
                .build()?;
            caller.call()?;
        }
        Varlociraptor::FilterCalls { method } => match method {
            FilterMethod::ControlFDR {
                ref calls,
                ref events,
                fdr,
                ref vartype,
                minlen,
                maxlen,
            } => {
                let events = events
                    .into_iter()
                    .map(|event| SimpleEvent {
                        name: event.to_owned(),
                    })
                    .collect_vec();
                let vartype = match (vartype, minlen, maxlen) {
                    (&VariantType::Insertion(None), Some(minlen), Some(maxlen)) => {
                        VariantType::Insertion(Some(minlen..maxlen))
                    }
                    (&VariantType::Deletion(None), Some(minlen), Some(maxlen)) => {
                        VariantType::Deletion(Some(minlen..maxlen))
                    }
                    (vartype @ _, _, _) => vartype.clone(),
                };

                filtration::fdr::control_fdr::<_, &PathBuf, &str>(
                    calls,
                    None,
                    &events,
                    &vartype,
                    LogProb::from(Prob::checked(fdr)?),
                )?;
            }
            FilterMethod::PosteriorOdds { ref events, odds } => {
                let events = events
                    .into_iter()
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
    }
    Ok(())
}

pub fn est_or_load_alignment_properites(
    alignment_properties_file: &Option<impl AsRef<Path>>,
    bam_file: impl AsRef<Path>,
) -> Result<AlignmentProperties, Box<Error>> {
    if let Some(alignment_properties_file) = alignment_properties_file {
        Ok(serde_json::from_reader(File::open(
            alignment_properties_file,
        )?)?)
    } else {
        estimate_alignment_properties(bam_file)
    }
}
