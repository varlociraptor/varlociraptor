// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::error::Error;
use std::path::PathBuf;

use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::bayesian::model::Model;
use bio::stats::Prob;
use rust_htslib::bam;

use itertools::Itertools;
use structopt::StructOpt;

use varlociraptor::call::CallerBuilder;
use varlociraptor::model::modes::common::FlatPrior;
use varlociraptor::model::modes::tumor::{
    TumorNormalLikelihood, TumorNormalPair, TumorNormalPosterior,
};
use varlociraptor::model::sample::{estimate_alignment_properties, SampleBuilder};
use varlociraptor::model::ContinuousAlleleFreqs;

#[derive(Debug, StructOpt)]
#[structopt(
    name = "varlociraptor",
    about = "A caller for SNVs and indels in tumor-normal pairs."
)]
#[structopt(raw(setting = "structopt::clap::AppSettings::ColoredHelp"))]
enum Varlociraptor {
    #[structopt(
        name = "call-tumor-normal",
        about = "Call somatic and germline variants from a tumor-normal sample pair and a VCF/BCF with candidate variants."
    )]
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
            help = "VCF/BCF file to process (if omitted, read from STDIN)."
        )]
        candidates: Option<PathBuf>,
        #[structopt(
            parse(from_os_str),
            help = "BCF file that shall contain the results (if omitted, write to STDOUT)."
        )]
        output: Option<PathBuf>,
        #[structopt(short, long, default_value = "1.0", help = "Purity of tumor sample.")]
        purity: f64,
        #[structopt(
            long,
            default_value = "2.8e-6",
            help = "Rate of spuriously inserted bases by the sequencer (Illumina: 2.8e-6, see Schirmer et al. BMC Bioinformatics 2016)."
        )]
        spurious_ins_rate: f64,
        #[structopt(
            long,
            default_value = "5.1e-6",
            help = "Rate of spuriosly deleted bases by the sequencer (Illumina: 5.1e-6, see Schirmer et al. BMC Bioinformatics 2016)."
        )]
        spurious_del_rate: f64,
        #[structopt(
            long,
            default_value = "0.0",
            help = "Extension rate of spurious insertions by the sequencer (Illumina: 0.0, see Schirmer et al. BMC Bioinformatics 2016)"
        )]
        spurious_insext_rate: f64,
        #[structopt(
            long,
            default_value = "0.0",
            help = "Extension rate of spurious deletions by the sequencer (Illumina: 0.0, see Schirmer et al. BMC Bioinformatics 2016)"
        )]
        spurious_delext_rate: f64,
        #[structopt(long, help = "Don't call SNVs.")]
        omit_snvs: bool,
        #[structopt(long, help = "Don't call Indels.")]
        omit_indels: bool,
        #[structopt(
            parse(from_os_str),
            help = "Optional path where read observations shall be written to. The resulting file contains a line for each observation with tab-separated values."
        )]
        observations: Option<PathBuf>,
        #[structopt(
            long,
            default_value = "1000",
            help = "Omit longer indels when calling."
        )]
        max_indel_len: u32,
        #[structopt(
            help = "Assume that the END tag is exclusive (i.e. it points to the position after the variant). This is needed, e.g., for DELLY."
        )]
        exclusive_end: bool,
        #[structopt(
            default_value = "100",
            help = "Number of bases to consider left and right of indel breakpoint when calculating read support. This number should not be too large in order to avoid biases caused by other close variants."
        )]
        indel_window: u32,
        #[structopt(
            default_value = "500",
            help = "Maximum number of observations to use for calling. If locus is exceeding this number, downsampling is performed."
        )]
        max_depth: usize,
    },
    #[structopt(
        name = "filter-calls",
        about = "Filter calls by either controlling the false discovery rate (FDR) at given level, or by posterior odds against the given events."
    )]
    FilterCalls {
        #[structopt(subcommand)]
        method: FilterMethod,
    },
}

#[derive(Debug, StructOpt)]
enum FilterMethod {
    #[structopt(name = "control-fdr")]
    ControlFDR {
        #[structopt(short = "e", help = "Events to consider.")]
        events: Vec<String>,
        #[structopt(short = "a", help = "FDR to control for.")]
        fdr: f64,
        #[structopt(long, help = "Variant type to consider.")]
        vartype: Vec<String>,
        #[structopt(help = "Minimum indel length to consider.")]
        min_len: Option<usize>,
        #[structopt(help = "Maximum indel length to consider (exclusive).")]
        max_len: Option<usize>,
    },
    #[structopt(name = "posterior-odds")]
    PosteriorOdds {
        #[structopt(
            raw(
                possible_values = "{use strum::IntoEnumIterator; &KassRaftery::iter().map(|v| v.into()).collect_vec()}"
            ),
            help = "Kass-Raftery score to filter against."
        )]
        odds: KassRaftery,
        #[structopt(short = "e", help = "Events to consider.")]
        events: Vec<String>,
    },
}

pub fn main() -> Result<(), Box<Error>> {
    let opt = Varlociraptor::from_args();

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
            ref observations,
            max_indel_len,
            exclusive_end,
            max_depth,
            ref reference,
            ref candidates,
            purity,
            ref output,
        } => {
            let tumor_alignment_properties = estimate_alignment_properties(tumor)?;
            let normal_alignment_properties = estimate_alignment_properties(normal)?;
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
                .exclusive_end(exclusive_end)
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
        Varlociraptor::FilterCalls { .. } => {}
    }
    Ok(())
}
