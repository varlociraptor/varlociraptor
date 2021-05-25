// Copyright 2019 Johannes Köster, Jan Forster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::{BTreeMap, HashMap};
use std::fs::create_dir_all;
use std::ops::RangeInclusive;
use std::path::{Path, PathBuf};

use anyhow::Result;
use bio::io::bed;
use bio::stats::bayesian::bayes_factors::{evidence, BayesFactor};
use bio::stats::{LogProb, PHREDProb, Prob};
use derive_builder::Builder;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;

use lp_modeler::dsl::*;
use lp_modeler::format::lp_format::LpFileFormat;
use lp_modeler::solvers::{GurobiSolver, SolverTrait};

use crate::utils::info_phred_to_log_prob;

#[derive(Builder)]
#[builder(pattern = "owned")]
pub(crate) struct Caller<'a> {
    #[builder(private)]
    bcf_reader: bcf::IndexedReader,
    bed_path: &'a PathBuf,
    #[builder(private)]
    contig_lens: HashMap<u32, usize>,
    #[builder(private)]
    alpha: Prob,
    control_local_fdr: bool,
    filter_bayes_factor_minimum_barely: bool,
    problems_folder: Option<PathBuf>,
}

impl CallerBuilder<'_> {
    /// Add alpha value for false discovery rate control of loss-of-heterozygosity calls
    pub(crate) fn add_and_check_alpha(mut self, alpha: f64) -> Result<Self> {
        self.alpha = Some(Prob::checked(alpha)?);
        Ok(self)
    }

    pub(crate) fn bcf<P: AsRef<Path>>(mut self, in_path: P) -> Result<Self> {
        self = self.bcf_reader(bcf::IndexedReader::from_path(in_path)?);

        let bcf_reader = self.bcf_reader.as_ref().unwrap();

        let mut contig_lens = HashMap::new();
        // collect reference sequences / contigs
        for rec in bcf_reader.header().header_records() {
            if let bcf::header::HeaderRecord::Contig { values, .. } = rec {
                let name = values.get("ID").unwrap();
                let len = values.get("length").unwrap();
                contig_lens.insert(bcf_reader.header().name2rid(name.as_bytes())?, len.parse()?);
            }
        }

        self = self.contig_lens(contig_lens);

        Ok(self)
    }
}

impl Caller<'_> {
    pub(crate) fn call(&mut self) -> Result<()> {
        let mut bed_writer = bed::Writer::to_file(self.bed_path)?;
        let write_problems = match &self.problems_folder {
            Some(output_path) => match create_dir_all(output_path) {
                Ok(()) => true,
                Err(err) => {
                    eprintln!("Could not create the directory specified via the --problems_folder. \
                                    Will not output the per contig linear programming problem formulations. \
                                    Error is: '{}'", err);
                    false
                }
            },
            None => false,
        };
        let contig_ids: Vec<u32> = self.contig_lens.keys().cloned().collect();
        let mut bed_record = bed::Record::new();
        for contig_id in contig_ids {
            let contig_name =
                String::from_utf8_lossy(self.bcf_reader.header().rid2name(contig_id)?).into_owned();
            bed_record.set_chrom(&contig_name);
            // Problem Data
            let contig: ContigLogPosteriorsLOH;
            if let Some(contig_length) = self.contig_lens.get(&contig_id) {
                contig =
                    ContigLogPosteriorsLOH::new(&mut self.bcf_reader, &contig_id, contig_length)?;
            } else {
                panic!(
                    "This should not have happened: cannot find contig_id '{}'.",
                    contig_id
                )
            };
            if contig.positions.is_empty() {
                continue;
            }
            let intervals = contig.create_all_intervals(
                self.alpha,
                &self.control_local_fdr,
                &self.filter_bayes_factor_minimum_barely,
            )?;
            if intervals.len() == 0 {
                eprintln!("No LOH candidate intervals found. BED file will be empty.");
                return Ok(());
            }

            // Define problem and objective sense
            let mut problem = LpProblem::new("LOH segmentation", LpObjective::Maximize);

            // Define Variables
            let interval_loh_indicator: HashMap<&RangeInclusive<usize>, LpBinary> = intervals
                .iter()
                .map(|(range, _)| {
                    let key = range;
                    let loh_indicator = LpBinary::new(&format!("i_{}_{}", key.start(), key.end()));
                    (key, loh_indicator)
                })
                .collect();

            // Define Objective Function: maximise length of selected LOH regions
            let lengths: Vec<LpExpression> = {
                interval_loh_indicator
                    .iter()
                    .map(|(&interval, loh_indicator)| {
                        let interval_length = (contig.positions[*interval.end()]
                            - contig.positions[*interval.start()]
                            + 1) as f32;
                        interval_length * loh_indicator
                    })
                    .collect()
            };
            problem += lengths.sum();

            // Constraint: no overlapping intervals for selected intervals
            for current_interval in intervals.keys() {
                let mut n_overlapping_selected: Vec<&LpBinary> = Vec::new();
                for (&interval, loh_indicator) in interval_loh_indicator.iter() {
                    // current interval will overlap itself, but should not add to the constraint
                    if current_interval != interval
                        && (current_interval.contains(interval.start())
                        // adjacency is equivalent to an overlap
                        || current_interval.contains(&(interval.start().saturating_sub(1)))
                        || current_interval.contains(interval.end())
                        // adjacency is equivalent to an overlap
                        || current_interval.contains(&(interval.end() + 1)))
                    {
                        n_overlapping_selected.push(loh_indicator);
                    }
                }
                if n_overlapping_selected.len() > 0 {
                    problem += n_overlapping_selected.sum().le(intervals.len() as f32
                        * (1 - interval_loh_indicator.get(&current_interval).unwrap()));
                }
            }

            // Constraint: control false discovery rate at alpha
            let alpha_f64 = f64::from(self.alpha);
            let selected_probs_vec: Vec<LpExpression> = {
                interval_loh_indicator
                    .iter()
                    .map(|(&interval, loh_indicator)| {
                        let interval_prob_no_loh = (f64::from(Prob::from(
                            intervals.get(&interval).unwrap().ln_one_minus_exp(),
                        )) - alpha_f64) as f32;
                        interval_prob_no_loh * loh_indicator
                    })
                    .collect()
            };
            problem += selected_probs_vec.sum().le(0);

            if write_problems {
                problem.write_lp(&*format!(
                    "{}/{}.problem.lp",
                    self.problems_folder
                        .as_ref()
                        .unwrap()
                        .as_os_str()
                        .to_str()
                        .unwrap(),
                    contig_name
                ))?;
            }

            // Specify solver
            let solver = GurobiSolver::new();

            // (terminate if error, or assign status & variable values)
            match solver.run(&problem) {
                Ok(solution) => {
                    println!("Status: {:?}", solution.status);
                    let mut sorted_records: BTreeMap<u64, bed::Record> = BTreeMap::new();
                    for (var_name, var_value) in solution.results.iter() {
                        let split: Vec<_> = var_name.split('_').collect();
                        if split.len() != 3 {
                            panic!(
                                format!("Current chromosome index: {}\n\
                                        The solver seems to have renamed the variables, probably due to some parsing problem.\n\
                                        Please use the `--problems-folder` option to save the problem as a *.problem.lp file and manually run it through the solver to check this.",
                                        contig_id)
                            );
                        } else {
                            let start_index: usize = split[1].parse()?;
                            let end_index: usize = split[2].parse()?;
                            let int_var_value = *var_value as u32;
                            if int_var_value == 1 {
                                debug!("{} is selected", var_name);
                                let score = f64::from(
                                    PHREDProb::from(
                                        *intervals.get(&(start_index..=end_index)).unwrap(),
                                    )
                                )
                                    .to_string();
                                // Write result to bed
                                bed_record.set_start(contig.positions[start_index]); // 0-based as in bcf::Record::pos()
                                bed_record.set_end(contig.positions[end_index] + 1); // 1-based
                                bed_record.set_score(score.as_str());
                                sorted_records.insert(bed_record.start(), bed_record.clone());
                            }
                        }
                    }
                    for (_, record) in sorted_records {
                        bed_writer.write(&record)?;
                    }
                },
                Err(error) => panic!("Could not open the temporary solution file. Probably, the solver did not find a solution. It returned: '{}'", error),
            }
        }
        Ok(())
    }
}

#[derive(Debug)]
pub(crate) struct ContigLogPosteriorsLOH {
    contig_id: u32,
    contig_length: usize,
    cum_loh_posteriors: Vec<LogProb>,
    positions: Vec<u64>,
    valid_start_ends: Vec<bool>,
}

fn site_posterior_not_no_loh(
    record: &mut bcf::Record,
    loh_field_name: &String,
    no_loh_field_name: &String,
    absent_field_name: &String,
    uninteresting_field_name: &String,
    artifact_field_name: &String,
    minimum_background_het: LogProb,
    valid_start_end_background_het: LogProb,
) -> Option<(LogProb, bool)> {
    let loh_log_prob = info_phred_to_log_prob(record, loh_field_name);
    let no_loh_log_prob = info_phred_to_log_prob(record, no_loh_field_name);
    let background_het = loh_log_prob.ln_add_exp(no_loh_log_prob);
    if background_het < minimum_background_het {
        return None;
    }
    // a LogProb(0.0) would set the cumulative sum of log posteriors to
    // LogProb(-inf), rendering them useless
    let checked_posterior = if no_loh_log_prob > LogProb(-0.0000001) {
        let absent_log_prob = info_phred_to_log_prob(record, absent_field_name);
        let uninteresting_log_prob = info_phred_to_log_prob(record, uninteresting_field_name);
        let artifact_log_prob = info_phred_to_log_prob(record, artifact_field_name);
        loh_log_prob
            .ln_add_exp(absent_log_prob)
            .ln_add_exp(uninteresting_log_prob)
            .ln_add_exp(artifact_log_prob)
    } else {
        no_loh_log_prob.ln_one_minus_exp()
    };
    let valid_start_end =
        loh_log_prob > no_loh_log_prob && background_het > valid_start_end_background_het;
    Some((checked_posterior, valid_start_end))
}

impl ContigLogPosteriorsLOH {
    pub(crate) fn new(
        bcf_reader: &mut bcf::IndexedReader,
        contig_id: &u32,
        contig_length: &usize,
    ) -> Result<ContigLogPosteriorsLOH> {
        let minimum_background_het: LogProb = LogProb::from(Prob(0.2));
        let valid_start_end_background_het: LogProb = LogProb::from(Prob(0.5));
        let mut record = bcf_reader.empty_record();
        let mut cum_loh_posteriors: Vec<LogProb> = Vec::new();
        let mut positions = Vec::new();
        let mut valid_start_ends: Vec<bool> = Vec::new();
        let loh_field_name = &String::from("PROB_LOH");
        let no_loh_field_name = &String::from("PROB_NO_LOH");
        let absent_field_name = &String::from("PROB_ABSENT");
        let uninteresting_field_name = &String::from("PROB_UNINTERESTING");
        let artifact_field_name = &String::from("PROB_ARTIFACT");
        bcf_reader.fetch(*contig_id, 0, (contig_length - 1) as u64)?;
        // put in 1st LOH probability
        while let Some(result) = bcf_reader.read(&mut record) {
            match result {
                Ok(()) => {
                    match site_posterior_not_no_loh(
                        &mut record,
                        loh_field_name,
                        no_loh_field_name,
                        absent_field_name,
                        uninteresting_field_name,
                        artifact_field_name,
                        minimum_background_het,
                        valid_start_end_background_het,
                    ) {
                        Some((posterior, valid_start_end)) => {
                            cum_loh_posteriors.push(posterior);
                            positions.push(record.pos() as u64);
                            valid_start_ends.push(valid_start_end);
                            break;
                        }
                        None => continue,
                    }
                }
                Err(err) => eprintln!(
                    "Error while trying to read records on contig with ID: {}\n Error is: {}",
                    contig_id, err
                ),
            }
        }
        if cum_loh_posteriors.len() == 0 {
            eprintln!(
                "Found no records with at least barely heterozygous evidence on contig with ID: {}",
                contig_id
            );
            return Ok(ContigLogPosteriorsLOH {
                contig_id: *contig_id,
                contig_length: *contig_length,
                cum_loh_posteriors,
                positions,
                valid_start_ends,
            });
        }
        // cumulatively add the following LOH probabilities
        while let Some(result) = bcf_reader.read(&mut record) {
            match result {
                Ok(()) => {
                    match site_posterior_not_no_loh(
                        &mut record,
                        loh_field_name,
                        no_loh_field_name,
                        absent_field_name,
                        uninteresting_field_name,
                        artifact_field_name,
                        minimum_background_het,
                        valid_start_end_background_het,
                    ) {
                        Some((posterior, valid_start_end)) => {
                            cum_loh_posteriors.push(cum_loh_posteriors.last().unwrap() + posterior);
                            positions.push(record.pos() as u64);
                            valid_start_ends.push(valid_start_end);
                        }
                        None => continue,
                    }
                }
                Err(err) => eprintln!(
                    "Error while trying to read records on contig with ID: {}\n Error is: {}",
                    contig_id, err
                ),
            }
        }
        Ok(ContigLogPosteriorsLOH {
            contig_id: *contig_id,
            contig_length: *contig_length,
            cum_loh_posteriors,
            positions,
            valid_start_ends,
        })
    }

    pub(crate) fn create_all_intervals(
        &self,
        alpha: Prob,
        control_local_fdr: &bool,
        filter_bayes_factor_minimum_barely: &bool,
    ) -> Result<HashMap<RangeInclusive<usize>, LogProb>> {
        let log_one_minus_alpha = LogProb::from(alpha).ln_one_minus_exp();
        let mut intervals = HashMap::new();
        for start in 0..(self.cum_loh_posteriors.len() - 1) {
            if !self.valid_start_ends[start] {
                continue;
            }
            for end in (start + 1)..self.cum_loh_posteriors.len() {
                if !self.valid_start_ends[end] {
                    continue;
                }
                let posterior_log_probability = if start > 0 {
                    self.cum_loh_posteriors[end] - self.cum_loh_posteriors[start - 1]
                } else {
                    self.cum_loh_posteriors[end]
                };
                // skipping stuff can make the problem formulation much smaller to start with
                if *control_local_fdr && posterior_log_probability < log_one_minus_alpha {
                    continue;
                }
                if *filter_bayes_factor_minimum_barely
                    && BayesFactor::new(
                        posterior_log_probability,
                        posterior_log_probability.ln_one_minus_exp(),
                    )
                    .evidence_kass_raftery()
                        == evidence::KassRaftery::None
                {
                    continue;
                }
                intervals.insert(start..=end, posterior_log_probability);
            }
        }
        Ok(intervals)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_loh_caller() {
        let test_input = PathBuf::from("tests/resources/test_loh/loh_no_loh.bcf");
        let test_output = PathBuf::from("tests/resources/test_loh/loh_no_loh.out.bed");
        let expected_bed: Vec<u8> = Vec::from("chr8\t240134\t243122\t\t0.0000012738993696486607\n");
        let alpha = 0.98;
        let mut caller = CallerBuilder::default()
            .bcf(&test_input)
            .unwrap()
            .bed_path(&test_output)
            .add_and_check_alpha(alpha)
            .unwrap()
            .control_local_fdr(false)
            .filter_bayes_factor_minimum_barely(false)
            .problems_folder(Some(PathBuf::from(
                "tests/resources/test_loh/test_loh_caller",
            )))
            .build()
            .unwrap();

        match caller.call() {
            Ok(_) => println!("Caller returned successfully"),
            Err(err) => panic!("Caller did not return successfully! Err: {}", err),
        }
        let produced_bed = fs::read(test_output).expect("Cannot open test output file.");
        assert_eq!(expected_bed, produced_bed);
    }

    #[test]
    fn test_slightly_loh_no_het_between_loh() {
        let test_input =
            PathBuf::from("tests/resources/test_loh/slightly_loh_no_het_between_loh.bcf");
        let test_output =
            PathBuf::from("tests/resources/test_loh/slightly_loh_no_het_between_loh.out.bed");
        let expected_bed: Vec<u8> =
            Vec::from("chr8\t7999999\t8002000\t\t0.0000006369733104886611\n");
        let alpha = 0.01;
        let mut caller = CallerBuilder::default()
            .bcf(&test_input)
            .unwrap()
            .bed_path(&test_output)
            .add_and_check_alpha(alpha)
            .unwrap()
            .control_local_fdr(false)
            .filter_bayes_factor_minimum_barely(false)
            .problems_folder(None)
            .build()
            .unwrap();

        match caller.call() {
            Ok(_) => println!("Caller returned successfully"),
            Err(err) => panic!("Caller did not return successfully! Err: {}", err),
        }
        let produced_bed = fs::read(test_output).expect("Cannot open test output file.");
        assert_eq!(expected_bed, produced_bed);
    }

    #[test]
    fn test_slightly_loh_artifact_between_loh() {
        let test_input =
            PathBuf::from("tests/resources/test_loh/slightly_loh_artifact_between_loh.bcf");
        let test_output =
            PathBuf::from("tests/resources/test_loh/slightly_loh_artifact_between_loh.out.bed");
        let expected_bed: Vec<u8> =
            Vec::from("chr8\t7999999\t8002000\t\t0.0000006369733104886611\n");
        let alpha = 0.01;
        let mut caller = CallerBuilder::default()
            .bcf(&test_input)
            .unwrap()
            .bed_path(&test_output)
            .add_and_check_alpha(alpha)
            .unwrap()
            .control_local_fdr(false)
            .filter_bayes_factor_minimum_barely(false)
            .problems_folder(None)
            .build()
            .unwrap();

        match caller.call() {
            Ok(_) => println!("Caller returned successfully"),
            Err(err) => panic!("Caller did not return successfully! Err: {}", err),
        }
        let produced_bed = fs::read(test_output).expect("Cannot open test output file.");
        assert_eq!(expected_bed, produced_bed);
    }

    #[test]
    fn test_slightly_no_loh_no_het_between_loh() {
        let test_input =
            PathBuf::from("tests/resources/test_loh/slightly_no_loh_no_het_between_loh.bcf");
        let test_output =
            PathBuf::from("tests/resources/test_loh/slightly_no_loh_no_het_between_loh.out.bed");
        let expected_bed: Vec<u8> =
            Vec::from("chr8\t7999999\t8002000\t\t0.0000006369733104886611\n");
        let alpha = 0.01;
        let mut caller = CallerBuilder::default()
            .bcf(&test_input)
            .unwrap()
            .bed_path(&test_output)
            .add_and_check_alpha(alpha)
            .unwrap()
            .control_local_fdr(false)
            .filter_bayes_factor_minimum_barely(false)
            .problems_folder(None)
            .build()
            .unwrap();

        match caller.call() {
            Ok(_) => println!("Caller returned successfully"),
            Err(err) => panic!("Caller did not return successfully! Err: {}", err),
        }
        let produced_bed = fs::read(test_output).expect("Cannot open test output file.");
        assert_eq!(expected_bed, produced_bed);
    }

    #[test]
    fn test_slightly_loh_no_het_between_no_loh() {
        let test_input =
            PathBuf::from("tests/resources/test_loh/slightly_loh_no_het_between_no_loh.bcf");
        let test_output =
            PathBuf::from("tests/resources/test_loh/slightly_loh_no_het_between_no_loh.out.bed");
        let expected_bed: Vec<u8> = vec![];
        let alpha = 0.01;
        let mut caller = CallerBuilder::default()
            .bcf(&test_input)
            .unwrap()
            .bed_path(&test_output)
            .add_and_check_alpha(alpha)
            .unwrap()
            .control_local_fdr(false)
            .filter_bayes_factor_minimum_barely(false)
            .problems_folder(None)
            .build()
            .unwrap();

        match caller.call() {
            Ok(_) => println!("Caller returned successfully"),
            Err(err) => panic!("Caller did not return successfully! Err: {}", err),
        }
        let produced_bed = fs::read(test_output).expect("Cannot open test output file.");
        assert_eq!(expected_bed, produced_bed);
    }

    #[test]
    fn test_slightly_loh_artifact_between_no_loh() {
        let test_input =
            PathBuf::from("tests/resources/test_loh/slightly_loh_artifact_between_no_loh.bcf");
        let test_output =
            PathBuf::from("tests/resources/test_loh/slightly_loh_artifact_between_no_loh.out.bed");
        let expected_bed: Vec<u8> = vec![];
        let alpha = 0.2;
        let mut caller = CallerBuilder::default()
            .bcf(&test_input)
            .unwrap()
            .bed_path(&test_output)
            .add_and_check_alpha(alpha)
            .unwrap()
            .control_local_fdr(false)
            .filter_bayes_factor_minimum_barely(false)
            .problems_folder(None)
            .build()
            .unwrap();

        match caller.call() {
            Ok(_) => println!("Caller returned successfully"),
            Err(err) => panic!("Caller did not return successfully! Err: {}", err),
        }
        let produced_bed = fs::read(test_output).expect("Cannot open test output file.");
        assert_eq!(expected_bed, produced_bed);
    }

    #[test]
    fn test_medium_no_loh_no_het_between_loh() {
        let test_input =
            PathBuf::from("tests/resources/test_loh/medium_no_loh_no_het_between_loh.bcf");
        let test_output_1 =
            PathBuf::from("tests/resources/test_loh/medium_no_loh_no_het_between_loh.out1.bed");
        let expected_bed_1: Vec<u8> = Vec::from("chr8\t7999999\t8002245\t\t0.7573838081615665\n");
        let alpha_1 = 0.2;
        let mut caller = CallerBuilder::default()
            .bcf(&test_input)
            .unwrap()
            .bed_path(&test_output_1)
            .add_and_check_alpha(alpha_1)
            .unwrap()
            .control_local_fdr(false)
            .filter_bayes_factor_minimum_barely(false)
            .problems_folder(None)
            .build()
            .unwrap();

        match caller.call() {
            Ok(_) => println!("Caller returned successfully"),
            Err(err) => panic!("Caller did not return successfully! Err: {}", err),
        }
        let produced_bed_1 = fs::read(test_output_1).expect("Cannot open test output file.");
        assert_eq!(expected_bed_1, produced_bed_1);
        let test_output_2 =
            PathBuf::from("tests/resources/test_loh/medium_no_loh_no_het_between_loh.out2.bed");
        let expected_bed_2: Vec<u8> = Vec::from(
            "chr8\t7999999\t8000611\t\t0.0000006369733104886611\n\
                       chr8\t8001999\t8002245\t\t0.0000006369733104871483\n",
        );
        let alpha_2 = 0.05;
        let mut caller = CallerBuilder::default()
            .bcf(&test_input)
            .unwrap()
            .bed_path(&test_output_2)
            .add_and_check_alpha(alpha_2)
            .unwrap()
            .control_local_fdr(false)
            .filter_bayes_factor_minimum_barely(false)
            .problems_folder(None)
            .build()
            .unwrap();

        match caller.call() {
            Ok(_) => println!("Caller returned successfully"),
            Err(err) => panic!("Caller did not return successfully! Err: {}", err),
        }
        let produced_bed_2 = fs::read(test_output_2).expect("Cannot open test output file.");
        assert_eq!(expected_bed_2, produced_bed_2);
    }

    #[test]
    fn test_medium_loh_no_het_between_no_loh() {
        let test_input =
            PathBuf::from("tests/resources/test_loh/medium_loh_no_het_between_no_loh.bcf");
        let test_output =
            PathBuf::from("tests/resources/test_loh/medium_loh_no_het_between_no_loh.out.bed");
        let expected_bed: Vec<u8> = vec![];
        let alpha = 0.05;
        let mut caller = CallerBuilder::default()
            .bcf(&test_input)
            .unwrap()
            .bed_path(&test_output)
            .add_and_check_alpha(alpha)
            .unwrap()
            .control_local_fdr(false)
            .filter_bayes_factor_minimum_barely(false)
            .problems_folder(None)
            .build()
            .unwrap();

        match caller.call() {
            Ok(_) => println!("Caller returned successfully"),
            Err(err) => panic!("Caller did not return successfully! Err: {}", err),
        }
        let produced_bed = fs::read(test_output).expect("Cannot open test output file.");
        assert_eq!(expected_bed, produced_bed);
    }

    #[test]
    fn test_contig_loh_probs() {
        let test_input = "tests/resources/test_loh/loh_no_loh.bcf";
        let mut loh_calls = bcf::IndexedReader::from_path(test_input).unwrap();
        let contig_id = loh_calls.header().name2rid(b"chr8").unwrap();
        let contig_length = 145138636;
        let alpha = Prob(0.01);
        let loh_log_posteriors =
            ContigLogPosteriorsLOH::new(&mut loh_calls, &contig_id, &contig_length).unwrap();

        assert_eq!(loh_log_posteriors.contig_length, 145138636);
        assert_eq!(loh_log_posteriors.contig_id, 0);
        assert_eq!(
            loh_log_posteriors.cum_loh_posteriors[0],
            LogProb(-0.0000001466630849263759)
        );
        assert_eq!(
            loh_log_posteriors.cum_loh_posteriors[2],
            LogProb(-15.735135536385188)
        );

        let intervals = loh_log_posteriors
            .create_all_intervals(alpha, &false, &false)
            .unwrap();
        println!("intervals: {:?}", intervals);
        assert_eq!(
            intervals.get(&(0..=1)).unwrap(),
            &LogProb(-0.0000002933261698527518)
        );
    }
}
