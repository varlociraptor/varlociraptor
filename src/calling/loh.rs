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
use lp_modeler::solvers::{CbcSolver, SolverTrait};

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
            let (intervals, log_probs) = contig.create_all_intervals(
                self.alpha,
                &self.control_local_fdr,
                &self.filter_bayes_factor_minimum_barely,
            )?;

            // Define problem and objective sense
            let mut problem = LpProblem::new("LOH segmentation", LpObjective::Maximize);

            // Define Variables
            let interval_loh_indicator: HashMap<&RangeInclusive<usize>, LpBinary> = intervals
                .iter()
                .enumerate()
                .map(|(index, range)| {
                    let key = range;
                    let loh_indicator = LpBinary::new(&format!("x{}", index));
                    (key, loh_indicator)
                })
                .collect();

            // Define Objective Function: maximise length of selected LOH regions
            let lengths: Vec<LpExpression> = {
                interval_loh_indicator
                    .iter()
                    .map(|(&interval, loh_indicator)| {
                        let interval_length = (interval.end() - interval.start() + 1) as f32;
                        interval_length * loh_indicator
                    })
                    .collect()
            };
            problem += lengths.sum();

            // Constraint: no overlapping intervals for selected intervals
            for current_interval in intervals.iter() {
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
            let selected_probs_vec: Vec<LpExpression> = {
                interval_loh_indicator
                    .iter()
                    .map(|(&interval, loh_indicator)| {
                        let index: usize = loh_indicator
                            .name
                            .strip_prefix('x')
                            .unwrap()
                            .parse()
                            .unwrap();
                        let interval_prob_no_loh =
                            f64::from(Prob::from(log_probs[index].ln_one_minus_exp())) as f32;
                        (interval_prob_no_loh - f64::from(self.alpha) as f32) * loh_indicator
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
            let solver = CbcSolver::new();

            // (terminate if error, or assign status & variable values)
            match solver.run(&problem) {
                Ok(solution) => {
                    debug!("Status: {:?}", solution.status);
                    let mut sorted_records: BTreeMap<u64, bed::Record> = BTreeMap::new();
                    for (var_name, var_value) in solution.results.iter() {
                        let int_var_value = *var_value as u32;
                        if int_var_value == 1 {
                            debug!("{} is selected", var_name);
                            if let Some(index_from_name) = var_name.strip_prefix('x') {
                                let interval_index: usize = index_from_name.parse()?;
                                let start_index = *intervals[interval_index].start();
                                let end_index = *intervals[interval_index].end();
                                let score = f64::from(
                                    PHREDProb::from(
                                        log_probs[interval_index],
                                    )
                                )
                                    .to_string();
                                // Write result to bed
                                bed_record.set_start(contig.positions[start_index]); // 0-based as in bcf::Record::pos()
                                bed_record.set_end(contig.positions[end_index] + 1); // 1-based
                                bed_record.set_score(score.as_str());
                                sorted_records.insert(bed_record.start(), bed_record.clone());
                            } else {
                                panic!("The solver changed the variable names into something unexpected. They should start with an 'x'.");
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
}

fn site_posterior_loh_or_hom(
    record: &mut bcf::Record,
    loh_field_name: &String,
    no_loh_field_name: &String,
    hom_field_name: &String,
    absent_field_name: &String,
) -> Option<LogProb> {
    let site_likelihood_loh = info_phred_to_log_prob(record, loh_field_name);
    let site_likelihood_no_loh = info_phred_to_log_prob(record, no_loh_field_name);
    let site_likelihood_hom = info_phred_to_log_prob(record, hom_field_name).ln_add_exp(info_phred_to_log_prob(record, absent_field_name));
    // Kass-Raftery evidence of at least barely for a heterozygous site
    // TODO: remove, once we have copy number estimation based on DP field
    if site_likelihood_hom > site_likelihood_loh.ln_add_exp(site_likelihood_no_loh) {
        None
    } else {
        let site_likelihood_loh_or_hom = site_likelihood_loh.ln_add_exp(site_likelihood_hom);
        Some(
            site_likelihood_loh_or_hom
                - (site_likelihood_loh_or_hom.ln_add_exp(site_likelihood_no_loh)),
        )
    }
}

impl ContigLogPosteriorsLOH {
    pub(crate) fn new(
        bcf_reader: &mut bcf::IndexedReader,
        contig_id: &u32,
        contig_length: &usize,
    ) -> Result<ContigLogPosteriorsLOH> {
        let mut record = bcf_reader.empty_record();
        let mut cum_loh_posteriors: Vec<LogProb> = Vec::new();
        let mut positions = Vec::new();
        let loh_field_name = &String::from("PROB_LOH");
        let no_loh_field_name = &String::from("PROB_NO_LOH");
        let hom_field_name = &String::from("PROB_UNINTERESTING");
        let absent_field_name = &String::from("PROB_ABSENT");
        bcf_reader.fetch(*contig_id, 0, (contig_length - 1) as u64)?;
        // put in 1st LOH probability
        if bcf_reader.read(&mut record)? {
            let mut posterior = site_posterior_loh_or_hom(
                &mut record,
                loh_field_name,
                no_loh_field_name,
                hom_field_name,
                absent_field_name,
            );
            while posterior.is_none() {
                bcf_reader.read(&mut record)?;
                posterior = site_posterior_loh_or_hom(
                    &mut record,
                    loh_field_name,
                    no_loh_field_name,
                    hom_field_name,
                    absent_field_name,
                );
            }
            match posterior {
                Some(p) => {
                    cum_loh_posteriors.push(p);
                    positions.push(record.pos() as u64);
                },
                None => eprintln!("Found no records with at least barely heterozygous evidence on contig with ID: {}", contig_id)
            }
        } else {
            // no records found
            eprintln!("Found no records on contig with ID: {}", contig_id);
        }
        // cumulatively add the following LOH probabilities
        while bcf_reader.read(&mut record)? {
            let mut posterior = site_posterior_loh_or_hom(
                &mut record,
                loh_field_name,
                no_loh_field_name,
                hom_field_name,
                absent_field_name,
            );
            while posterior.is_none() {
                bcf_reader.read(&mut record)?;
                posterior = site_posterior_loh_or_hom(
                    &mut record,
                    loh_field_name,
                    no_loh_field_name,
                    hom_field_name,
                    absent_field_name,
                );
            }
            match posterior {
                Some(p) => {
                    cum_loh_posteriors.push(cum_loh_posteriors.last().unwrap() + p);
                    positions.push(record.pos() as u64);
                },
                None => eprintln!("Found only one record with at least barely heterozygous evidence on contig with ID: {}", contig_id)
            }
        }
        Ok(ContigLogPosteriorsLOH {
            contig_id: *contig_id,
            contig_length: *contig_length,
            cum_loh_posteriors,
            positions,
        })
    }

    pub(crate) fn create_all_intervals(
        &self,
        alpha: Prob,
        control_local_fdr: &bool,
        filter_bayes_factor_minimum_barely: &bool,
    ) -> Result<(Vec<RangeInclusive<usize>>, Vec<LogProb>)> {
        let log_one_minus_alpha = LogProb::from(alpha).ln_one_minus_exp();
        let mut intervals = Vec::new();
        let mut log_probs = Vec::new();
        for start in 0..self.cum_loh_posteriors.len() {
            for end in start..self.cum_loh_posteriors.len() {
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
                intervals.push(start..=end);
                log_probs.push(posterior_log_probability);
            }
        }
        Ok((intervals, log_probs))
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
        let expected_bed: Vec<u8> = Vec::from("chr8\t240134\t240135\t\t0.0000006369496848265032\n");
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

        caller.call().unwrap();
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
            Vec::from("chr8\t7999999\t8002000\t\t0.0000006369733111102624\n");
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

        caller.call().unwrap();
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
            Vec::from("chr8\t7999999\t8002000\t\t0.0000013527459286120604\n");
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

        caller.call().unwrap();
        let produced_bed = fs::read(test_output).expect("Cannot open test output file.");
        assert_eq!(expected_bed, produced_bed);
    }

    #[test]
    fn test_slightly_no_loh_no_het_between_loh() {
        let test_input =
            PathBuf::from("tests/resources/test_loh/slightly_no_loh_no_het_between_loh.bcf");
        let test_output =
            PathBuf::from("tests/resources/test_loh/slightly_no_loh_no_het_between_loh.out.bed");
        let expected_bed: Vec<u8> = Vec::from("chr8\t7999999\t8000000\t\t0.000000000023626283759240104\n\
                                                    chr8\t8001999\t8002000\t\t0.0000006369496861517909\n");
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

        caller.call().unwrap();
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

        caller.call().unwrap();
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

        caller.call().unwrap();
        let produced_bed = fs::read(test_output).expect("Cannot open test output file.");
        assert_eq!(expected_bed, produced_bed);

        let test_output2 =
            PathBuf::from("tests/resources/test_loh/slightly_loh_artifact_between_no_loh.out2.bed");
        let expected_bed2: Vec<u8> = Vec::from("chr8\t8000999\t8001000\t\t2.691023572732213\n");
        let alpha2 = 0.55;
        let mut caller2 = CallerBuilder::default()
            .bcf(&test_input)
            .unwrap()
            .bed_path(&test_output2)
            .add_and_check_alpha(alpha2)
            .unwrap()
            .control_local_fdr(false)
            .filter_bayes_factor_minimum_barely(false)
            .problems_folder(None)
            .build()
            .unwrap();

        caller2.call().unwrap();
        let produced_bed2 = fs::read(test_output2).expect("Cannot open test output file.");
        assert_eq!(expected_bed2, produced_bed2);
    }

    #[test]
    fn test_medium_no_loh_no_het_between_loh() {
        let test_input =
            PathBuf::from("tests/resources/test_loh/medium_no_loh_no_het_between_loh.bcf");
        let test_output =
            PathBuf::from("tests/resources/test_loh/medium_no_loh_no_het_between_loh.out.bed");
        let expected_bed: Vec<u8> =
            Vec::from("chr8\t7999999\t8002000\t\t0.0000006369733111102624\n");
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

        caller.call().unwrap();
        let produced_bed = fs::read(test_output).expect("Cannot open test output file.");
        assert_eq!(expected_bed, produced_bed);
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

        caller.call().unwrap();
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
        // The two records in loh_no_loh.bcf are mirrored regarding PROB_LOH and PROB_NO_LOH, thus
        // the Prob(LOH = 1) ( .last() in the vector) and Prob(LOH  = 0) ( .first() in the vector)
        // should be equal cumulative sums at the second ( .last() ) position.
        assert_eq!(
            loh_log_posteriors.cum_loh_posteriors[0],
            LogProb(-0.0000001466630849268762)
        );
        assert_eq!(
            loh_log_posteriors.cum_loh_posteriors[1],
            LogProb(-15.735129302014165)
        );

        let (intervals, log_probs) = loh_log_posteriors
            .create_all_intervals(alpha, &false, &false)
            .unwrap();
        println!("intervals: {:?}", intervals);
        assert_eq!(intervals[2], 1..=1);
        assert_eq!(log_probs[2], LogProb(-15.73512915535108));
        assert_eq!(intervals[1], 0..=1);
        assert_eq!(log_probs[1], LogProb(-15.735129302014165));
        assert_eq!(intervals[0], 0..=0);
        assert_eq!(log_probs[0], LogProb(-0.0000001466630849268762));

        let (intervals_filter_bayes_factor, _) = loh_log_posteriors
            .create_all_intervals(alpha, &false, &true)
            .unwrap();
        println!(
            "intervals after Bayes factor filtering: {:?}",
            intervals_filter_bayes_factor
        );
        assert_eq!(intervals_filter_bayes_factor.len(), 1);

        let (intervals_control_local_fdr, _) = loh_log_posteriors
            .create_all_intervals(alpha, &true, &false)
            .unwrap();
        println!(
            "intervals after local FDR filtering: {:?}",
            intervals_control_local_fdr
        );
        assert_eq!(intervals_control_local_fdr.len(), 1);
    }
}
