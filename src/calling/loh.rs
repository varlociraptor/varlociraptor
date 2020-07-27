// Copyright 2019 Johannes Köster, Jan Forster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::{BTreeMap, HashMap};
use std::ops::RangeInclusive;
use std::path::{Path, PathBuf};

use anyhow::Result;
use bio::io::bed;
use bio::stats::{LogProb, PHREDProb, Prob};
use derive_builder::Builder;
use itertools::iproduct;
use ordered_float::NotNan;
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
        let contig_ids: Vec<u32> = self.contig_lens.keys().cloned().collect();
        let mut bed_record = bed::Record::new();
        for contig_id in contig_ids {
            let contig_name =
                String::from_utf8_lossy(self.bcf_reader.header().rid2name(contig_id)?);
            bed_record.set_chrom(&contig_name);
            // Problem Data
            let contig: ContigLogProbsLOH;
            if let Some(contig_length) = self.contig_lens.get(&contig_id) {
                contig = ContigLogProbsLOH::new(&mut self.bcf_reader, &contig_id, contig_length)?;
            } else {
                panic!(
                    "This should not have happened: cannot find contig_id '{}'.",
                    contig_id
                )
            };
            if contig.positions.is_empty() { continue; }
            let intervals = contig.create_all_intervals();
            let mut intervals_overlap_or_adjacent: HashMap<
                (&RangeInclusive<usize>, &RangeInclusive<usize>),
                bool,
            > = HashMap::new();
            for (interval1, interval2) in iproduct!(intervals.keys(), intervals.keys()) {
                let key = (interval1, interval2);
                let start2_minus_one = if interval2.start() >= &1 {
                    interval2.start() - 1
                } else {
                    *interval2.start()
                };
                let overlap_indicator = interval1.contains(interval2.start())
                    | interval1.contains(&start2_minus_one)
                    | interval1.contains(interval2.end())
                    | interval1.contains(&(interval2.end() + 1));
                intervals_overlap_or_adjacent.insert(key, overlap_indicator);
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
            let contig_log_prob: Vec<LpExpression> = {
                interval_loh_indicator
                    .iter()
                    .map(|(&interval, loh_indicator)| {
                        let log_probs = intervals.get(&interval).unwrap();
                        let loh = NotNan::from(log_probs.loh);
                        let no_loh = NotNan::from(log_probs.no_loh);
                        let abs_log_probs_diff: NotNan<f64>;
                        let log_diffs: f32;
                        if loh >= no_loh {
                            abs_log_probs_diff = loh - no_loh;
                            log_diffs = abs_log_probs_diff.into_inner() as f32;
                        } else {
                            abs_log_probs_diff = no_loh - loh;
                            log_diffs = -abs_log_probs_diff.into_inner() as f32;
                        };

                        debug!("{:?} log_diffs: {}", interval, log_diffs);
                        log_diffs * loh_indicator
                    })
                    .collect()
            };
            problem += contig_log_prob.sum();

            // Constraint: no overlapping intervals for selected intervals
            for current_interval in intervals.keys() {
                let n_overlapping_selected: Vec<LpExpression> = {
                    let mut interval_loh_indicator_without_ci = interval_loh_indicator.clone();
                    interval_loh_indicator_without_ci.remove(current_interval);
                    interval_loh_indicator_without_ci
                        .iter()
                        .map(|(&interval, loh_indicator)| {
                            let interval_overlaps: i32 = if *intervals_overlap_or_adjacent
                                .get(&(current_interval, interval))
                                .unwrap()
                            {
                                1
                            } else {
                                0
                            };
                            interval_overlaps * loh_indicator
                        })
                        .collect()
                };
                problem += n_overlapping_selected.sum().le(intervals.len() as f32
                    * (1 - interval_loh_indicator.get(&current_interval).unwrap()));
            }

            // Constraint: control false discovery rate at alpha
            let selected_probs_vec: Vec<LpExpression> = {
                interval_loh_indicator
                    .iter()
                    .map(|(&interval, loh_indicator)| {
                        let interval_prob_no_loh = f64::from(Prob::from(
                            intervals.get(&interval).unwrap().loh.ln_one_minus_exp(),
                        )) as f32;
                        debug!(
                            "interval: {:?}, interval_prob_no_loh: {}",
                            interval, interval_prob_no_loh
                        );
                        (interval_prob_no_loh - f64::from(self.alpha) as f32) * loh_indicator
                    })
                    .collect()
            };
            problem += selected_probs_vec.sum().le(0);

            #[cfg(debug_assertions)]
            problem.write_lp("tests/resources/test_loh/problem.lp")?;

            // Specify solver
            let solver = CbcSolver::new();
            let result = solver.run(&problem);

            // (terminate if error, or assign status & variable values)
            assert!(result.is_ok(), result.unwrap_err());
            let (status, results) = result.unwrap();

            debug!("Status: {:?}", status);
            let mut sorted_records: BTreeMap<u64, bed::Record> = BTreeMap::new();
            for (var_name, var_value) in &results {
                let split: Vec<_> = var_name.split('_').collect();
                let start_index: usize = split[1].parse()?;
                let end_index: usize = split[2].parse()?;
                let int_var_value = *var_value as u32;
                if int_var_value == 1 {
                    debug!("{} is selected", var_name);
                    let score = f64::from(PHREDProb::from(
                        intervals.get(&(start_index..=end_index)).unwrap().loh,
                    ))
                    .to_string();
                    // Write result to bed
                    bed_record.set_start(contig.positions[start_index]); // 0-based as in bcf::Record::pos()
                    bed_record.set_end(contig.positions[end_index] + 1); // 1-based
                    bed_record.set_score(score.as_str());
                    sorted_records.insert(bed_record.start(), bed_record.clone());
                }
            }
            for (_, record) in sorted_records {
                bed_writer.write(&record)?;
            }
        }
        Ok(())
    }
}

#[derive(Debug)]
pub(crate) struct LogProbsLOHnoLOH {
    loh: LogProb,
    no_loh: LogProb,
}

#[derive(Debug)]
pub(crate) struct ContigLogProbsLOH {
    contig_id: u32,
    contig_length: usize,
    cum_loh: Vec<LogProb>,
    cum_no_loh: Vec<LogProb>,
    positions: Vec<u64>,
}

impl ContigLogProbsLOH {
    pub(crate) fn new(
        bcf_reader: &mut bcf::IndexedReader,
        contig_id: &u32,
        contig_length: &usize,
    ) -> Result<ContigLogProbsLOH> {
        let mut record = bcf_reader.empty_record();
        let mut cum_log_prob_loh = Vec::with_capacity(*contig_length);
        let mut cum_log_prob_no_loh = Vec::with_capacity(*contig_length);
        let mut positions = Vec::with_capacity(*contig_length);
        bcf_reader.fetch(*contig_id, 0, (contig_length - 1) as u64)?;
        // put in 1st LOH probability
        if bcf_reader.read(&mut record)? {
            let log_probs = log_probs_loh_no_loh(&mut record);
            cum_log_prob_loh.push(log_probs.loh);
            cum_log_prob_no_loh.push(log_probs.no_loh);
            positions.push(record.pos() as u64);
        }
        // cumulatively add the following LOH probabilities
        while bcf_reader.read(&mut record)? {
            let log_probs = log_probs_loh_no_loh(&mut record);
            cum_log_prob_loh.push(cum_log_prob_loh.last().unwrap() + log_probs.loh);
            cum_log_prob_no_loh.push(cum_log_prob_no_loh.last().unwrap() + log_probs.no_loh);
            positions.push(record.pos() as u64);
        }
        Ok(ContigLogProbsLOH {
            contig_id: *contig_id,
            contig_length: *contig_length,
            cum_loh: cum_log_prob_loh,
            cum_no_loh: cum_log_prob_no_loh,
            positions,
        })
    }

    pub(crate) fn create_all_intervals(&self) -> HashMap<RangeInclusive<usize>, LogProbsLOHnoLOH> {
        let mut intervals = HashMap::new();
        for start in 0..self.cum_loh.len() {
            for end in start..self.cum_loh.len() {
                let (start_log_prob_loh, start_log_prob_no_loh) = if start > 0 {
                    (
                        self.cum_loh[(start - 1) as usize],
                        self.cum_no_loh[(start - 1) as usize],
                    )
                } else {
                    (LogProb::ln_one(), LogProb::ln_one())
                };
                intervals.insert(
                    start..=end,
                    LogProbsLOHnoLOH {
                        loh: self.cum_loh[end as usize] - start_log_prob_loh,
                        no_loh: self.cum_no_loh[end as usize] - start_log_prob_no_loh,
                    },
                );
            }
        }
        intervals
    }
}

fn log_probs_loh_no_loh(record: &mut bcf::Record) -> LogProbsLOHnoLOH {
    let log_prob_loh = info_phred_to_log_prob(record, &String::from("PROB_LOH"));
    let log_prob_no_loh = info_phred_to_log_prob(record, &String::from("PROB_NO_LOH"));
    let log_prob_germline_het = log_prob_loh
        .ln_add_exp(log_prob_no_loh)
        .cap_numerical_overshoot(0.000001);
    let loh =
        (log_prob_germline_het + log_prob_loh).ln_add_exp(log_prob_germline_het.ln_one_minus_exp());
    let no_loh = (log_prob_germline_het + log_prob_no_loh)
        .ln_add_exp(log_prob_germline_het.ln_one_minus_exp());
    LogProbsLOHnoLOH { loh, no_loh }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_log_probs_loh_no_loh() {
        // set up test input
        let test_file = "tests/resources/test_loh/loh.bcf";
        let mut loh_calls = bcf::Reader::from_path(test_file).unwrap();
        let mut record = loh_calls.empty_record();

        let log_probs_loh = [
            LogProb(-0.0000000001113695252626908),
            LogProb(-15.735129611157593),
            LogProb(-0.000026217939537273105),
            LogProb(-0.000026656113726260797),
            LogProb(-0.000026218615627557916),
            LogProb(-7.757570481522067),
            LogProb(-0.00002621896814339697),
        ];

        let mut index = 0;
        while loh_calls.read(&mut record).unwrap() {
            let log_probs = log_probs_loh_no_loh(&mut record);
            println!(
                "LogProb LOH: {:?}, Prob LOH`: {:?}",
                log_probs.loh,
                Prob::from(log_probs.loh)
            );
            assert_eq!(log_probs.loh, log_probs_loh[index]);
            index += 1;
        }
    }

    #[test]
    fn test_loh_caller() {
        let test_input = PathBuf::from("tests/resources/test_loh/loh_no_loh.bcf");
        let test_output = PathBuf::from("tests/resources/test_loh/loh_no_loh.out.bed");
        let expected_bed: Vec<u8> = Vec::from("chr8\t249134\t249135\t\t0.0000006369450033776132\n");
        let alpha = 0.2;
        let mut caller = CallerBuilder::default()
            .bcf(&test_input)
            .unwrap()
            .bed_path(&test_output)
            .add_and_check_alpha(alpha)
            .unwrap()
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
        let expected_bed: Vec<u8> = Vec::from("chr8\t7999999\t8002000\t\t0.00022836455119424118\n");
        let alpha = 0.01;
        let mut caller = CallerBuilder::default()
            .bcf(&test_input)
            .unwrap()
            .bed_path(&test_output)
            .add_and_check_alpha(alpha)
            .unwrap()
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
        let expected_bed: Vec<u8> = Vec::from("chr8\t7999999\t8002000\t\t0.00022836753668196763\n");
        let alpha = 0.01;
        let mut caller = CallerBuilder::default()
            .bcf(&test_input)
            .unwrap()
            .bed_path(&test_output)
            .add_and_check_alpha(alpha)
            .unwrap()
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
            .build()
            .unwrap();

        caller.call().unwrap();
        let produced_bed = fs::read(test_output).expect("Cannot open test output file.");
        assert_eq!(expected_bed, produced_bed);
    }

    #[test]
    fn test_medium_no_loh_no_het_between_loh() {
        let test_input =
            PathBuf::from("tests/resources/test_loh/medium_no_loh_no_het_between_loh.bcf");
        let test_output =
            PathBuf::from("tests/resources/test_loh/medium_no_loh_no_het_between_loh.out.bed");
        let expected_bed : Vec<u8> = Vec::from(
            "chr8\t7999999\t8000000\t\t0.00011386306467910803\nchr8\t8001999\t8002000\t\t0.0000006369450033776132\n"
        );
        let alpha = 0.05;
        let mut caller = CallerBuilder::default()
            .bcf(&test_input)
            .unwrap()
            .bed_path(&test_output)
            .add_and_check_alpha(alpha)
            .unwrap()
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
        let expected_bed: Vec<u8> = Vec::from("chr8\t8000999\t8001000\t\t0.798912553509115\n");
        let alpha = 0.2;
        let mut caller = CallerBuilder::default()
            .bcf(&test_input)
            .unwrap()
            .bed_path(&test_output)
            .add_and_check_alpha(alpha)
            .unwrap()
            .build()
            .unwrap();

        caller.call().unwrap();
        let produced_bed = fs::read(test_output).expect("Cannot open test output file.");
        assert_eq!(expected_bed, produced_bed);
    }

    #[test]
    fn test_contig_loh_probs() {
        let test_input = "tests/resources/test_loh/loh.bcf";
        let mut loh_calls = bcf::IndexedReader::from_path(test_input).unwrap();
        let contig_id = loh_calls.header().name2rid(b"chr8").unwrap();
        let contig_length = 145138636;
        let loh_probs = ContigLogProbsLOH::new(&mut loh_calls, &contig_id, &contig_length).unwrap();

        assert_eq!(loh_probs.contig_length, 145138636);
        assert_eq!(loh_probs.contig_id, 0);
        for log_prob in loh_probs.cum_loh.clone() {
            println!(
                "contig '{:?}' (length '{:?}'), cum_log_prob_loh: '{:?}'",
                loh_probs.contig_id, loh_probs.contig_length, log_prob
            );
        }
        assert_eq!(
            loh_probs.cum_loh.last().unwrap(),
            &LogProb(-23.492805404428065)
        );

        let intervals = loh_probs.create_all_intervals();
        assert_eq!(
            intervals.get(&(0..=4)).unwrap().loh,
            LogProb(-15.735208703937852)
        );
        assert_eq!(
            intervals.get(&(0..=0)).unwrap().loh,
            LogProb(-0.0000000001113695252626908)
        );
        assert_relative_eq!(
            f64::from(*intervals.get(&(4..=4)).unwrap().loh),
            (-0.000026218615627557916 as f64),
            epsilon = 0.000000000000001
        );
        assert_relative_eq!(
            f64::from(*intervals.get(&(1..=4)).unwrap().loh),
            (-15.73520870382648 as f64),
            epsilon = 0.0000000000001
        );
    }
}
