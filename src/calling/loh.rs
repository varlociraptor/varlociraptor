// Copyright 2019 Johannes Köster, Jan Forster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::HashMap;
use std::ops::RangeInclusive;
use std::path::Path;

use anyhow::Result;
use itertools::iproduct;
use bio::stats::{LogProb, Prob};
use derive_builder::Builder;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;

use lp_modeler::dsl::*;
use lp_modeler::solvers::{CbcSolver, SolverTrait};

use crate::utils::info_phred_to_log_prob;
use std::cmp::max;

#[derive(Builder)]
#[builder(pattern = "owned")]
pub(crate) struct Caller {
    #[builder(private)]
    bcf_reader: bcf::IndexedReader,
    #[builder(private)]
    bcf_writer: bcf::Writer,
    #[builder(private)]
    contig_lens: HashMap<u32, u64>,
    alpha: Prob
}

impl CallerBuilder {

    /// Add alpha value for false discovery rate control of loss-of-heterozygosity calls
    pub(crate) fn add_and_check_alpha(mut self, alpha: f64) -> Result<Self> {
        match Prob::checked(alpha) {
            Ok(correct_alpha) => {
                self.alpha = Some(correct_alpha);
                Ok(self)
            }
            Err(_err) => {
                panic!("Incorrect alpha specified: {}. Must be 0 <= alpha <= 1].", alpha);
            }
        }
    }

    pub(crate) fn bcfs<P: AsRef<Path>>(mut self, in_path: P, out_path: Option<P>) -> Result<Self> {
        self = self.bcf_reader( bcf::IndexedReader::from_path( in_path )? );

        let bcf_reader = self.bcf_reader.as_ref().unwrap();

        let mut header = bcf::Header::new();
        for sample in bcf_reader.header().samples() {
            header.push_sample(sample);
        }

        header.push_record(
            "##INFO=<ID=LOHEND,Number=1,Type=Integer,Description=\"Last variant position supporting loss-of-heterozygosity region.\">"
                .as_bytes(),
        );

        let mut contig_lens = HashMap::new();
        // register sequences
        for rec in bcf_reader.header().header_records() {
            if let bcf::header::HeaderRecord::Contig { values, .. } = rec {
                let name = values.get("ID").unwrap();
                let len = values.get("length").unwrap();
                contig_lens.insert(bcf_reader.header().name2rid( name.as_bytes() )?, len.parse()?);
                header.push_record(format!("##contig=<ID={},length={}>", name, len).as_bytes());
            }
        }

        self = self.contig_lens(contig_lens);

        Ok(self.bcf_writer(if let Some(path) = out_path {
            bcf::Writer::from_path(path, &header, false, bcf::Format::BCF)?
        } else {
            bcf::Writer::from_stdout(&header, false, bcf::Format::BCF)?
        }))
    }
}


impl Caller {
    pub(crate) fn call(&mut self) -> Result<()> {
        let contig_ids: Vec<u32> = self.contig_lens.keys().cloned().collect();
        for contig_id in contig_ids {

            // Problem Data
            let contig: ContigLOHProbs;
            if let Some(contig_length) = self.contig_lens.get(&contig_id) {
                contig = ContigLOHProbs::new(&mut self.bcf_reader, &contig_id, contig_length)?;
            } else {
                panic!(
                    "This should not have happened: cannot find contig_id '{}'.",
                    contig_id
                )
            };
            let intervals = contig.create_all_intervals();
            let n_intervals = (contig.contig_length * (contig.contig_length + 1) / 2) as i32;
            let mut intervals_overlap_indicator: HashMap<(&RangeInclusive<usize>, &RangeInclusive<usize>), bool> = HashMap::new();
            for (interval1, interval2) in iproduct!(intervals.keys(), intervals.keys()) {
                let key = (interval1, interval2);
                let overlap_indicator = interval1.contains(interval2.start()) | interval1.contains(interval2.end());
                intervals_overlap_indicator.insert(key, overlap_indicator);
            }
            // Define problem and objective sense
            let mut problem = LpProblem::new("LOH segmentation", LpObjective::Maximize);

            // Define Variables
            let interval_loh_indicator: HashMap< &RangeInclusive<usize>, LpBinary > =
                intervals.iter()
                    .map(|(range, _)| {
                        let key = range;
                        let loh_indicator = LpBinary::new(&format!("{:?}", key));
                        (key, loh_indicator)
                    })
                    .collect();

            // Define Objective Function: maximise length of selected LOH regions
            let selected_lengths_vec: Vec<LpExpression> = {
                interval_loh_indicator.iter()
                    .map(|(&interval, loh_indicator)| {
                        let interval_length = (interval.end() - interval.start() + 1) as i32;
                        interval_length * loh_indicator
                    }).collect()
            };
            problem += selected_lengths_vec.sum();

            // Constraint: no overlapping intervals
            for (current_interval, _) in &intervals {
                let n_overlapping_selected: Vec<LpExpression> = {
                    interval_loh_indicator.iter()
                        .map(|(&interval, loh_indicator)| {
                            let interval_overlaps: i32 = if *intervals_overlap_indicator.get(&(current_interval, interval)).unwrap() {
                                1
                            } else {
                                0
                            };
                            interval_overlaps * loh_indicator
                        }
                        ).collect()
                };
                problem += n_overlapping_selected.sum()
                    .le(n_intervals * (1 - interval_loh_indicator.get(&current_interval).unwrap() ) );
            }

            // Constraint: control false discovery rate at alpha
            let selected_probs_vec: Vec<LpExpression> = {
                interval_loh_indicator.iter()
                    .map(|(&interval, loh_indicator)| {
                        let interval_prob = f64::from( Prob::from(intervals.get(&interval).unwrap().ln_one_minus_exp() ) ) as f32;
                        interval_prob * loh_indicator
                    }).collect()
            };
            let selected: Vec<&LpBinary> = interval_loh_indicator.values().collect();
            problem += selected_probs_vec.sum().le(f64::from( self.alpha ) as f32 * selected.sum());

            // Specify solver
            let solver = CbcSolver::new();

            // Run optimisation and process output hashmap
            // Write result to bcf
            match solver.run(&problem) {
                Ok((status, var_values)) => {
                    println!("Status {:?}", status);
                    for (name, value) in var_values.iter() {
                        println!("value of {} = {}", name, value);
                    }
                },
                Err(msg) => println!("{}", msg),
            }
        }
        Ok(())
    }
}


#[derive(Debug) ]
pub(crate) struct ContigLOHProbs {
    contig_id: u32,
    contig_length: u64,
    cum_log_prob_loh: Vec<LogProb>,
}

impl ContigLOHProbs {
    pub(crate) fn new(
        bcf_reader: &mut bcf::IndexedReader,
        contig_id: &u32,
        contig_length: &u64
    ) -> Result<ContigLOHProbs> {
        let mut record = bcf_reader.empty_record();
        let mut cum_log_prob_loh = Vec::with_capacity(*contig_length as usize);
        bcf_reader.fetch(*contig_id, 0, (contig_length - 1).into())?;
        // put in 1st LOH probability
        if bcf_reader.read(&mut record)? {
            cum_log_prob_loh.push(log_prob_loh_given_germ_het(&mut record));
        }
        // cumulatively add the following LOH probabilities
        while bcf_reader.read(&mut record)? {
            cum_log_prob_loh.push(cum_log_prob_loh.last().unwrap() + log_prob_loh_given_germ_het(&mut record));
        }
        Ok(
            ContigLOHProbs {
                contig_id: *contig_id,
                contig_length: *contig_length as u64,
                cum_log_prob_loh: cum_log_prob_loh,
            }
        )
    }

    pub(crate) fn create_all_intervals(
        &self
    ) -> HashMap<RangeInclusive<usize>, LogProb> {
        let mut intervals= HashMap::new();
        for start in 0..self.cum_log_prob_loh.len() {
            for end in start..self.cum_log_prob_loh.len() {
                let start_log_prob = if start > 0 {
                    self.cum_log_prob_loh[ (start - 1) as usize ]
                } else {
                    LogProb::ln_one()
                };
                intervals.insert(
                    start..=end,
                    self.cum_log_prob_loh[ end as usize ] - start_log_prob
                );
            }
        }
        intervals
    }
}

fn log_prob_loh_given_germ_het(
    record: &mut bcf::Record,
) -> LogProb {
    let log_prob_loh = info_phred_to_log_prob(record, &String::from("PROB_LOH") );
    let log_prob_no_loh = info_phred_to_log_prob(record, &String::from("PROB_NO_LOH") );
    let log_prob_germline_het = log_prob_loh
        .ln_add_exp(log_prob_no_loh)
        .cap_numerical_overshoot(0.000001);
    let log_prob_loh_given_germ_het = log_prob_germline_het.ln_one_minus_exp().ln_add_exp(log_prob_loh + log_prob_germline_het);
    log_prob_loh_given_germ_het
}

//#[cfg(test)]
mod tests {
    use super::*;
    use std::path::{Path, PathBuf};

    #[test]
    fn test_log_prob_loh_given_germ_het() {
        // set up test input
        let test_file = "tests/resources/loh.bcf";
        let mut loh_calls = bcf::Reader::from_path(test_file).unwrap();
        let mut record = loh_calls.empty_record();

        let log_probs = [
            LogProb(-15.735129611157593),
            LogProb(-0.0000000001113695252626908),
            LogProb(-0.00002621896814339697),
            LogProb(-0.0004275831483212109),
            LogProb(-0.000026656113726260797),
        ];

        let mut index = 0;
        while loh_calls.read(&mut record).unwrap() {
            assert_eq!(
                log_prob_loh_given_germ_het(&mut record),
                log_probs[index]
            );
            index += 1;
        }
    }

//    #[test]
//    fn test_loh_caller() {
//        let test_input = PathBuf("tests/resources/loh.vcf");
//        let test_output = Some(PathBuf::from("tests/resources/loh.out.vcf"));
//        let alpha = 0.2;
//        let mut caller = CallerBuilder::default()
//            .bcfs(&test_input, test_output.as_ref())?
//            .add_and_check_alpha(alpha)?
//            .build()
//            .unwrap();
//
//        caller.call()?;
//    }

    #[test]
    fn test_contig_loh_probs() {
        let test_input = "tests/resources/loh.bcf";
        let mut loh_calls = bcf::IndexedReader::from_path(test_input).unwrap();
        let contig_id = loh_calls.header().name2rid( b"chr8" ).unwrap();
        let contig_length = 145138636;
        let loh_probs = ContigLOHProbs::new(&mut loh_calls, &contig_id, &contig_length).unwrap();

        assert_eq!(
            loh_probs.contig_length,
            145138636
        );
        assert_eq!(
            loh_probs.contig_id,
            0
        );
        for log_prob in loh_probs.cum_log_prob_loh.clone() {
            println!("contig '{:?}' (length '{:?}'), cum_log_prob_loh: '{:?}'", loh_probs.contig_id, loh_probs.contig_length, log_prob);
        }
        assert_eq!(
            loh_probs.cum_log_prob_loh.last().unwrap(),
            &LogProb(-15.735610069499153)
        );

        let intervals = loh_probs.create_all_intervals();
        assert_eq!(
            intervals.get( &(0..=4) ).unwrap(),
            &LogProb(-15.735610069499153)
        );
        assert_eq!(
            intervals.get( &(0..=0) ).unwrap(),
            &LogProb(-15.735129611157593)
        );
        assert_relative_eq!(
            f64::from( *intervals.get( &(4..=4) ).unwrap() ),
            &(-0.000026656113726260797 as f64),
            epsilon = 0.000000000000001
        );
        assert_relative_eq!(
            f64::from( *intervals.get( &(1..=4) ).unwrap() ),
            &(-0.00048045834156 as f64),
            epsilon = 0.0000000000001
        );
    }

}
