// Copyright 2019 Johannes Köster, Jan Forster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::HashMap;
use std::ops::RangeInclusive;
use std::path::Path;

use anyhow::Result;
use itertools::iproduct;
use bio::stats::{LogProb, PHREDProb, Prob};
use derive_builder::Builder;
use rust_htslib::bcf;
use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::Read;

use lp_modeler::dsl::*;
use lp_modeler::solvers::{CbcSolver, SolverTrait};

#[derive(Builder)]
#[builder(pattern = "owned")]
pub(crate) struct Caller {
    #[builder(private)]
    bcf_reader: bcf::IndexedReader,
    #[builder(private)]
    bcf_writer: bcf::Writer,
    #[builder(private)]
    contig_lens: HashMap<u32, u32>,
    alpha: Prob
}

impl CallerBuilder {

    /// Add alpha value for false discovery rate control of loss-of-heterozygosity calls
    pub(crate) fn add_and_check_alpha(&self, alpha: f64) -> Result<Self> {
        if let correct_alpha= Prob::checked(alpha)? {
            self.alpha = Some(correct_alpha);
            Ok(*self)
        } else {
            panic!("Incorrect alpha specified: {}. Must be 0 <= alpha <= 1].", alpha);
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
        for (contig_id, contig_length) in self.contig_lens {

            // Problem Data
            let intervals = ContigLOHProbs::new(&mut self, &contig_id)?.create_all_intervals();
            let n_intervals = contig_length * (contig_length + 1) / 2;
            let mut intervals_overlap_indicator: HashMap<(RangeInclusive<u64>, RangeInclusive<u64>), Bool> = HashMap::new();
            for (interval1, interval2) in iproduct!(intervals.keys(), intervals.keys()) {
                let key = (interval1, interval2);
                let overlap_indicator = interval1.contains(interval2.start()) | interval1.contains(interval2.end());
                intervals_overlap_indicator.insert(key, overlap_indicator);
            }
            // Define problem and objective sense
            let mut problem = LpProblem::new("LOH segmentation", LpObjective::Maximize);

            // Define Variables
            let interval_loh_indicator: HashMap<RangeInclusive<u64>, LpBinary> =
                intervals.iter()
                    .map(| (&range, &val)| {
                        let key = range;
                        let loh_indicator = LpBinary::new(&format!("{:?}", key));
                        (key, loh_indicator)
                    } )
                    .collect();

            // Define Objective Function: maximise length of selected LOH regions
            let selected_lengths_vec: Vec<LpExpression> = {
                interval_loh_indicator.iter()
                    .map( |(&interval, loh_indicator)| {
                        loh_indicator * (interval.end() - interval.start() + 1)
                } )
            }.collect();
            problem += selected_lengths_vec.sum();

            // Constraint: no overlapping intervals
            for ( current_interval, _) in intervals {
                let n_overlapping_selected: Vec<LpExpression> = {
                    interval_loh_indicator.iter()
                        .map(| (&interval, loh_indicator) |
                            loh_indicator * intervals_overlap_indicator.get(&(current_interval, interval) )
                        ).collect()
                };
                let interval_bound = (1 - interval_loh_indicator.get(&current_interval).unwrap()) * n_intervals;
                problem += n_overlapping_selected.sum().le( interval_bound );
            }
        }

            // Constraint: control false discovery rate at alpha
            let selected_probs_vec: Vec<LpExpression> = {
                interval_loh_indicator.iter()
                    .map( |(&interval, loh_indicator)| {
                        loh_indicator * Prob::from( intervals.get(*interval ).unwrap().ln_one_minus_exp() )
                    } )
            }.collect();
            problem += selected_probs_vec.sum().le(interval_loh_indicator.values().collect().sum() * self.alpha );

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

pub(crate) fn info_phred_to_log_prob(
    record: &mut bcf::Record,
    info_field_name: &[u8]
) -> LogProb {
    if let Ok(prob) = record.info(info_field_name).float() {
        if let Some(_prob) = prob {
            if !_prob[0].is_missing() && !_prob[0].is_nan() {
                let log_prob = LogProb::from(PHREDProb(_prob[0] as f64));
                assert!(
                    log_prob.is_valid(),
                    "invalid PHRED probability '{:?}': {}, at pos: {:?}",
                    info_field_name,
                    _prob[0],
                    record.pos()
                );
                log_prob
            } else {
                panic!(
                    "PHRED probability '{:?}' at pos '{:?}' is missing or NaN",
                    info_field_name,
                    record.pos()
                )
            }
        } else {
            panic!(
                "Expected PHRED probability value in field '{:?}' at pos '{:?}', got None.",
                info_field_name,
                record.pos()
            )
        }
    } else {
        panic!(
            "error unpacking PHRED probability INFO field '{:?}' at pos '{:?}",
            info_field_name,
            record.pos()
        )
    }
}

pub(crate) struct Interval {
    range: RangeInclusive<u64>,
    posterior_prob_loh: LogProb,
}


#[derive(Clone, Debug, Hash, Eq, PartialEq, PartialOrd, Ord)]
pub(crate) struct ContigLOHProbs {
    contig_id: u32,
    length: u64,
    cum_log_prob_loh: Vec<LogProb>,
}

impl ContigLOHProbs {
    pub(crate) fn new(
        caller: &mut Caller,
        contig_id: &u32
    ) -> Result<ContigLOHProbs> {
        let mut record = caller.bcf_reader.empty_record();
        if let Some(contig_length) = caller.contig_lens.get(contig_id) {
            let mut cum_log_prob_loh = Vec::with_capacity(*contig_length as usize + 1);
            caller.bcf_reader.fetch(*contig_id, 0, (contig_length - 1).into());
            cum_log_prob_loh.push(LogProb::ln_one());
            while caller.bcf_reader.read(&mut record)? {
                cum_log_prob_loh.push(cum_log_prob_loh.last().unwrap() + log_prob_loh_given_germ_het(&mut record));
            }
            Ok(
                ContigLOHProbs {
                    contig_id: *contig_id,
                    length: *contig_length as u64,
                    cum_log_prob_loh: cum_log_prob_loh,
                }
            )
        } else {
            panic!(
                "This should not have happened: cannot find contig_id '{}'.",
                contig_id
            )
        }
    }

    pub(crate) fn create_all_intervals(
        &self
    ) -> HashMap<RangeInclusive<u64>, LogProb> {
        let mut intervals= HashMap::new();
        for start in 1..=self.length {
            for end in start..=self.length {
                intervals.insert(
                    start..=end,
                    self.cum_log_prob_loh[ end as usize ] - self.cum_log_prob_loh[ (start - 1) as usize ]
                );
            }
        }
        intervals
    }
}

fn log_prob_loh_given_germ_het(
    record: &mut bcf::Record,
) -> LogProb {
    let log_prob_loh = info_phred_to_log_prob(record, b"PROB_LOH");
    let log_prob_no_loh = info_phred_to_log_prob(record, b"PROB_NO_LOH");
    let log_prob_germline_het = log_prob_loh.ln_add_exp(log_prob_no_loh);
    let log_prob_loh_given_germ_het = log_prob_germline_het.ln_one_minus_exp().ln_add_exp(log_prob_loh + log_prob_germline_het);
    log_prob_loh_given_germ_het
}

// //#[cfg(test)]
// mod tests {
//     use super::*;
//
//     #[test]
//     fn test_allele_freq_pdf() {
//         assert_eq!(
//             allele_freq_pdf(AlleleFreq(0.64), AlleleFreq(1.0), 10),
//             LogProb::ln_zero()
//         );
//         assert_eq!(
//             allele_freq_pdf(AlleleFreq(0.1), AlleleFreq(0.0), 10),
//             LogProb::ln_zero()
//         );
//     }
//
// }
