// Copyright 2019 Johannes Köster, Jan Forster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::{BTreeMap, HashMap};
use std::error::Error;
use std::iter;
use std::mem;
use std::path::Path;

use bio::stats::{hmm, LogProb, PHREDProb};
use derive_builder::Builder;
use itertools::Itertools;
use itertools_num::linspace;
use rayon::prelude::*;
use rgsl::randist::binomial::binomial_pdf;
use rgsl::randist::poisson::poisson_pdf;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;

use crate::model::modes::tumor::TumorNormalPairView;
use crate::model::AlleleFreq;
use crate::utils;

const MIN_DEPTH: u32 = 10;
const MAX_GAIN: i32 = 21;

pub fn depth_pmf(observed_depth: u32, true_depth: f64) -> LogProb {
    LogProb(poisson_pdf(observed_depth, true_depth).ln())
}

pub fn allele_freq_pdf(
    observed_allele_freq: AlleleFreq,
    true_allele_freq: AlleleFreq,
    depth: u32,
) -> LogProb {
    let k = (*observed_allele_freq * depth as f64).round() as u32;
    LogProb(binomial_pdf(k, *true_allele_freq, depth).ln())
}

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    #[builder(private)]
    bcf_reader: bcf::Reader,
    #[builder(private)]
    bcf_writer: bcf::Writer,
    prior: LogProb,
    purity: f64,
}

impl CallerBuilder {
    pub fn bcfs<P: AsRef<Path>>(
        mut self,
        in_path: Option<P>,
        out_path: Option<P>,
    ) -> Result<Self, Box<Error>> {
        self = self.bcf_reader(if let Some(path) = in_path {
            bcf::Reader::from_path(path)?
        } else {
            bcf::Reader::from_stdin()?
        });

        let bcf_reader = self.bcf_reader.as_ref().unwrap();

        let mut header = bcf::Header::new();
        for sample in bcf_reader.header().samples() {
            header.push_sample(sample);
        }

        header.push_record(
            "##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number in tumor sample\">"
                .as_bytes(),
        );
        header.push_record(
            "##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Subclone fraction affected by \
             the CNV.\">"
                .as_bytes(),
        );
        header.push_record(
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End of copy number variation.\">"
                .as_bytes(),
        );
        header.push_record(
            "##INFO=<ID=LOCI,Number=1,Type=Integer,Description=\"Number of contained loci.\">"
                .as_bytes(),
        );
        header.push_record(
            "##FORMAT=<ID=LOCI_DP,Number=.,Type=Integer,Description=\"Depths of contained loci.\">"
                .as_bytes(),
        );
        header.push_record(
            "##FORMAT=<ID=LOCI_VAF,Number=.,Type=Integer,Description=\"VAFs of contained loci.\">"
                .as_bytes(),
        );

        // register sequences
        for rec in bcf_reader.header().header_records() {
            match rec {
                bcf::header::HeaderRecord::Contig { values, .. } => {
                    let name = values.get("ID").unwrap();
                    let len = values.get("length").unwrap();
                    header.push_record(format!("##contig=<ID={},length={}>", name, len).as_bytes());
                }
                _ => (),
            }
        }

        Ok(self.bcf_writer(if let Some(path) = out_path {
            bcf::Writer::from_path(path, &header, false, false)?
        } else {
            bcf::Writer::from_stdout(&header, false, false)?
        }))
    }
}

impl Caller {
    pub fn call(&mut self) -> Result<(), Box<Error>> {
        let min_prob_germline_het = LogProb(0.95_f64.ln());

        // obtain records
        let mut calls = HashMap::new();
        for record in self.bcf_reader.records() {
            let mut record = record?;
            let call = Call::new(&mut record)?.unwrap();
            if call.prob_germline_het >= min_prob_germline_het && call.depth_normal >= MIN_DEPTH {
                calls
                    .entry(call.rid)
                    .or_insert_with(|| Vec::new())
                    .push(call);
            }
        }

        // normalization
        let mean_depth = |filter: &Fn(&Call) -> u32| {
            calls.values().flatten().map(filter).sum::<u32>() as f64 / calls.len() as f64
        };
        let mean_depth_tumor = mean_depth(&|call: &Call| call.depth_tumor);
        let mean_depth_normal = mean_depth(&|call: &Call| call.depth_normal);
        let depth_norm_factor = mean_depth_tumor / mean_depth_normal;

        let prior = self.prior;
        let purity = self.purity;
        let cnv_calls: BTreeMap<_, _> = calls
            .par_iter()
            .map(|(rid, calls)| {
                let hmm = HMM::new(depth_norm_factor, prior, purity);

                let (states, _prob) = hmm::viterbi(&hmm, calls);

                (
                    *rid,
                    states
                        .iter()
                        .map(|state| hmm.states[**state])
                        .zip(calls.iter())
                        .group_by(|item| item.0)
                        .into_iter()
                        .filter_map(|(cnv, group)| {
                            if cnv.gain == 0 {
                                return None;
                            }
                            let group = group.into_iter().map(|item| item.1).collect_vec();
                            let first_call = group[0];
                            if group.len() > 1 {
                                let last_call = group[group.len() - 1];

                                // calculate posterior probability of no CNV
                                let prob_no_cnv = hmm.prob_no_cnv(&group);

                                Some(CNVCall {
                                    rid: *rid,
                                    pos: first_call.start,
                                    end: last_call.start + 1,
                                    cnv: cnv,
                                    prob_no_cnv,
                                    calls: group,
                                })
                            } else {
                                None
                            }
                        })
                        .collect_vec(),
                )
            })
            .collect();

        let mut record = self.bcf_writer.empty_record();
        for calls in cnv_calls.values() {
            for call in calls {
                call.write(&mut record, depth_norm_factor)?;
                self.bcf_writer.write(&record)?;
            }
        }
        Ok(())
    }
}

pub struct CNVCall<'a> {
    rid: u32,
    pos: u32,
    end: u32,
    cnv: CNV,
    prob_no_cnv: LogProb,
    calls: Vec<&'a Call>,
}

impl<'a> CNVCall<'a> {
    pub fn write(
        &self,
        record: &mut bcf::Record,
        depth_norm_factor: f64,
    ) -> Result<(), Box<Error>> {
        record.set_rid(&Some(self.rid));
        record.set_pos(self.pos as i32);
        record.push_info_integer(b"END", &[self.end as i32])?;
        record.set_alleles(&[b".", b"<CNV>"])?;
        record.push_info_integer(b"CN", &[2 + self.cnv.gain])?;
        record.push_info_float(b"VAF", &[*self.cnv.allele_freq as f32])?;
        record.push_info_integer(b"LOCI", &[self.calls.len() as i32])?;

        let mut loci_dp = Vec::new();
        loci_dp.extend(self.calls.iter().map(|call| call.depth_tumor as i32));
        loci_dp.extend(
            self.calls
                .iter()
                .map(|call| (call.depth_normal as f64 * depth_norm_factor).round() as i32),
        );
        record.push_format_integer(b"LOCI_DP", &loci_dp)?;

        let mut loci_vaf = Vec::new();
        loci_vaf.extend(self.calls.iter().map(|call| *call.allele_freq_tumor as f32));
        loci_vaf.extend(
            self.calls
                .iter()
                .map(|call| *call.allele_freq_normal as f32),
        );
        record.push_format_float(b"LOCI_VAF", &loci_vaf)?;
        record.set_qual(*PHREDProb::from(self.prob_no_cnv) as f32);

        Ok(())
    }
}

pub struct HMM {
    states: Vec<CNV>,
    state_by_gain: HashMap<i32, Vec<hmm::State>>,
    depth_norm_factor: f64,
    prob_cnv: LogProb,
    prob_no_cnv: LogProb,
}

impl HMM {
    fn new(depth_norm_factor: f64, prob_cnv: LogProb, purity: f64) -> Self {
        let mut states = Vec::new();
        let mut state_by_gain = HashMap::new();
        for allele_freq in linspace(0.1, 1.0, 10) {
            for gain in -2..MAX_GAIN {
                if gain != 0 || allele_freq == 1.0 {
                    let cnv = CNV {
                        gain: gain,
                        allele_freq: AlleleFreq(allele_freq),
                        purity,
                    };
                    state_by_gain
                        .entry(gain)
                        .or_insert_with(Vec::new)
                        .push(hmm::State(states.len() - 1));
                    states.push(cnv);
                }
            }
        }
        let prob_no_cnv = prob_cnv.ln_one_minus_exp();
        // we have states.len() different possible CNVs (TODO, this is likely wrong)
        let prob_cnv = prob_cnv - LogProb((states.len() as f64 - 1.0).ln());

        HMM {
            states,
            state_by_gain,
            depth_norm_factor,
            prob_cnv: prob_cnv,
            prob_no_cnv,
        }
    }

    pub fn prob_no_cnv(&self, observations: &[&Call]) -> LogProb {
        let likelihood_no_cnv = likelihood(
            self,
            iter::repeat(self.state_by_gain.get(&0).unwrap()[0]),
            observations.iter().cloned(),
        );
        let mut likelihoods = vec![likelihood_no_cnv];
        for gain in -2..MAX_GAIN {
            if gain != 0 {
                let af_spectrum = self.state_by_gain.get(&gain).unwrap();
                likelihoods.push(LogProb::ln_simpsons_integrate_exp(
                    |i, _| {
                        let state = af_spectrum[i];
                        likelihood(self, iter::repeat(state), observations.iter().cloned())
                    },
                    0.0,
                    1.0,
                    af_spectrum.len() - 1,
                ));
            }
        }

        LogProb::ln_sum_exp(&likelihoods)
    }
}

impl hmm::Model<Call> for HMM {
    fn num_states(&self) -> usize {
        self.states.len()
    }

    fn states(&self) -> hmm::StateIter {
        hmm::StateIter::new(self.num_states())
    }

    fn transitions(&self) -> hmm::StateTransitionIter {
        hmm::StateTransitionIter::new(self.num_states())
    }

    fn transition_prob(&self, from: hmm::State, to: hmm::State) -> LogProb {
        let from = self.states[*from];
        let to = self.states[*to];
        if from == to {
            LogProb(0.5f64.ln()) + self.prob_no_cnv
        } else if to.gain == 0 {
            LogProb(0.5f64.ln()) + self.prob_no_cnv
        } else {
            self.prob_cnv
        }
    }

    fn initial_prob(&self, state: hmm::State) -> LogProb {
        let state = self.states[*state];
        if state.gain == 0 {
            self.prob_no_cnv
        } else {
            self.prob_cnv
        }
    }

    fn observation_prob(&self, state: hmm::State, call: &Call) -> LogProb {
        let cnv = self.states[*state];
        let prob05 = LogProb(0.5f64.ln());

        // handle allele freq changes
        let prob_af = if let Some(alt_af) = cnv.expected_allele_freq_alt_affected() {
            let p = (prob05 + call.prob_allele_freq_tumor(alt_af) + call.prob_germline_het)
                .ln_add_exp(
                    prob05
                        + call.prob_allele_freq_tumor(
                            cnv.expected_allele_freq_ref_affected().unwrap(),
                        ),
                );
            (call.prob_germline_het + p).ln_add_exp(call.prob_germline_het.ln_one_minus_exp())
        } else {
            LogProb::ln_one()
        };

        // handle depth changes
        let prob_depth = call.prob_depth_tumor(
            call.depth_normal as f64 * self.depth_norm_factor * cnv.expected_depth_factor(),
        );

        prob_af + prob_depth
    }
}

pub fn likelihood<'a, O: 'a>(
    hmm: &hmm::Model<O>,
    states: impl IntoIterator<Item = hmm::State>,
    observations: impl Iterator<Item = &'a O>,
) -> LogProb {
    let mut from = None;
    let mut p = LogProb::ln_one();
    for (state, obs) in states.into_iter().zip(observations) {
        if let Some(from) = from {
            p += hmm.transition_prob(from, state);
        } else {
            p = hmm.initial_prob(state);
        }
        p += hmm.observation_prob(state, obs);
        from = Some(state);
    }

    p
}

pub fn marginal<'a, O: 'a>(
    hmm: &hmm::Model<O>,
    observations: impl IntoIterator<Item = &'a O>,
) -> LogProb {
    let mut prev = vec![LogProb::ln_zero(); hmm.num_states()];
    let mut curr = prev.clone();

    for (i, obs) in observations.into_iter().enumerate() {
        for to in hmm.states() {
            let prob_obs = hmm.observation_prob(to, obs);
            curr[*to] = if i == 0 {
                hmm.initial_prob(to)
            } else {
                prob_obs
                    + LogProb::ln_sum_exp(
                        &hmm.states()
                            .map(|from| prev[*from] + hmm.transition_prob(from, to))
                            .collect_vec(),
                    )
            };
        }
        mem::swap(&mut prev, &mut curr);
    }

    LogProb::ln_sum_exp(&prev.into_iter().collect_vec())
        .cap_numerical_overshoot(utils::NUMERICAL_EPSILON)
}

#[derive(Debug)]
pub struct Call {
    prob_germline_het: LogProb,
    allele_freq_tumor: AlleleFreq,
    allele_freq_normal: AlleleFreq,
    depth_tumor: u32,
    depth_normal: u32,
    start: u32,
    rid: u32,
}

impl Call {
    pub fn new(record: &mut bcf::Record) -> Result<Option<Self>, Box<Error>> {
        let prob_germline_het = record.info(b"PROB_GERMLINE_HET").float()?;
        if let Some(prob_germline_het) = prob_germline_het {
            let prob_germline_het = LogProb::from(PHREDProb(prob_germline_het[0] as f64));
            let depths = record
                .format(b"DP")
                .integer()?
                .into_iter()
                .map(|d| d[0] as u32)
                .collect_vec();
            let allele_freqs = record.format(b"AF").float()?;

            Ok(Some(Call {
                allele_freq_tumor: AlleleFreq(allele_freqs.tumor()[0] as f64),
                allele_freq_normal: AlleleFreq(allele_freqs.normal()[0] as f64),
                depth_tumor: *depths.tumor(),
                depth_normal: *depths.normal(),
                prob_germline_het: prob_germline_het,
                start: record.pos(),
                rid: record.rid().unwrap(),
            }))
        } else {
            Ok(None)
        }
    }

    pub fn prob_allele_freq_tumor(&self, true_allele_freq: AlleleFreq) -> LogProb {
        allele_freq_pdf(self.allele_freq_tumor, true_allele_freq, self.depth_tumor)
    }

    pub fn prob_depth_tumor(&self, true_depth: f64) -> LogProb {
        depth_pmf(self.depth_tumor, true_depth)
    }
}

#[derive(PartialEq, Copy, Clone, Debug)]
pub struct CNV {
    gain: i32,
    allele_freq: AlleleFreq,
    purity: f64,
}

impl CNV {
    pub fn expected_allele_freq_alt_affected(&self) -> Option<AlleleFreq> {
        if self.gain > -2 {
            Some(AlleleFreq(
                *self.allele_freq * (1.0 + self.gain as f64) / (2.0 + self.gain as f64)
                    + (1.0 - *self.allele_freq) * 0.5,
            ))
        } else if self.purity < 1.0 {
            // gain = -2: all lost in tumor cells, hence 100% normal cells at this locus.
            // Therefore VAF=0.5.
            Some(AlleleFreq(0.5))
        } else {
            None
        }
    }

    pub fn expected_allele_freq_ref_affected(&self) -> Option<AlleleFreq> {
        self.expected_allele_freq_alt_affected()
            .map(|af| AlleleFreq(1.0) - af)
    }

    pub fn expected_depth_factor(&self) -> f64 {
        self.purity * (*self.allele_freq * (2.0 + self.gain as f64) / 2.0 + 1.0 - *self.allele_freq)
            + (1.0 - self.purity)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_allele_freq_pdf() {
        assert_eq!(
            allele_freq_pdf(AlleleFreq(0.64), AlleleFreq(1.0), 10),
            LogProb::ln_zero()
        );
        assert_eq!(
            allele_freq_pdf(AlleleFreq(0.1), AlleleFreq(0.0), 10),
            LogProb::ln_zero()
        );
    }
}
