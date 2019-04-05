use std::error::Error;
use std::path::Path;
use std::collections::BTreeMap;

use bio::stats::{hmm, LogProb, PHREDProb};
use derive_builder::Builder;
use itertools::Itertools;
use itertools_num::linspace;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;

use crate::model::modes::tumor::TumorNormalPairView;
use crate::model::AlleleFreq;

lazy_static! {
    // Any changes here would lead to incompatibilities when reading in older VCFs.
    // Hence be careful!
    pub static ref AFS: Vec<AlleleFreq> = linspace(0.0, 1.0, 11)
        .map(|af| AlleleFreq(af))
        .collect_vec();
}

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    #[builder(private)]
    bcf_reader: bcf::Reader,
    #[builder(private)]
    bcf_writer: bcf::Writer,
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
            "##INFO=<ID=HET_DEVIATION,Number=1,Type=Float,Description=\"Systematic deviation from \
             VAF=0.5. Can be either a subclonal loss or gain.\">"
                .as_bytes(),
        );
        header.push_record(
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End of copy number variation.\">"
                .as_bytes(),
        );

        Ok(self.bcf_writer(if let Some(path) = out_path {
            bcf::Writer::from_path(path, &header, false, false)?
        } else {
            bcf::Writer::from_stdout(&header, false, false)?
        }))
    }
}

impl Caller {
    pub fn call(&mut self) -> Result<(), Box<Error>> {
        let mut calls = Vec::new();
        for record in self.bcf_reader.records() {
            let mut record = record?;
            calls.push(Call::new(&mut record)?.unwrap());
        }

        for (rid, calls) in calls.into_iter().group_by(|call| call.rid).into_iter() {
            let hmm = HMM::new();
            let calls = calls.into_iter().collect_vec();

            let (states, _prob) = hmm::viterbi(&hmm, &calls);

            let mut record = self.bcf_writer.empty_record();

            for (cnv, group) in states
                .iter()
                .map(|s| hmm.states[**s])
                .zip(&calls)
                .group_by(|item| item.0)
                .into_iter()
            {
                let mut group = group.into_iter();
                let first_call = group.next().unwrap().1;
                let pos = first_call.start;
                let end = group.last().unwrap().1.start + 1;
                record.set_rid(&Some(rid));
                record.set_pos(pos as i32);
                record.push_info_integer(b"END", &[end as i32])?;
                record.set_alleles(&[b".", b"<CNV>"])?;
                match cnv {
                    CNV::Loss => {
                        record.push_info_integer(b"CN", &[0])?;
                    },
                    CNV::LOH => {
                        record.push_info_integer(b"CN", &[1])?;
                    },
                    CNV::Partial(af_delta) => {
                        record.push_info_float(b"HET_DEVIATION", &[*af_delta as f32])?;
                    }
                }



                self.bcf_writer.write(&record)?;
            }
        }
        Ok(())
    }
}

pub struct HMM {
    states: Vec<CNV>,
    afs_idx: BTreeMap<AlleleFreq, usize>
}

impl HMM {
    fn new() -> Self {
        let states = vec![CNV::Loss, CNV::LOH, CNV::Partial(AlleleFreq(0.1)), CNV::Partial(AlleleFreq(0.2)), CNV::Partial(AlleleFreq(0.3)), CNV::Partial(AlleleFreq(0.4))];
        let afs_idx = AFS.iter().enumerate().map(|(i, af)| (*af, i)).collect();

        HMM { states, afs_idx }
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

    fn transition_prob(&self, _from: hmm::State, _to: hmm::State) -> LogProb {
        LogProb(0.0001_f64.ln())
    }

    fn initial_prob(&self, _state: hmm::State) -> LogProb {
        LogProb((1.0 / self.num_states() as f64).ln())
    }

    fn observation_prob(&self, state: hmm::State, call: &Call) -> LogProb {
        let cnv = self.states[*state];
        let prob05 = LogProb(0.5f64.ln());
        match cnv {
            CNV::Loss => {
                // no coverage in tumor, coverage in normal
                call.prob_nocov_normal
                    .ln_add_exp(call.prob_nocov_tumor + call.prob_nocov_normal.ln_one_minus_exp())
            },
            CNV::LOH => {
                // complete loss of heterozygosity in tumor
                let prob_loss_alt = call.prob_afs_het[0];
                let prob_loss_ref = call.prob_afs_het[AFS.len() - 1];
                LogProb::ln_sum_exp(&[
                    prob05 + prob_loss_alt,
                    prob05 + prob_loss_ref,
                    call.prob_germline_not_het,
                ])
            },
            CNV::Partial(af_delta) => {
                let prob = |af| call.prob_afs_het[*self.afs_idx.get(&af).unwrap()];
                let het = AlleleFreq(0.5);
                let prob_gain_alt = prob(het + af_delta);
                let prob_gain_ref = prob(het - af_delta);
                LogProb::ln_sum_exp(&[
                    prob05 + prob_gain_alt,
                    prob05 + prob_gain_ref,
                    call.prob_germline_not_het,
                ])
            }
        }
    }
}

pub struct Call {
    prob_afs_het: Vec<LogProb>,
    prob_nocov_tumor: LogProb,
    prob_nocov_normal: LogProb,
    prob_germline_not_het: LogProb,
    start: u32,
    rid: u32,
}

impl Call {
    pub fn new(record: &mut bcf::Record) -> Result<Option<Self>, Box<Error>> {
        let prob_germline_het = record.info(b"PROB_GERMLINE_HET").float()?;
        if let Some(prob_germline_het) = prob_germline_het {
            let logprob = |p: f32| LogProb::from(PHREDProb(p as f64));
            let prob_germline_het = logprob(prob_germline_het[0]);
            let prob_afs_het = record
                .format(b"PROB_AFS_HET")
                .float()?
                .tumor()
                .into_iter()
                .cloned()
                .map(&logprob)
                .collect();
            let prob_nocov = record.format(b"PROB_NOCOV").float()?;

            Ok(Some(Call {
                prob_afs_het,
                prob_nocov_tumor: logprob(prob_nocov.tumor()[0]),
                prob_nocov_normal: logprob(prob_nocov.normal()[0]),
                prob_germline_not_het: prob_germline_het.ln_one_minus_exp(),
                start: record.pos(),
                rid: record.rid().unwrap(),
            }))
        } else {
            Ok(None)
        }
    }
}

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum CNV {
    Loss,
    LOH,
    Partial(AlleleFreq),
}
