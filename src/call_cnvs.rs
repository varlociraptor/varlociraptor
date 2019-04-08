use std::error::Error;
use std::path::Path;
use std::collections::BTreeMap;

use bio::stats::{hmm, LogProb, PHREDProb};
use derive_builder::Builder;
use itertools::Itertools;
use itertools_num::linspace;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;
use rgsl::randist::poisson::poisson_pdf;

use crate::model::modes::tumor::TumorNormalPairView;
use crate::model::AlleleFreq;

lazy_static! {
    // Any changes here would lead to incompatibilities when reading in older VCFs.
    // Hence be careful!
    pub static ref AFS: Vec<AlleleFreq> = linspace(0.0, 1.0, 11)
        .map(|af| AlleleFreq(af))
        .collect_vec();
}

pub fn depth_pmf(depth: u32, mean: f64) -> LogProb {
    LogProb(poisson_pdf(depth, mean))
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
        let mean_depth_tumor = calls.iter().map(|call| call.depth_tumor).mean();
        let mean_depth_normal = calls.iter().map(|call| call.depth_normal).mean();
        let depth_norm_factor = mean_depth_tumor / mean_depth_normal;

        for (rid, calls) in calls.into_iter().group_by(|call| call.rid).into_iter() {
            let hmm = HMM::new(depth_norm_factor);
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
                record.push_info_integer(b"CN", &[2 + cnv.gain]);
                record.push_info_float(b"VAF", &[*cnv.allele_freq as f32]);

                self.bcf_writer.write(&record)?;
            }
        }
        Ok(())
    }
}

pub struct HMM {
    states: Vec<CNV>,
    depth_norm_factor: f64,
}

impl HMM {
    fn new(depth_norm_factor: f64) -> Self {
        let mut states = Vec::new();
        for allele_freq in linspace(0.0, 1.0, 10) {
            for gain in 0..20 {
                states.push(CNV { gain: gain, allele_freq: AlleleFreq(allele_freq)});
            }
        }

        HMM { states, depth_norm_factor }
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

        // handle allele freq changes
        let prob_af = LogProb::ln_sum_exp(&[
            prob05 + call.prob_allele_freq_tumor(cnv.expected_allele_freq_alt_affected()),
            prob05 + call.prob_allele_freq_tumor(cnv.expected_allele_freq_ref_affected()),
            call.prob_germline_not_het
        ]);

        // handle depth changes
        let prob_depth = call.prob_depth_tumor(call.depth_normal * self.depth_norm_factor * cnv.expected_depth_factor());

        prob_af + prob_depth
    }
}

pub struct Call {
    prob_afs_het: Vec<LogProb>,
    prob_nocov_tumor: LogProb,
    prob_nocov_normal: LogProb,
    prob_germline_not_het: LogProb,
    depth_tumor: u32,
    depth_normal: u32,
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

    pub fn prob_allele_freq_tumor(&self, allele_freq: AlleleFreq) -> LogProb {
        match AFS.binary_search(&allele_freq) {
            Ok(i) => self.prob_afs_het[i],
            Err(i) => {
                // need to interpolate
                let p0 = *self.prob_afs_het[i - 1];
                let p1 = *self.prob_afs_het[i];
                let af0 = *AFS[i - 1];
                let af1 = *AFS[i];

                LogProb(p0 * (af1 - *allele_freq) / (af1 - af0) + p1 * (*allele_freq - af0) / (af1 - af0))
            }
        }
    }

    pub fn prob_depth_tumor(&self, depth: f64) -> LogProb {
        depth_pmf(self.depth_tumor, depth)
    }
}

#[derive(PartialEq, Copy, Clone, Debug)]
pub struct CNV {
    gain: i32,
    allele_freq: AlleleFreq,
}


impl CNV {
    pub fn expected_allele_freq_alt_affected(&self) -> AlleleFreq {
        AlleleFreq(*self.allele_freq * (1.0 + self.gain as f64) / (2.0 + self.gain as f64) + (1.0 - *self.allele_freq) * 0.5)
    }

    pub fn expected_allele_freq_ref_affected(&self) -> AlleleFreq {
        AlleleFreq(1.0) - self.expected_allele_freq_alt_affected()
    }

    pub fn expected_depth_factor(&self) -> f64 {
        (*self.allele_freq * (2.0 + self.gain as f64) / 2.0 + 1.0 - *self.allele_freq)
    }
}
