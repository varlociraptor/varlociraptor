use std::path::{PathBuf, Path};
use std::error::Error;

use bio::stats::{PHREDProb, LogProb, hmm};
use derive_builder::Builder;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;
use itertools::Itertools;
use itertools_num::linspace;

use crate::model::AlleleFreq;
use crate::model::modes::tumor::TumorNormalPairView;

lazy_static! {
    static ref AFS: Vec<AlleleFreq> = linspace(0.0, 1.0, 11).map(|af| AlleleFreq(af)).collect_vec();
}

const prob05: LogProb = LogProb(0.5f64.ln());

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Caller {
    regions: PathBuf,
    #[builder(private)]
    bcf_reader: bcf::Reader,
    #[builder(private)]
    bcf_writer: bcf::Writer,
}

impl CallerBuilder {
    pub fn bcfs<P: AsRef<Path>>(mut self, in_path: Option<P>, out_path: Option<P>) -> Result<Self, Box<Error>> {
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
            "##INFO=<ID=CN,Number=1,Type=Float,Description=\"Copy number gain or loss in tumor sample \
             (negative=loss, positive=gain). Fractional because CNV can be subclonal.\">".as_bytes()
        );
        header.push_record(
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End of copy number variation.\">".as_bytes()
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
        let calls = Vec::new();
        for record in self.bcf_reader.records() {
            let record = record?;
            calls.push(Call::new(&record)?);
        }

        let hmm = HMM::new();

        let (states, prob) = hmm::viterbi(&hmm, calls);

        let mut record = self.bcf_writer.empty_record();

        for (gain, group) in states.iter().map(|s| hmm.states[s]).zip(&calls).group_by(|item| item.0) {
            let pos = group.iter().first().1.start;
            let end = group.iter().last().1.start + 1;
            record.set_pos(pos);
            record.push_info_integer(b"END", &[end])?;
            record.push_info_float(b"CN", &[2.0 + gain])?;
            record.set_alleles(&[".", "<CNV>"])?;
        }


    }
}

pub struct HMM {
    states: Vec<Gain>,
    start_state: State
}

impl HMM {
    fn new() -> Self {
        let states = AFS.iter().map(|tumor_af| Gain::from(tumor_af)).collect_vec();
        let start_state = Gain::from(AlleleFreq(0.5));

        HMM{states, start_state}
    }
}

impl Model<Call> for HMM {
    fn num_states(&self) -> usize {
        self.states.len()
    }

    fn states(&self) -> StateIter {
        StateIter::new(self.num_states())
    }

    fn transitions(&self) -> StateTransitionIter {
        StateTransitionIter::new(self.num_states())
    }

    fn transition_prob(&self, from: State, to: State) -> LogProb {
        LogProb(0.0001_f64.ln())
    }

    fn initial_prob(&self, state: State) -> LogProb {
        if state == self.start_state {
            LogProb::ln_one()
        }
        else {
            LogProb::ln_zero()
        }
    }

    fn observation_prob(&self, state: State, call: &Call) -> LogProb {

        let gain = self.states[state];

        if gain > -2.0 {
                let prob_gain_alt = call.prob_afs_het[state];
                let prob_gain_ref = call.prob_afs_het[self.num_states() - 1 - state];
                LogProb::ln_sum_exp(&[prob05 + prob_gain_alt, prob05 + prob_gain_ref, call.prob_germline_het.ln_one_minus_exp()])
        }

        else {
                call.prob_nocov_normal.ln_add_exp(call.prob_nocov_tumor + call.prob_nocov_normal.ln_one_minus_exp())
        }

    }

}


pub struct Call {
    prob_afs_het: Vec<LogProb>,
    prob_nocov_tumor: LogProb,
    prob_nocov_normal: LogProb,
    prob_germline_not_het: LogProb,
    start: u32,
    rid: u32
}

impl Call {
    pub fn new(record: &bcf::Record) -> Result<Option<Self>, Box<Error>> {
        let prob_germline_het = record.info(b"PROB_GERMLINE_HET").float()?;
        if let Some(prob_germline_het) = prob_germline_het {
            let logprob = |p| LogProb::from(PHREDProb(p as f64));
            let prob_afs_het = record.format(b"PROB_AFS_HET").float()?;
            let prob_nocov = record.format(b"PROB_NOCOV").float()?;


            Some(Call {
                prob_afs_het: prob_afs_het.tumor().into_iter().map(&logprob).collect(),
                prob_nocov_tumor: logprob(prob_nocov.tumor()),
                prob_nocov_normal: logprob(prob_nocov.normal()),
                prob_germline_not_het: logprob(prob_germline_het).ln_one_minus_exp(),
                start: record.pos(),
                rid: record.rid()
            })
        } else {
            None
        }


    }
}

pub struct Gain(f64);

impl From<AlleleFreq> for Gain {
    fn from(af: AlleleFreq) -> Gain {
        Gain((2.0 * af - 1.0) / (1.0 - af))
    }
}



