use std::path::{PathBuf, Path};
use std::error::Error;

use bio::stats::{PHREDProb, LogProb};
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
        // TODO read regions from BED, fetch records from bcf_reader, write results
    }
}


pub struct Gain {
    value: f64,
    prob: LogProb,
}


pub fn prob_gains(record: &bcf::Record) -> Result<Vec<Gain>, Box<Error>> {
    let logprob = |p| LogProb::from(PHREDProb(p as f64));
    let prob_afs_het = record.format(b"PROB_AFS_HET").float()?;
    let prob_nocov = record.format(b"PROB_NOCOV").float()?;

    let mut gains = Vec::new();

    // partial loss or gain
    for prob_germline_het in record.info(b"PROB_GERMLINE_HET").float()?.unwrap() {
        let prob_germline_het = logprob(*prob_germline_het);

        for i in 0..AFS.len() {
            let tumor_af = *AFS[i];
            let gain = (2.0 * tumor_af - 1.0) / (1.0 - tumor_af);
            assert!(gain > -1.0);

            let prob_gain_alt = logprob(prob_afs_het.tumor()[i]);
            let prob_gain_ref = logprob(prob_afs_het.tumor()[AFS.len() - 1 - i]);

            let prob = prob_germline_het + LogProb::ln_sum_exp(&[prob05 + prob_gain_alt, prob05 + prob_gain_ref]) + prob_germline_het.ln_one_minus_exp();

            gains.push(Gain {
                value: gain,
                prob: prob,
            })
        }
    }

    // complete loss (gain=-2)
    gains.push(Gain {
        value: -2.0,
        prob: LogProb::ln_sum_exp(&[
                logprob(prob_nocov.tumor()[0]) + logprob(prob_nocov.normal()[0]).ln_one_minus_exp(),
                logprob(prob_nocov.normal()[0])
            ])
        });

    Ok(gains)
}

const prob05: LogProb = LogProb(0.5f64.ln());
