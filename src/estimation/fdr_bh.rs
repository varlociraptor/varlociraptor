use std::error::Error;
use std::io;

use csv;
use itertools::Itertools;
use rust_htslib::bcf;
use bio::stats::{LogProb, PHREDProb, Prob};
use ordered_float::NotNaN;

use utils;
use Event;
use model;


const ALPHAS: [f64; 21] = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.25];


fn collect_dist<E: Event>(calls: &bcf::Reader, event: &E, vartype: &model::VariantType) -> Result<Vec<NotNaN<f64>>, Box<Error>> {
    let mut record = bcf::Record::new();
    let mut prob_dist = Vec::new();
    loop {
        if let Err(e) = calls.read(&mut record) {
            if e.is_eof() {
                break;
            } else {
                return Err(Box::new(e));
            }
        }

        let variants = try!(utils::collect_variants(&mut record, false, false, None));
        let tag = event.tag_name("PROB");
        let event_probs = try!(record.info(tag.as_bytes()).float());
        if let Some(event_probs) = event_probs {
            // tag present
            for (variant, event_prob) in variants.into_iter().zip(event_probs.into_iter()) {
                if let Some(variant) = variant {
                    if !variant.is_type(vartype) || event_prob.is_nan() {
                        continue;
                    }
                    let event_prob = LogProb::from(PHREDProb(*event_prob as f64));
                    prob_dist.push(try!(NotNaN::new(*event_prob)));
                }
            }
        }
    }
    prob_dist.sort();
    Ok(prob_dist)
}


fn pval(x: NotNaN<f64>, dist: &[NotNaN<f64>]) -> LogProb {
    let i = dist.binary_search(&x).unwrap_or_else(|i| i);
    let f0 = LogProb::from(Prob(i as f64 / dist.len() as f64));
    f0.ln_one_minus_exp()
}


#[derive(RustcEncodable)]
struct Record {
    alpha: f64,
    gamma: PHREDProb
}


/// Print thresholds to control FDR of given calls at multiple levels.
///
/// # Arguments
///
/// * `calls` - BCF reader with prosic calls
/// * `null_calls` - calls under the null model, e.g. obtained by swapping tumor and normal sample
/// * `writer` - writer for resulting thresholds
/// * `event` - the event to control
/// * `vartype` - the variant type to consider
pub fn control_fdr<E: Event, W: io::Write>(
    calls: &mut bcf::Reader,
    null_calls: &mut bcf::Reader,
    writer: &mut W,
    event: &E,
    vartype: &model::VariantType) -> Result<(), Box<Error>> {
    let mut writer = csv::Writer::from_writer(writer).delimiter(b'\t');

    let null_dist = try!(collect_dist(null_calls, event, vartype));
    let prob_dist = try!(collect_dist(calls, event, vartype));
    debug!("{} observations in null distribution.", null_dist.len());
    debug!("{} observations in call distribution.", prob_dist.len());
    let pvals = prob_dist.iter().map(|&p| pval(p, &null_dist)).collect_vec();
    let m = pvals.len() as f64;
    let mk_pvals = pvals.iter().enumerate().map(|(k, &p)| (*p) + m.ln() - (m - k as f64 + 1.0).ln()).collect_vec(); // p * m / (m - k + 1)

    try!(writer.write(["FDR", "max-prob"].into_iter()));

    for &alpha in &ALPHAS {
        let mut record = Record { alpha: alpha, gamma: PHREDProb::from(Prob(1.0)) };
        for (&mkp, &event_prob) in mk_pvals.iter().zip(prob_dist.iter()) {
            // the pvalues will be monotolically decreasing
            if mkp <= alpha {
                // the highest pval that is below the threshold
                // record corresponding event probability
                record.gamma = PHREDProb::from(LogProb(*event_prob));
                break;
            }
        }
        try!(writer.encode(&record));
    }
    try!(writer.flush());

    Ok(())
}
