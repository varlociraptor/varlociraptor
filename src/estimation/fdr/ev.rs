//! This module implements fully bayesian FDR control as presented by
//! Müller, Parmigiani, and Rice, "FDR and Bayesian Multiple Comparisons Rules" (July 2006).
//! Johns Hopkin's University, Dept. of Biostatistics Working Papers. Working Paper 115.
//! Basically, the expected FDR is calculated directly from the posterior error probabilities.

use std::error::Error;
use std::io;

use bio::stats::{LogProb, PHREDProb, bayesian};
use itertools::Itertools;
use rust_htslib::bcf;
use csv;

use Event;
use estimation::fdr::{Record, ALPHAS};
use model;
use utils;


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
    writer: &mut W,
    event: &E,
    vartype: &model::VariantType) -> Result<(), Box<Error>> {
    let mut writer = csv::Writer::from_writer(writer).delimiter(b'\t');
    try!(writer.write(["FDR", "max-prob"].into_iter()));

    let prob_dist = utils::collect_prob_dist(calls, event, vartype)?;

    if prob_dist.is_empty() {
        for &alpha in &ALPHAS {
            writer.write([&format!("{}", alpha), ""].iter())?;
        }
        return Ok(());
    }

    // estimate FDR
    let pep_dist = prob_dist.into_iter().map(|p| LogProb(*p).ln_one_minus_exp()).collect_vec();
    let fdrs = bayesian::expected_fdr(&pep_dist);
    debug!("FDRs {:?}", &fdrs[..200]);
    debug!("PEPs {:?}", &pep_dist[..200]);

    for &alpha in ALPHAS.iter().rev() {
        let ln_alpha = LogProb(alpha.ln());
        // find the largest pep for which fdr <= alpha
        // do not let peps with the same value cross the boundary
        for i in (0..fdrs.len()).rev() {
            if fdrs[i] <= ln_alpha && (i == 0 || pep_dist[i] != pep_dist[i - 1]) {
                writer.encode(&Record { alpha: alpha, gamma: PHREDProb::from(pep_dist[i]) })?;
                break;
            }
        }
    }
    writer.flush()?;

    Ok(())
}
