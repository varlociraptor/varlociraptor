// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! This module implements fully bayesian FDR control as presented by
//! Müller, Parmigiani, and Rice, "FDR and Bayesian Multiple Comparisons Rules" (July 2006).
//! Johns Hopkin's University, Dept. of Biostatistics Working Papers. Working Paper 115.
//! Basically, the expected FDR is calculated directly from the posterior error probabilities.

use std::error::Error;
use std::path::Path;

use bio::stats::{bayesian, LogProb};
use itertools::Itertools;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;

use crate::model;
use crate::utils;
use crate::Event;
use crate::utils::is_phred_scaled;

/// Print thresholds to control FDR of given calls at multiple levels.
///
/// # Arguments
///
/// * `inbcf` - path to BCF with varlociraptor calls (None for stdin)
/// * `outbcf` - path to BCF with filtered varlociraptor calls (None for stdout)
/// * `null_calls` - calls under the null model, e.g. obtained by swapping tumor and normal sample
/// * `writer` - writer for resulting thresholds
/// * `events` - the set of events to control (sum of the probabilities of the individual events at a site)
/// * `vartype` - the variant type to consider
/// * `alpha` - the FDR threshold to control for
pub fn control_fdr<E: Event, R, W>(
    inbcf: R,
    outbcf: Option<W>,
    events: &[E],
    vartype: &model::VariantType,
    alpha: LogProb,
) -> Result<(), Box<Error>>
where
    R: AsRef<Path>,
    W: AsRef<Path>,
{
    // first pass on bcf file
    let mut inbcf_reader = bcf::Reader::from_path(&inbcf)?;

    if !is_phred_scaled(&inbcf_reader) {
        panic!("Event probabilities are not PHRED scaled, aborting.")
    }

    // setup output file
    let header = bcf::Header::from_template(inbcf_reader.header());
    let mut outbcf = match outbcf {
        Some(p) => bcf::Writer::from_path(p, &header, false, false)?,
        None => bcf::Writer::from_stdout(&header, false, false)?,
    };

    let mut threshold = None;

    if alpha != LogProb::ln_one() {
        // do not filter by FDR if alpha is 1.0
        // TODO: remove hits where another event has a higher probability
        // Otherwise, if there are just enough calls, events like PROB_SOMATIC=8, PROB_ABSENT=2
        // can end up in the filtered results.
        let prob_dist = utils::collect_prob_dist(&mut inbcf_reader, events, vartype)?
            .into_iter()
            .rev()
            .map(|p| LogProb(*p))
            .collect_vec();

        // estimate FDR
        let pep_dist = prob_dist.iter().map(|p| p.ln_one_minus_exp()).collect_vec();
        let fdrs = bayesian::expected_fdr(&pep_dist);

        if fdrs.is_empty() {
            threshold = None;
        } else if fdrs[0] > alpha {
            threshold = Some(LogProb::ln_one());
        } else {
            // find the largest pep for which fdr <= alpha
            // do not let peps with the same value cross the boundary
            for i in (0..fdrs.len()).rev() {
                if fdrs[i] <= alpha && (i == 0 || pep_dist[i] != pep_dist[i - 1]) {
                    let prob = prob_dist[i];

                    threshold = Some(prob);
                    break;
                }
            }
        }
    }

    // second pass on bcf file
    let mut inbcf_reader = bcf::Reader::from_path(&inbcf)?;
    utils::filter_by_threshold(&mut inbcf_reader, threshold, &mut outbcf, events, vartype)?;

    Ok(())
}
