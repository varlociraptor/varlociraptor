// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! This module implements fully bayesian FDR control as presented by
//! Müller, Parmigiani, and Rice, "FDR and Bayesian Multiple Comparisons Rules" (July 2006).
//! Johns Hopkin's University, Dept. of Biostatistics Working Papers. Working Paper 115.
//! Basically, the expected FDR is calculated directly from the posterior error probabilities.

use std::path::Path;

use anyhow::{bail, Result};
use bio::stats::{bayesian, LogProb};
use itertools::Itertools;
use ordered_float::NotNan;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;

use crate::utils::{event_to_tag, is_phred_scaled};
use crate::variants::model;
use crate::{errors, Event};
use crate::{utils, SimpleEvent};

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
pub fn control_fdr<R, W>(
    inbcf: R,
    outbcf: Option<W>,
    events: &[SimpleEvent],
    vartype: Option<&model::VariantType>,
    alpha: LogProb,
    local: bool,
    smart: bool,
) -> Result<()>
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
        Some(p) => bcf::Writer::from_path(p, &header, false, bcf::Format::Bcf)?,
        None => bcf::Writer::from_stdout(&header, false, bcf::Format::Bcf)?,
    };

    let mut threshold = None;

    // skip events that are not present in the header
    let cleaned_events = events
        .into_iter()
        .filter(|event| {
            inbcf_reader
                .header()
                .info_type(event_to_tag(*event).as_bytes())
                .is_ok()
        })
        .cloned()
        .collect_vec();
    if cleaned_events.is_empty() {
        bail!(errors::Error::InvalidFDRControlEvents)
    } else if cleaned_events != events {
        warn!(
            "Skipping events not present in BCF header: {}",
            events
                .iter()
                .filter_map(|event| if !cleaned_events.contains(event) {
                    Some(event.name())
                } else {
                    None
                })
                .join(", ")
        );
    }

    if local {
        threshold = Some(alpha.ln_one_minus_exp());
    } else if alpha != LogProb::ln_one() {
        let dist_events = if smart {
            vec![SimpleEvent::new("ABSENT"), SimpleEvent::new("ARTIFACT")]
        } else {
            cleaned_events.clone()
        };

        // do not filter by FDR if alpha is 1.0
        // TODO: remove hits where another event has a higher probability
        // Otherwise, if there are just enough calls, events like PROB_SOMATIC=8, PROB_ABSENT=2
        // can end up in the filtered results.
        let convert_prob = if smart {
            |p: NotNan<f64>| LogProb(*p).ln_one_minus_exp()
        } else {
            |p: NotNan<f64>| LogProb(*p)
        };

        let prob_dist = utils::collect_prob_dist(&mut inbcf_reader, &dist_events, vartype)?
            .into_iter()
            .rev()
            .map(convert_prob)
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

        // second pass on bcf file
        inbcf_reader = bcf::Reader::from_path(&inbcf)?;
    }

    utils::filter_by_threshold(
        &mut inbcf_reader,
        threshold,
        &mut outbcf,
        &cleaned_events,
        vartype,
        smart,
    )?;

    Ok(())
}
