use std::collections::BTreeMap;
use std::str;

use anyhow::Result;
use rust_htslib::bcf::{self, Read};
use serde_json::{json, Value};

use crate::errors;
use crate::variants::model::AlleleFreq;

#[derive(Debug, Clone, Serialize)]
struct ScatterProbs {
    chromosome: String,
    log_prob: f32,
    pos: i64,
}

pub(crate) fn variant_probabilities() -> Result<()> {
    let mut bcf = bcf::Reader::from_stdin()?;

    let mut plot_data = Vec::new();
    let header = bcf.header().to_owned();
    for record in bcf.records() {
        let rec = record.unwrap();
        if let Some(rid) = rec.rid() {
            let chr = String::from_utf8(header.rid2name(rid)?.to_vec())?;
            let pos = rec.pos();
            let log_prob = rec.info(b"PROB_PRESENT").float()?.unwrap()[0].log10();
            plot_data.push(ScatterProbs {
                chromosome: chr,
                log_prob,
                pos,
            });
        }
    }
    dbg!(plot_data);
    Ok(())
}
