use std::str;

use anyhow::Result;
use rust_htslib::bcf::{self, Read};
use serde_json::{json, Value};

use crate::errors;

#[derive(Debug, Clone, Serialize)]
struct ScatterProbs {
    chromosome: String,
    prob: f32,
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
            let pos = rec.pos() - rec.rlen();
            let prob = rec.info(b"PROB_PRESENT").float()?.unwrap()[0];
            plot_data.push(ScatterProbs {
                chromosome: chr,
                prob,
                pos,
            });
        }
    }
    if plot_data.is_empty() {
        return Err(errors::Error::NoRecordsFound.into());
    }

    let print_plot =
        |data: serde_json::Value, blueprint: &str| -> Result<()> {
            let mut blueprint = serde_json::from_str(blueprint)?;
            if let Value::Object(ref mut blueprint) = blueprint {
                blueprint["data"][0]["values"] = data;
                println!("{}", serde_json::to_string_pretty(blueprint)?);
                Ok(())
            } else {
                unreachable!();
            }
        };

    print_plot(
        json!(plot_data),
        //serde_json::value::Value::String(tumor.to_string()),
        include_str!("../../templates/plots/probabilities_scatter_contour.json"),
    )
}
