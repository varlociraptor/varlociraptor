use crate::errors;
use anyhow::Result;
use bcf::header::HeaderRecord;
use rust_htslib::bcf::{self, Read};
use serde_json::{json, Value};
use std::str;

use bio::stats::{PHREDProb, Prob};
#[derive(Debug, Clone, Serialize)]
struct ScatterProbs {
    chromosome: String,
    event: String,
    prob: f32,
    pos: i64,
}

pub(crate) fn variant_probabilities() -> Result<()> {
    let mut bcf = bcf::Reader::from_stdin()?;

    let mut plot_data = Vec::new();
    let header = bcf.header().to_owned();
    let events: Vec<String> = header
        .header_records()
        .into_iter()
        .filter_map(|x| match x {
            HeaderRecord::Info { key: _, values } => Some(values["ID"].clone()),
            _ => None,
        })
        .filter(|x| x.starts_with("PROB_"))
        .collect();

    for record in bcf.records() {
        let rec = record.unwrap();
        if let Some(rid) = rec.rid() {
            for event in events.iter() {
                let chr = String::from_utf8(header.rid2name(rid)?.to_vec())?;
                let pos = rec.pos() - rec.rlen();
                let phred_prob = rec.info(event.as_bytes()).float()?.unwrap()[0];
                let prob = *Prob::from(PHREDProb(phred_prob as f64)) as f32;
                if prob > 0.0 {
                    plot_data.push(ScatterProbs {
                        chromosome: chr,
                        event: event[5..].to_string(),
                        prob,
                        pos,
                    });
                }
            }
        }
    }
    if plot_data.is_empty() {
        return Err(errors::Error::NoRecordsFound.into());
    }

    let print_plot = |data: serde_json::Value, blueprint: &str| -> Result<()> {
        let mut blueprint = serde_json::from_str(blueprint)?;
        if let Value::Object(ref mut blueprint) = blueprint {
            blueprint["data"]["values"] = data;
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
