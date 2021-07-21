use std::collections::BTreeMap;
use std::str;

use anyhow::Result;
use bio::stats::{LogProb, PHREDProb};
use rust_htslib::bcf::{self, Read};
use serde_json::{json, Value};

use crate::errors;
use crate::variants::model::AlleleFreq;
use crate::{Event, SimpleEvent};

#[derive(Debug, Clone, Serialize)]
struct Scattervaf {
    sample: String,
    normal_vaf: f64,
    tumor_vaf: f64,
}

pub(crate) fn vaf_scatter(
    mutational_events: &[String],
    sample_x: &str,
    sample_y: &[String],
) -> Result<()> {
    let mut bcf = bcf::Reader::from_stdin()?;
    let header = bcf.header().to_owned();

    let mut plot_data = Vec::new();

    let id_x = bcf
        .header()
        .sample_id(sample_x.as_bytes())
        .unwrap_or_else(|| panic!("Sample {} not found", sample_x));

    let mut ids_y = BTreeMap::new();

    for s in sample_y {
        let id_y = bcf
            .header()
            .sample_id(s.as_bytes())
            .unwrap_or_else(|| panic!("Sample {} not found", s));
        ids_y.insert(s, id_y);
    }

    for record in bcf.records() {
        let mut rec = record.unwrap();
        let contig = str::from_utf8(header.rid2name(rec.rid().unwrap()).unwrap())?;
        let vcfpos = rec.pos() + 1;

        // obtain VAF estimates (do it here already to work around a segfault in htslib)
        let x_vafs = rec.format(b"AF").float()?[id_x].to_owned();

        let alt_allele_count = (rec.allele_count() - 1) as usize;

        // collect allele probabilities for given events
        let mut allele_probs = vec![LogProb::ln_zero(); alt_allele_count];
        for e in mutational_events {
            let e = SimpleEvent { name: e.to_owned() };
            let tag_name = e.tag_name("PROB");
            if let Some(probs) = rec.info(tag_name.as_bytes()).float()? {
                for i in 0..alt_allele_count {
                    allele_probs[i] =
                        allele_probs[i].ln_add_exp(LogProb::from(PHREDProb(probs[i] as f64)));
                }
            } else {
                info!(
                    "Skipping variant {}:{} because it does not contain the required INFO tag {}.",
                    contig, vcfpos, tag_name
                );
                continue;
            }
        }

        // push into MB function

        for i in 0..alt_allele_count {
            for (y, id_y) in &ids_y {
                let y_vafs = rec.format(b"AF").float()?[*id_y].to_owned();
                // if all alt_alleles are NaN, the list will only contain one NaN, so check for size
                if (i == x_vafs.len() && x_vafs[0].is_nan())
                    || (i == y_vafs.len() && y_vafs[0].is_nan())
                {
                    continue;
                }
                let x_vaf = x_vafs[i] as f64;
                let y_vaf = y_vafs[i] as f64;
                if x_vaf.is_nan() || y_vaf.is_nan() {
                    continue;
                }
                let y_allele_freq = AlleleFreq(y_vaf);
                let x_allele_freq = AlleleFreq(x_vaf);

                plot_data.push(Scattervaf {
                    sample: y.to_string(),
                    normal_vaf: *x_allele_freq,
                    tumor_vaf: *y_allele_freq,
                });
            }
        }
    }

    if plot_data.is_empty() {
        return Err(errors::Error::NoRecordsFound.into());
    }

    let print_plot =
        |data: serde_json::Value, ylabel: serde_json::Value, blueprint: &str| -> Result<()> {
            let mut blueprint = serde_json::from_str(blueprint)?;
            if let Value::Object(ref mut blueprint) = blueprint {
                blueprint["data"][0]["values"] = data;
                //blueprint["axes"][0]["title"] = xlabel;
                blueprint["axes"][1]["title"] = ylabel;
                // print to STDOUT
                println!("{}", serde_json::to_string_pretty(blueprint)?);
                Ok(())
            } else {
                unreachable!();
            }
        };

    print_plot(
        json!(plot_data),
        //serde_json::value::Value::String(tumor.to_string()),
        serde_json::value::Value::String(sample_x.to_string()),
        include_str!("../../templates/plots/vaf_scatter_contour.json"),
    )
}
