use std::collections::BTreeMap;
use std::str;
use std::str::FromStr;

use anyhow::Result;
use bio::stats::{LogProb, PHREDProb};
use itertools::Itertools;
use itertools_num::linspace;
use rust_htslib::bcf::{self, Read};
use serde_json::{json, Value};
use kodama::{Method, linkage};

use crate::errors;
use crate::variants::model::AlleleFreq;
use crate::{Event, SimpleEvent};


/// Consider only variants in coding regions.
/// We rely on the ANN field for this.
fn is_valid_variant(rec: &mut bcf::Record) -> Result<bool> {
    for ann in rec
        .info(b"ANN")
        .string()?
        .expect("ANN field not found. Annotate VCF with e.g. snpEff.")
        .iter()
    {
        let mut coding = false;
        for (i, entry) in ann.split(|c| *c == b'|').enumerate() {
            if i == 7 {
                coding = entry == b"protein_coding";
            }
            if i == 13 {
                coding &= entry != b"";
            }
        }
        if coding {
            return Ok(true);
        }
    }
    Ok(false)
}

#[derive(Debug, Clone, Serialize)]
struct SCATTERVAF {
    sample: String,
    normal_vaf: f64,
    tumor_vaf: f64,
}


pub(crate) fn vaf_scatter(
    mutational_events: &[String],
    normal: &String,
    tumor: &[String],
) -> Result<()> {
    let mut bcf = bcf::Reader::from_stdin()?;
    let header = bcf.header().to_owned();

    
    let mut plot_data = Vec::new();
    
    let normal_id = bcf
    .header()
    .sample_id(normal.as_bytes())
    .unwrap_or_else(|| panic!("Sample {} not found", normal));

    let mut tumor_ids = BTreeMap::new();

    for s in tumor {
        let tumor_id = bcf
            .header()
            .sample_id(s.as_bytes())
            .unwrap_or_else(|| panic!("Sample {} not found", s));
        tumor_ids.insert(s, tumor_id);
    }

    'records: loop {
        let mut rec = bcf.empty_record();
        match bcf.read(&mut rec) {
            None => break,
            Some(res) => res?,
        }

        let contig = str::from_utf8(header.rid2name(rec.rid().unwrap()).unwrap())?;
        let vcfpos = rec.pos() + 1;

        if !is_valid_variant(&mut rec)? {
            info!(
                "Skipping variant {}:{} because it is not coding.",
                contig, vcfpos
            );
            continue;
        }

        // obtain VAF estimates (do it here already to work around a segfault in htslib)
        let normal_vafs = rec.format(b"AF").float()?[normal_id].to_owned();



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
                continue 'records;
            }
        }

        // push into MB function

        for i in 0..alt_allele_count {
            for (t, tumor_id) in &tumor_ids {
                let tumor_vafs = rec.format(b"AF").float()?[*tumor_id].to_owned();
                // if all alt_alleles are NaN, the list will only contain one NaN, so check for size
                if (i == normal_vafs.len() && normal_vafs[0].is_nan()) || (i == tumor_vafs.len() && tumor_vafs[0].is_nan())
                {
                    continue;
                }
                let normal_vaf = normal_vafs[i] as f64;
                let tumor_vaf = tumor_vafs[i] as f64;
                if normal_vaf.is_nan() || tumor_vaf.is_nan() {
                    continue;
                }
                let tumor_allele_freq = AlleleFreq(tumor_vaf);
                let normal_allele_freq = AlleleFreq(normal_vaf);            

                plot_data.push(SCATTERVAF {
                    sample: t.to_string(),
                    normal_vaf: *normal_allele_freq,
                    tumor_vaf: *tumor_allele_freq,
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
        serde_json::value::Value::String(normal.to_string()),
        include_str!("../../templates/plots/vaf_scatter_contour.json"),
    )
}