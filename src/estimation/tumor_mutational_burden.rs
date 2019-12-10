use std::collections::BTreeMap;
use std::error::Error;
use std::str;

use bio::stats::{LogProb, PHREDProb};
use itertools::Itertools;
use rust_htslib::bcf::{self, Read};
use serde_json::{json, Value};

use crate::errors;
use crate::model::AlleleFreq;
use crate::{Event, SimpleEvent};

/// Consider only variants in coding regions.
/// We rely on the ANN field for this.
fn is_valid_variant(rec: &mut bcf::Record) -> Result<bool, Box<dyn Error>> {
    for ann in rec
        .info(b"ANN")
        .string()?
        .expect("ANN field not found. Annotate VCF with e.g. snpEff.")
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
struct TMB {
    min_vaf: f64,
    tmb: f64,
}

/// Estimate tumor mutational burden based on Varlociraptor calls from STDIN and print result to STDOUT.
pub fn estimate(
    somatic_tumor_events: &[String],
    tumor_name: &str,
    coding_genome_size: u64,
) -> Result<(), Box<dyn Error>> {
    let mut bcf = bcf::Reader::from_stdin()?;
    let header = bcf.header().to_owned();

    let tumor_id = bcf
        .header()
        .sample_id(tumor_name.as_bytes())
        .expect(&format!("Sample {} not found", tumor_name));

    let mut tmb = BTreeMap::new();
    'records: loop {
        let mut rec = bcf.empty_record();
        if !bcf.read(&mut rec)? {
            break;
        }

        let contig = str::from_utf8(header.rid2name(rec.rid().unwrap()).unwrap())?;
        let vcfpos = rec.pos() + 1;
        // obtain VAF estimates (do it here already to work around a segfault in htslib)
        let vafs = rec.format(b"AF").float()?[tumor_id].to_owned();

        if !is_valid_variant(&mut rec)? {
            info!(
                "Skipping variant {}:{} because it is not coding.",
                contig, vcfpos
            );
            continue;
        }

        let alt_allele_count = (rec.allele_count() - 1) as usize;

        // collect allele probabilities for given events
        let mut allele_probs = vec![LogProb::ln_zero(); alt_allele_count];
        for e in somatic_tumor_events {
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

        // push into TMB function
        for i in 0..alt_allele_count {
            let vaf = AlleleFreq(vafs[i] as f64);
            let vaf_probs = tmb.entry(vaf).or_insert_with(|| Vec::new());
            vaf_probs.push(allele_probs[i]);
        }
    }

    if tmb.is_empty() {
        return Err(errors::Error::NoRecordsFound)?;
    }

    let mut plot_data = Vec::new();
    // calculate TMB function (expected number of somatic variants per minimum allele frequency)
    for min_vaf in tmb.keys() {
        let probs = tmb
            .range(min_vaf..)
            .map(|(_, probs)| probs)
            .flatten()
            .cloned()
            .collect_vec();
        let count = LogProb::ln_sum_exp(&probs).exp();
        // Expected number of variants with VAF>=min_vaf per megabase.
        let tmb = (count / coding_genome_size as f64) * 1000000.0;
        plot_data.push(TMB {
            min_vaf: **min_vaf,
            tmb,
        });
    }

    let mut blueprint = serde_json::from_str(include_str!("../../templates/plots/tmb.json"))?;
    if let Value::Object(ref mut blueprint) = blueprint {
        blueprint["data"]["values"] = json!(plot_data);
        // print to STDOUT
        println!("{}", serde_json::to_string_pretty(blueprint)?);
    } else {
        unreachable!();
    }

    Ok(())
}
