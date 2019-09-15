use std::collections::BTreeMap;
use std::str;
use std::error::Error;

use bio::stats::{LogProb,PHREDProb};
use rust_htslib::bcf::{self, Read};
use itertools::Itertools;
use serde_json::{Value,json};

use crate::model::AlleleFreq;

/// Consider only variants in coding regions.
/// We rely on the ANN field for this.
fn is_valid_variant(rec: &mut bcf::Record) -> Result<bool, Box<Error>> {
    
    for ann in rec.info(b"ANN").string()?.expect("ANN field not found. Annotate VCF with e.g. snpEff.") {
        for entry in ann.split(|c| *c == b'|') {
            if entry == b"protein_coding" {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

#[derive(Debug, Clone, Serialize)]
struct TMB {
    min_vaf: f64,
    tmb: f64,
}

/// Estimate tumor mutational burden and print result to STDOUT.
pub fn estimate(bcf: &mut bcf::Reader, somatic_tumor_events: Vec<&[u8]>, tumor_name: &[u8], coding_genome_size: u64) -> Result<(), Box<Error>> {
    let tumor_id = bcf.header().sample_id(tumor_name).expect(&format!("Sample {} not found", str::from_utf8(tumor_name).unwrap()));

    let mut tmb = BTreeMap::new();
    for res in bcf.records() {
        let mut rec = res?;

        if !is_valid_variant(&mut rec)? {
            continue;
        }

        let alt_allele_count = (rec.allele_count() - 1) as usize;

        // collect allele probabilities for given events
        let allele_probs = vec![LogProb::ln_zero(); alt_allele_count];
        for e in &somatic_tumor_events {
            let probs = rec.info(e).float()?.expect(&format!("Event {} not found in VCF/BCF file.", str::from_utf8(e).unwrap()));
            for i in 0..alt_allele_count {
                allele_probs[i].ln_add_exp(LogProb::from(PHREDProb(probs[i] as f64)));
            }
        }

        // obtain VAF estimates
        let vafs = rec.format(b"AF").float()?[tumor_id];

        // push into TMB function
        for i in 0..alt_allele_count {
            let vaf = AlleleFreq(vafs[i] as f64);
            let vaf_probs = tmb.entry(vaf).or_insert_with(|| Vec::new());
            vaf_probs.push(allele_probs[i]);
        }
    }

    let mut plot_data = Vec::new();
    // calculate TMB function (expected number of somatic variants per minimum allele frequency)
    for min_vaf in tmb.keys() {
        let probs = tmb.range(min_vaf..).map(|(_, probs)| probs).flatten().cloned().collect_vec();
        let count = *LogProb::ln_sum_exp(&probs);
        // Expected number of variants with VAF>=min_vaf per megabase.
        let tmb = count / coding_genome_size as f64 * 1000000.0;
        plot_data.push(TMB { min_vaf: **min_vaf, tmb });
    }

    let mut blueprint = serde_json::from_str(include_str!("../../templates/plots/tmb.json"))?;
    if let Value::Object(ref mut blueprint) = blueprint {
        blueprint["data"] = json!(plot_data);
        // print to STDOUT
        println!("{}", serde_json::to_string_pretty(blueprint)?);
    } else {
        unreachable!();
    }

    Ok(())
}