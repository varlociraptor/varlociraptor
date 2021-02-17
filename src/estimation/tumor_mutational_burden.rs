use std::collections::BTreeMap;
use std::str;
use std::str::FromStr;

use anyhow::Result;
use bio::stats::{LogProb, PHREDProb};
use itertools::Itertools;
use itertools_num::linspace;
use rust_htslib::bcf::{self, Read};
use serde_json::{json, Value};

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
struct TMB {
    min_vaf: f64,
    tmb: f64,
}

#[derive(Debug, Clone, Serialize)]
struct TMBStrat {
    min_vaf: f64,
    tmb: f64,
    vartype: Vartype,
}

#[derive(Debug, Clone, Serialize)]
struct TMBBin {
    vaf: f64,
    tmb: f64,
    vartype: Vartype,
}

struct Record {
    prob: LogProb,
    vartype: Vartype,
}

#[derive(
    Display,
    Debug,
    Clone,
    Copy,
    Serialize,
    Deserialize,
    EnumString,
    EnumIter,
    IntoStaticStr,
    EnumVariantNames,
)]
#[strum(serialize_all = "kebab_case")]
pub enum PlotMode {
    Hist,
    Curve,
}

/// Estimate tumor mutational burden based on Varlociraptor calls from STDIN and print result to STDOUT.
pub(crate) fn estimate(
    somatic_tumor_events: &[String],
    tumor_name: &str,
    coding_genome_size: u64,
    mode: PlotMode,
) -> Result<()> {
    let mut bcf = bcf::Reader::from_stdin()?;
    let header = bcf.header().to_owned();

    let tumor_id = bcf
        .header()
        .sample_id(tumor_name.as_bytes())
        .unwrap_or_else(|| panic!("Sample {} not found", tumor_name)); // TODO throw a proper error

    let mut tmb = BTreeMap::new();
    'records: loop {
        let mut rec = bcf.empty_record();
        match bcf.read(&mut rec) {
            None => break,
            Some(res) => res?,
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
        let vartypes = vartypes(&rec);

        // push into TMB function
        for i in 0..alt_allele_count {
            let vaf = AlleleFreq(vafs[i] as f64);
            let entry = tmb.entry(vaf).or_insert_with(Vec::new);
            entry.push(Record {
                prob: allele_probs[i],
                vartype: vartypes[i],
            });
        }
    }

    if tmb.is_empty() {
        return Err(errors::Error::NoRecordsFound.into());
    }

    let calc_tmb = |probs: &[LogProb]| -> f64 {
        let count = LogProb::ln_sum_exp(probs).exp();
        // Expected number of variants with VAF>=min_vaf per megabase.
        (count / coding_genome_size as f64) * 1000000.0
    };

    let print_plot =
        |data: serde_json::Value, blueprint: &str, cutpoint_tmb: f64, max_tmb: f64| -> Result<()> {
            let mut blueprint = serde_json::from_str(blueprint)?;
            if let Value::Object(ref mut blueprint) = blueprint {
                blueprint["data"]["values"] = data;
                blueprint["vconcat"][0]["encoding"]["y"]["scale"]["domain"] =
                    json!([cutpoint_tmb, max_tmb]);
                blueprint["vconcat"][1]["encoding"]["y"]["scale"]["domain"] =
                    json!([0.0, cutpoint_tmb]);
                // print to STDOUT
                println!("{}", serde_json::to_string_pretty(blueprint)?);
                Ok(())
            } else {
                unreachable!();
            }
        };

    let min_vafs = linspace(0.0, 1.0, 100).map(AlleleFreq);

    match mode {
        PlotMode::Hist => {
            let mut plot_data = Vec::new();
            // perform binning for histogram
            let mut max_tmbs = Vec::new();
            let mut cutpoint_tmbs = Vec::new();
            for (i, center_vaf) in linspace(0.05, 0.95, 19).enumerate() {
                let groups = tmb
                    .range(AlleleFreq(center_vaf - 0.05)..AlleleFreq(center_vaf + 0.05))
                    .map(|(_, records)| records)
                    .flatten()
                    .map(|record| (record.vartype, record.prob))
                    .into_group_map();
                for (vartype, probs) in groups {
                    let tmb = calc_tmb(&probs);
                    if i == 0 {
                        max_tmbs.push(tmb);
                    }
                    // cutpoint beyond 15%
                    if i == 2 {
                        cutpoint_tmbs.push(tmb);
                    }
                    plot_data.push(TMBBin {
                        vaf: center_vaf,
                        tmb,
                        vartype,
                    });
                }
            }

            let max_tmb: f64 = max_tmbs.iter().sum();
            let cutpoint_tmb: f64 = cutpoint_tmbs.iter().sum();

            print_plot(
                json!(plot_data),
                include_str!("../../templates/plots/vaf_hist.json"),
                cutpoint_tmb,
                max_tmb,
            )
        }
        PlotMode::Curve => {
            let mut plot_data = Vec::new();
            let mut max_tmbs = Vec::new();
            let mut cutpoint_tmbs = Vec::new();
            // calculate TMB function (expected number of somatic variants per minimum allele frequency)
            for (i, min_vaf) in min_vafs.enumerate() {
                let groups = tmb
                    .range(min_vaf..)
                    .map(|(_, records)| records)
                    .flatten()
                    .map(|record| (record.vartype, record.prob))
                    .into_group_map();

                for (vartype, probs) in groups {
                    let tmb = calc_tmb(&probs);
                    if i == 0 {
                        max_tmbs.push(tmb);
                    }
                    if i == 10 {
                        cutpoint_tmbs.push(tmb);
                    }
                    plot_data.push(TMBStrat {
                        min_vaf: *min_vaf,
                        tmb,
                        vartype,
                    });
                }
            }

            let max_tmb: f64 = max_tmbs.iter().sum();
            let cutpoint_tmb: f64 = cutpoint_tmbs.iter().sum();

            print_plot(
                json!(plot_data),
                include_str!("../../templates/plots/vaf_curve_strat.json"),
                cutpoint_tmb,
                max_tmb,
            )
        }
    }
}

#[derive(
    Display,
    Debug,
    Clone,
    Copy,
    Serialize,
    Deserialize,
    EnumString,
    EnumIter,
    IntoStaticStr,
    PartialEq,
    Hash,
    Eq,
)]
pub(crate) enum Vartype {
    DEL,
    INS,
    INV,
    DUP,
    BND,
    MNV,
    Complex,
    #[strum(serialize = "A>C")]
    #[serde(rename = "A>C")]
    AC,
    #[strum(serialize = "A>G")]
    #[serde(rename = "A>G")]
    AG,
    #[strum(serialize = "A>T")]
    #[serde(rename = "A>T")]
    AT,
    #[strum(serialize = "C>A")]
    #[serde(rename = "C>A")]
    CA,
    #[strum(serialize = "C>G")]
    #[serde(rename = "C>G")]
    CG,
    #[strum(serialize = "C>T")]
    #[serde(rename = "C>T")]
    CT,
    #[strum(serialize = "G>A")]
    #[serde(rename = "G>A")]
    GA,
    #[strum(serialize = "G>C")]
    #[serde(rename = "G>C")]
    GC,
    #[strum(serialize = "G>T")]
    #[serde(rename = "G>T")]
    GT,
    #[strum(serialize = "T>A")]
    #[serde(rename = "T>A")]
    TA,
    #[strum(serialize = "T>C")]
    #[serde(rename = "T>C")]
    TC,
    #[strum(serialize = "T>G")]
    #[serde(rename = "T>G")]
    TG,
}

pub(crate) fn vartypes(record: &bcf::Record) -> Vec<Vartype> {
    let ref_allele = record.alleles()[0];
    record.alleles()[1..]
        .iter()
        .map(|alt_allele| {
            if alt_allele == b"<DEL>" {
                Vartype::DEL
            } else if alt_allele == b"<INV>" {
                Vartype::INV
            } else if alt_allele == b"<DUP>" {
                Vartype::DUP
            } else if alt_allele == b"<BND>" {
                Vartype::BND
            } else if ref_allele.len() == 1 && alt_allele.len() == 1 {
                Vartype::from_str(&format!(
                    "{}>{}",
                    str::from_utf8(ref_allele).unwrap(),
                    str::from_utf8(alt_allele).unwrap()
                ))
                .unwrap()
            } else if ref_allele.len() > 1 && alt_allele.len() == 1 {
                Vartype::DEL
            } else if ref_allele.len() == 1 && alt_allele.len() > 1 {
                Vartype::INS
            } else if ref_allele.len() == alt_allele.len() && ref_allele.len() > 1 {
                Vartype::MNV
            } else {
                Vartype::Complex
            }
        })
        .collect_vec()
}
