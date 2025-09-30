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
fn is_valid_variant(rec: &mut bcf::Record, header: &bcf::header::HeaderView) -> Result<bool> {
    if let Some(ann) = rec.info(b"ANN").string()? {
        for ann in ann.iter() {
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
    } else {
        warn!(
            "No ANN field found in record at {}:{}.",
            std::str::from_utf8(header.rid2name(rec.rid().unwrap())?)?,
            rec.pos() + 1
        );
        Ok(false)
    }
}

#[derive(Debug, Clone, Serialize)]
struct MBStrat {
    min_vaf: f64,
    mb: f64,
    vartype: Signature,
}

#[derive(Debug, Clone, Serialize)]
struct MBBin {
    vaf: f64,
    mb: f64,
    vartype: Signature,
}

#[derive(Debug, Clone, Serialize)]
struct MBMultiBar {
    vaf: f64,
    mb: f64,
    vartype: Signature,
    sample: String,
}

struct RecordSig {
    prob: LogProb,
    vartype: Signature,
    sample: String,
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
    VariantNames,
)]
#[strum(serialize_all = "kebab_case")]
pub enum Mode {
    Hist,
    Curve,
    Multibar,
    Table,
}

pub(crate) fn collect_estimates(
    mutational_events: &[String],
    sample_names: &[String],
    coding_genome_size: u64,
    mode: Mode,
    cutoff: f64,
) -> Result<()> {
    let mut bcf = bcf::Reader::from_stdin()?;
    let header = bcf.header().to_owned();

    let mut sample_ids = BTreeMap::new();

    for s in sample_names {
        let sample_id = bcf
            .header()
            .sample_id(s.as_bytes())
            .unwrap_or_else(|| panic!("Sample {} not found", s));
        sample_ids.insert(s, sample_id);
    }

    let mut mb = BTreeMap::new();
    'records: loop {
        let mut rec: bcf::Record = bcf.empty_record();
        match bcf.read(&mut rec) {
            None => break,
            Some(res) => res?,
        }

        let contig = str::from_utf8(header.rid2name(rec.rid().unwrap()).unwrap())?;
        let vcfpos = rec.pos() + 1;
        // obtain VAF estimates (do it here already to work around a segfault in htslib)
        let mut vafmap = BTreeMap::new();
        for (name, id) in &sample_ids {
            let vafs = rec.format(b"AF").float()?[*id].to_owned();
            vafmap.insert(name, vafs);
        }
        if !is_valid_variant(&mut rec, &header)? {
            info!(
                "Skipping variant {}:{} because it is not coding.",
                contig, vcfpos
            );
            continue;
        }

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
        let vartypes = signatures(&rec);

        // push into MB function
        for sample_name in sample_ids.keys() {
            for i in 0..alt_allele_count {
                // if all alt_alleles are NaN, the list will only contain one NaN, so check for size
                if i == vafmap.get(sample_name).unwrap().len()
                    && vafmap.get(sample_name).unwrap()[0].is_nan()
                {
                    continue;
                }
                let vaf = vafmap.get(sample_name).unwrap()[i] as f64;
                if vaf.is_nan() {
                    continue;
                }
                let allele_freq = AlleleFreq(vaf);
                let entry = mb.entry(allele_freq).or_insert_with(Vec::new);
                entry.push(RecordSig {
                    prob: allele_probs[i],
                    vartype: vartypes[i],
                    sample: sample_name.to_string(),
                });
            }
        }
    }
    if mb.is_empty() {
        return Err(errors::Error::NoRecordsFound.into());
    }

    let calc_mb = |probs: &[LogProb]| -> f64 {
        let count = LogProb::ln_sum_exp(probs).exp();
        // Expected number of variants with VAF>=min_vaf per megabase.
        (count / coding_genome_size as f64) * 1000000.0
    };

    let print_plot =
        |data: serde_json::Value, blueprint: &str, cutpoint_mb: f64, max_mb: f64| -> Result<()> {
            let mut blueprint = serde_json::from_str(blueprint)?;
            if let Value::Object(ref mut blueprint) = blueprint {
                blueprint["data"]["values"] = data;
                if relative_eq!(cutpoint_mb, max_mb) {
                    blueprint["vconcat"][0]["encoding"]["y"]["scale"]["domain"] =
                        json!([0.0, max_mb]);
                } else {
                    blueprint["vconcat"][0]["encoding"]["y"]["scale"]["domain"] =
                        json!([cutpoint_mb, max_mb]);
                    blueprint["vconcat"][1]["encoding"]["y"]["scale"]["domain"] =
                        json!([0.0, cutpoint_mb]);
                }
                // print to STDOUT
                println!("{}", serde_json::to_string_pretty(blueprint)?);
                Ok(())
            } else {
                unreachable!();
            }
        };

    let min_vafs = linspace(0.0, 1.0, 100).map(AlleleFreq);

    match mode {
        Mode::Multibar => {
            let mut plot_data = Vec::new();
            let mut max_mb = 0.0;
            // calculate mutational burden (mb) function (expected number of variants per minimum allele frequency)
            let groups = mb
                .range(AlleleFreq(cutoff)..AlleleFreq(1.0))
                .flat_map(|(_, records)| records)
                .map(|record| ((record.vartype, record.sample.clone()), record.prob))
                .into_group_map();

            for ((vartype, sample), probs) in groups {
                let mb = calc_mb(&probs);
                if mb > max_mb {
                    max_mb = mb;
                }

                plot_data.push(MBMultiBar {
                    vaf: cutoff,
                    mb,
                    vartype,
                    sample,
                });
            }

            print_plot(
                json!(plot_data),
                include_str!("../../templates/plots/vaf_multi_bar.json"),
                max_mb,
                max_mb,
            )
        }
        Mode::Hist => {
            let mut plot_data = Vec::new();
            // perform binning for histogram
            let mut max_mbs = Vec::new();
            let mut cutpoint_mbs = Vec::new();
            for (i, center_vaf) in linspace(0.05, 0.95, 19).enumerate() {
                let groups = mb
                    .range(AlleleFreq(center_vaf - 0.05)..AlleleFreq(center_vaf + 0.05))
                    .flat_map(|(_, records)| records)
                    .map(|record| (record.vartype, record.prob))
                    .into_group_map();
                for (vartype, probs) in groups {
                    let mb = calc_mb(&probs);
                    if i == 0 {
                        max_mbs.push(mb);
                    }
                    // cutpoint beyond 15%
                    if i == 2 {
                        cutpoint_mbs.push(mb);
                    }
                    plot_data.push(MBBin {
                        vaf: center_vaf,
                        mb,
                        vartype,
                    });
                }
            }

            let max_mb: f64 = max_mbs.iter().sum();
            let cutpoint_mb: f64 = cutpoint_mbs.iter().sum();

            print_plot(
                json!(plot_data),
                include_str!("../../templates/plots/vaf_hist.json"),
                cutpoint_mb,
                max_mb,
            )
        }
        Mode::Curve => {
            let mut plot_data = Vec::new();
            let mut max_mbs = Vec::new();
            let mut cutpoint_mbs = Vec::new();
            // calculate MB function (expected number of variants per minimum allele frequency)
            for (i, min_vaf) in min_vafs.enumerate() {
                let groups = mb
                    .range(min_vaf..)
                    .flat_map(|(_, records)| records)
                    .map(|record| (record.vartype, record.prob))
                    .into_group_map();

                for (vartype, probs) in groups {
                    let mb = calc_mb(&probs);
                    if i == 0 {
                        max_mbs.push(mb);
                    }
                    if i == 10 {
                        cutpoint_mbs.push(mb);
                    }
                    plot_data.push(MBStrat {
                        min_vaf: *min_vaf,
                        mb,
                        vartype,
                    });
                }
            }

            let max_mb: f64 = max_mbs.iter().sum();
            let cutpoint_mb: f64 = cutpoint_mbs.iter().sum();

            print_plot(
                json!(plot_data),
                include_str!("../../templates/plots/vaf_curve_strat.json"),
                cutpoint_mb,
                max_mb,
            )
        }
        Mode::Table => {
            let mut writer = csv::WriterBuilder::new()
                .delimiter(b'\t')
                .from_writer(std::io::stdout());
            // calculate mutational burden (mb) function (expected number of variants per minimum allele frequency)
            for min_vaf in min_vafs {
                let groups = mb
                    .range(min_vaf..)
                    .flat_map(|(_, records)| records)
                    .map(|record| (record.vartype, record.prob))
                    .into_group_map();

                for (vartype, probs) in groups {
                    let mb = calc_mb(&probs);
                    writer.serialize(MBStrat {
                        min_vaf: *min_vaf,
                        mb,
                        vartype,
                    })?;
                }
            }
            Ok(())
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
    #[strum(serialize = "DEL")]
    #[serde(rename = "DEL")]
    Del,
    #[strum(serialize = "INS")]
    #[serde(rename = "INS")]
    Ins,
    #[strum(serialize = "METH")]
    #[serde(rename = "METH")]
    Meth,
    #[strum(serialize = "INV")]
    #[serde(rename = "INV")]
    Inv,
    #[strum(serialize = "DUP")]
    #[serde(rename = "DUP")]
    Dup,
    #[strum(serialize = "BND")]
    #[serde(rename = "BND")]
    Bnd,
    #[strum(serialize = "MNV")]
    #[serde(rename = "MNV")]
    Mnv,
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
pub(crate) enum Signature {
    #[strum(serialize = "DEL")]
    #[serde(rename = "DEL")]
    Del,
    #[strum(serialize = "METH")]
    #[serde(rename = "METH")]
    Meth,
    #[strum(serialize = "INS")]
    #[serde(rename = "INS")]
    Ins,
    #[strum(serialize = "INV")]
    #[serde(rename = "INV")]
    Inv,
    #[strum(serialize = "DUP")]
    #[serde(rename = "DUP")]
    Dup,
    #[strum(serialize = "BND")]
    #[serde(rename = "BND")]
    Bnd,
    #[strum(serialize = "MNV")]
    #[serde(rename = "MNV")]
    Mnv,
    Complex,
    #[strum(serialize = "C>A", serialize = "G>T")]
    #[serde(rename = "C>A")]
    CA,
    #[strum(serialize = "C>G", serialize = "G>C")]
    #[serde(rename = "C>G")]
    CG,
    #[strum(serialize = "C>T", serialize = "G>A")]
    #[serde(rename = "C>T")]
    CT,
    #[strum(serialize = "T>A", serialize = "A>T")]
    #[serde(rename = "T>A")]
    TA,
    #[strum(serialize = "T>C", serialize = "A>G")]
    #[serde(rename = "T>C")]
    TC,
    #[strum(serialize = "T>G", serialize = "A>C")]
    #[serde(rename = "T>G")]
    TG,
}

pub(crate) fn signatures(record: &bcf::Record) -> Vec<Signature> {
    let ref_allele = record.alleles()[0];
    record.alleles()[1..]
        .iter()
        .map(|alt_allele| {
            if alt_allele == b"<DEL>" {
                Signature::Del
            } else if alt_allele == b"<INV>" {
                Signature::Inv
            } else if alt_allele == b"<DUP>" {
                Signature::Dup
            } else if alt_allele == b"<BND>" {
                Signature::Bnd
            } else if alt_allele == b"<METH>" {
                Signature::Meth
            } else if ref_allele.len() == 1 && alt_allele.len() == 1 {
                Signature::from_str(&format!(
                    "{}>{}",
                    str::from_utf8(ref_allele).unwrap(),
                    str::from_utf8(alt_allele).unwrap()
                ))
                .unwrap()
            } else if ref_allele.len() > 1 && alt_allele.len() == 1 {
                Signature::Del
            } else if ref_allele.len() == 1 && alt_allele.len() > 1 {
                Signature::Ins
            } else if ref_allele.len() == alt_allele.len() && ref_allele.len() > 1 {
                Signature::Mnv
            } else {
                Signature::Complex
            }
        })
        .collect_vec()
}

#[allow(dead_code)]
pub(crate) fn vartypes(record: &bcf::Record) -> Vec<Vartype> {
    let ref_allele = record.alleles()[0];
    record.alleles()[1..]
        .iter()
        .map(|alt_allele| {
            if alt_allele == b"<DEL>" {
                Vartype::Del
            } else if alt_allele == b"<INV>" {
                Vartype::Inv
            } else if alt_allele == b"<DUP>" {
                Vartype::Dup
            } else if alt_allele == b"<BND>" {
                Vartype::Bnd
            } else if alt_allele == b"<METH>" {
                Vartype::Meth
            } else if ref_allele.len() == 1 && alt_allele.len() == 1 {
                Vartype::from_str(&format!(
                    "{}>{}",
                    str::from_utf8(ref_allele).unwrap(),
                    str::from_utf8(alt_allele).unwrap()
                ))
                .unwrap()
            } else if ref_allele.len() > 1 && alt_allele.len() == 1 {
                Vartype::Del
            } else if ref_allele.len() == 1 && alt_allele.len() > 1 {
                Vartype::Ins
            } else if ref_allele.len() == alt_allele.len() && ref_allele.len() > 1 {
                Vartype::Mnv
            } else {
                Vartype::Complex
            }
        })
        .collect_vec()
}
