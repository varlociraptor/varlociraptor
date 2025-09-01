use anyhow::{Context, Result};
use bio_types::genome::{AbstractInterval, AbstractLocus, Interval, Locus};
use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::{record::Numeric, Format, Header, Read, Reader, Record, Writer};
use std::collections::BTreeMap;
use std::ops::Range;
use std::path::PathBuf;

use crate::candidates;
use crate::utils::aux_info::AuxInfoCollector;
use crate::variants::model::VariantPrecision;
use crate::variants::types::breakends::Breakend;

fn write_cnv_records(
    cnv_candidates: Vec<BreakendInterval>,
    header: &Header,
    outbcf: Option<PathBuf>,
) -> Result<()> {
    let mut bcf_writer = match outbcf {
        Some(path) => Writer::from_path(path, header, true, Format::Bcf)
            .with_context(|| "Error opening BCF writer".to_string())?,
        None => Writer::from_stdout(header, true, Format::Bcf)
            .with_context(|| "Error opening BCF writer".to_string())?,
    };
    for candidate in cnv_candidates {
        let mut cnv_record = bcf_writer.empty_record();
        cnv_record.set_pos(candidate.interval.range().start as i64);
        cnv_record.set_qual(f32::missing());
        // TODO: What is the REF allele?
        cnv_record.set_alleles(&[b"N", b"<CNV>"])?;
        cnv_record
            .push_info_integer(b"END", &[candidate.interval.range().end as i32])
            .with_context(|| "Failed to push END info string")?;
        cnv_record
            .push_info_string(b"SVTYPE", &[b"CNV"])
            .with_context(|| "Failed to push SVTYPE info string")?;
        for (k, v) in &candidate.info_probs {
            if k == "AF" {
                cnv_record
                    .push_format_float(k.as_bytes(), &[(10.0_f64).powf(-*v / 10.0) as f32])
                    .with_context(|| format!("Failed to push {} format string", k))?;
            } else {
                cnv_record
                    .push_info_float(k.as_bytes(), &[*v as f32])
                    .with_context(|| format!("Failed to push {} info string", k))?;
            }
        }

        bcf_writer
            .write(&cnv_record)
            .with_context(|| "Failed to write BCF record".to_string())?;
    }
    Ok(())
}

// /// Builds a new BCF header based on the input header
// fn create_header_from_existing(old_header: &HeaderView) -> Result<Header> {
//     let mut header = Header::new();

//     for rid in 0..old_header.contig_count() {
//         let name = std::str::from_utf8(old_header.rid2name(rid)?)
//             .context("Invalid UTF-8 in contig name")?;
//         let header_contig_line = format!(r#"##contig=<ID={}>"#, name);
//         header.push_record(header_contig_line.as_bytes());
//     }
//     header.push_record(
//         b"##INFO=<ID=END,Number=1,Type=Integer,Description=\"Ending position of breakend\">",
//     );
//     header.push_record(
//         b"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
//     );
//     header.push_record(b"##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">");
//     Ok(header)
// }

fn create_header_from_existing(old_header: &HeaderView) -> anyhow::Result<Header> {
    // kompletten Header aus der View klonen
    let mut header = Header::from_template(old_header);
    header.push_record(
        b"##INFO=<ID=END,Number=1,Type=Integer,Description=\"Ending position of breakend\">",
    );
    Ok(header)
}

fn parse_breakend(
    record: &Record,
    header: &HeaderView,
    aux_info_collector: &AuxInfoCollector,
) -> Result<Option<Breakend>> {
    let contig = String::from_utf8(
        header
            .rid2name(record.rid().context("Missing rid")?)
            .context("Invalid RID name")?
            .to_vec(),
    )
    .context("Could not create String from rid")?;

    let mut mateid = None;
    if let Ok(Some(values)) = record.info(b"MATEID").string() {
        mateid = values.first().map(|v| v.to_vec());
    }
    let breakend = Breakend::new(
        Locus::new(contig, record.pos() as u64),
        record.alleles()[0],
        record.alleles()[1],
        &record.id(),
        mateid,
        VariantPrecision::Precise,
        aux_info_collector
            .collect(record)
            .context("Auxiliary Info could not be collected from record")?,
    )
    .context("Breakend could not be created")?;
    Ok(breakend)
}

fn breakends_to_intervals(candidates: Vec<CnvCandidate>) -> Vec<BreakendInterval> {
    #[derive(Clone, Debug, Default)]
    struct BorderProbs {
        prob_interval_start: f64,
        prob_interval_end: f64,
        count_prob_zero: i32,
    }

    impl BorderProbs {
        fn update_start(&mut self, value: f64) {
            if value != f64::NEG_INFINITY {
                self.prob_interval_start += value;
            } else {
                self.count_prob_zero += 1;
            }
        }

        fn update_end(&mut self, value: f64) {
            if value != f64::NEG_INFINITY {
                self.prob_interval_end += value;
            } else {
                self.count_prob_zero += 1;
            }
        }
    }

    #[derive(Clone, Debug, Default)]
    struct BorderInfo {
        delta: i32,
        info_probs_borders: BTreeMap<String, BorderProbs>,
    }

    #[derive(Clone, Debug, Default)]
    struct ProbTracking {
        current_prob_product: f64,
        count_prob_zero: i32,
    }

    let mut all_intervals: BTreeMap<Locus, BorderInfo> = BTreeMap::new();

    // 1) Sammle alle Borders
    for c in candidates {
        if let Some(join) = c.breakend.join() {
            let (bp_start, bp_end) = (c.breakend.locus(), join.locus());
            let (interval_start, interval_end) = if bp_start.pos() < bp_end.pos() {
                (bp_start.to_owned(), bp_end.to_owned())
            } else {
                (bp_end.to_owned(), bp_start.to_owned())
            };

            // Start-Border
            let entry_start = all_intervals.entry(interval_start).or_default();
            entry_start.delta += 1;
            for (k, &v) in &c.info_probs {
                entry_start
                    .info_probs_borders
                    .entry(k.clone())
                    .or_default()
                    .update_start(v);
            }

            // End-Border
            let entry_end = all_intervals.entry(interval_end).or_default();
            entry_end.delta += 1;
            for (k, &v) in &c.info_probs {
                entry_end
                    .info_probs_borders
                    .entry(k.clone())
                    .or_default()
                    .update_end(v);
            }
        }
    }

    // 2) Erzeuge Intervalle
    let mut cnv_intervals = Vec::new();
    let mut last_pos: Option<Locus> = None;
    let mut current_cov = 0;
    let mut interval_probs: BTreeMap<String, ProbTracking> = BTreeMap::new();
    let mut current_info_probs: BTreeMap<String, f64> = BTreeMap::new();
    for (locus, border_info) in all_intervals.iter() {
        // Update zero counts
        for (k, v) in &border_info.info_probs_borders {
            if v.count_prob_zero > 0 {
                interval_probs.entry(k.clone()).or_default().count_prob_zero += v.count_prob_zero;
            }
        }

        // Erzeuge Intervalle
        if let Some(prev_locus) = last_pos {
            if prev_locus.contig() == locus.contig()
                && prev_locus.pos() < locus.pos()
                && current_cov != 0
            {
                let interval = Interval::new(
                    prev_locus.contig().to_owned(),
                    Range {
                        start: prev_locus.pos() + 1,
                        end: locus.pos(),
                    },
                );
                cnv_intervals.push(BreakendInterval {
                    interval,
                    info_probs: current_info_probs.clone(),
                });
            }
        }
        // Update Wahrscheinlichkeiten
        for (k, v) in &border_info.info_probs_borders {
            let old_prob = current_info_probs.get(k).copied().unwrap_or(0.0);
            let new_prob = if v.count_prob_zero > 0 {
                0.0
            } else {
                old_prob + (v.prob_interval_start - v.prob_interval_end)
            };
            current_info_probs.insert(k.clone(), new_prob);
        }

        current_cov += border_info.delta;
        last_pos = Some(locus.to_owned());
    }

    cnv_intervals
}

#[derive(Builder, Getters, Debug)]
// Computed interval out of the breakends
pub struct BreakendInterval {
    #[getset(get = "pub")]
    interval: Interval,
    #[getset(get = "pub")]
    info_probs: BTreeMap<String, f64>,
}

#[derive(Builder, Getters, Debug)]
// Representative for every breakend in the call file
pub struct CnvCandidate {
    #[getset(get = "pub")]
    breakend: Breakend,
    #[getset(get = "pub")]
    info_probs: BTreeMap<String, f64>,
}

/// Reads an input BCF file, processes each record, and writes out modified CNV records
pub fn find_candidates(breakends_bcf: PathBuf, outbcf: Option<PathBuf>) -> Result<()> {
    let mut bcf_reader = Reader::from_path(&breakends_bcf)
        .with_context(|| format!("Error opening input file: {:?}", breakends_bcf))?;

    // Create header for interval bcf
    let header = create_header_from_existing(bcf_reader.header())?;
    let header_orig = bcf_reader.header().clone();

    let mut breakends = Vec::new();

    // TODO: Compute dynamicall (problem: Cant extract INFO information from Header(View)
    let prob_keys = vec![
        "PROB_HIGH".to_string(),
        "PROB_LOW".to_string(),
        "PROB_ABSENT".to_string(),
        "PROB_ARTIFACT".to_string(),
    ];
    let aux_info = AuxInfoCollector::new(&[], &bcf_reader)
        .with_context(|| "Failed to create AuxInfoCollector")?;

    // Exract breakends + information out of every record in this for loop
    for record_result in bcf_reader.records() {
        let record = record_result.with_context(|| "Error reading record")?;

        let breakend = match parse_breakend(&record, &header_orig, &aux_info)? {
            Some(b) if b.is_left_to_right() && b.join().is_some() => b,
            _ => continue,
        };

        let af = {
            let af_values = record.format(b"AF").float()?;
            if af_values.is_empty() || af_values[0].is_empty() {
                f64::NAN
            } else {
                -10.0 * (af_values[0][0] as f64).log10()
            }
        };

        let mut probs = BTreeMap::new();
        for key in &prob_keys {
            if let Ok(Some(value)) = record.info(key.as_bytes()).float() {
                if let Some(&val) = value.first() {
                    // Probs in the PROB tags are PHRED encoded, convert to linear scale

                    probs.insert(key.clone(), val as f64);
                }
            }
        }
        probs.insert("AF".to_string(), af);
        // Create new breakend candidate out of record
        let candidate = CnvCandidateBuilder::default()
            .breakend(breakend)
            .info_probs(probs)
            .build()?;

        breakends.push(candidate);
    }
    dbg!(&breakends);
    // Compute intervals out of all breakends
    let intervals: Vec<BreakendInterval> = breakends_to_intervals(breakends);
    dbg!(&intervals);
    write_cnv_records(intervals, &header, outbcf)?;

    Ok(())
}
