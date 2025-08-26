use anyhow::{Context, Result};
use bio_types::genome::{AbstractInterval, AbstractLocus, Interval, Locus};
use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::{record::Numeric, Format, Header, Read, Reader, Record, Writer};
use std::collections::BTreeMap;
use std::ops::Range;
use std::path::PathBuf;

use crate::utils::aux_info::AuxInfoCollector;
use crate::variants::model::VariantPrecision;
use crate::variants::types::breakends::Breakend;

fn write_cnv_records(
    cnv_candidates: Vec<BreakendIterval>,
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
        cnv_record
            .push_info_float(b"AF", &[candidate.af as f32])
            .with_context(|| "Failed to push AF info string")?;

        bcf_writer
            .write(&cnv_record)
            .with_context(|| "Failed to write BCF record".to_string())?;
    }
    Ok(())
}

/// Builds a new BCF header based on the input header
fn create_header_from_existing(old_header: &HeaderView) -> Result<Header> {
    let mut header = Header::new();

    for rid in 0..old_header.contig_count() {
        let name = std::str::from_utf8(old_header.rid2name(rid)?)
            .context("Invalid UTF-8 in contig name")?;
        let header_contig_line = format!(r#"##contig=<ID={}>"#, name);
        header.push_record(header_contig_line.as_bytes());
    }
    header.push_record(
        b"##INFO=<ID=END,Number=1,Type=Integer,Description=\"Ending position of breakend\">",
    );
    header.push_record(
        b"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
    );
    header.push_record(b"##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">");
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

fn breakends_to_intervals(candidates: Vec<CnvCandidate>) -> Vec<BreakendIterval> {
    #[derive(Default, Clone, Debug)]
    struct BorderInfo {
        delta: i32,
        af_new_interval: f64,
        af_old_interval: f64,
        zero_af_count: i32,
        event_probs_prod: BTreeMap<String, f64>,
    }

    impl BorderInfo {
        fn new() -> Self {
            Self {
                delta: 0,
                af_new_interval: 1.0,
                af_old_interval: 1.0,
                zero_af_count: 0,
                event_probs_prod: BTreeMap::new(),
            }
        }
    }

    let mut all_intervals: BTreeMap<Locus, BorderInfo> = BTreeMap::new();

    for c in candidates {
        if let Some(join) = c.breakend.join() {
            let (breakpoint_start, breakpoint_end) = (c.breakend.locus(), join.locus());

            let (interval_start, interval_end) = if breakpoint_start.pos() < breakpoint_end.pos() {
                (breakpoint_start.to_owned(), breakpoint_end.to_owned())
            } else {
                (breakpoint_end.to_owned(), breakpoint_start.to_owned())
            };

            let entry_start = all_intervals.entry(interval_start).or_default();
            entry_start.delta += 1;
            if c.af > 0.0 {
                entry_start.af_new_interval *= c.af;
            } else {
                entry_start.zero_af_count += 1;
            }
            for (k, v) in &c.event_probs {
                *entry_start.event_probs_prod.entry(k.clone()).or_insert(1.0) *= v;
            }

            let entry_end = all_intervals.entry(interval_end).or_default();
            entry_end.delta -= 1;
            if c.af > 0.0 {
                entry_end.af_old_interval *= c.af;
            } else {
                entry_end.zero_af_count -= 1;
            }
            for (k, v) in &c.event_probs {
                *entry_end.event_probs_prod.entry(k.clone()).or_insert(1.0) *= v;
            }
        }
    }
    // Compute intervals from borders
    let mut cnv_intervals = Vec::new();
    let mut last_pos: Option<Locus> = None;
    let mut current_cov = 0;
    let mut af_product = 1.0;
    let mut zero_af = 0;
    for (locus, border_info) in all_intervals.iter() {
        zero_af += border_info.zero_af_count;
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
                if zero_af > 0 || border_info.zero_af_count < 0 {
                    // at least one zero AF in the product, so AF is zero
                    cnv_intervals.push(BreakendIterval {
                        interval,
                        af: 0.0,
                        event_probs: border_info.event_probs_prod.clone(),
                    });
                } else {
                    // all AFs > 0, so product is valid
                    cnv_intervals.push(BreakendIterval {
                        interval,
                        af: af_product,
                        event_probs: border_info.event_probs_prod.clone(),
                    });
                }
            }
        }
        af_product *= border_info.af_new_interval;

        current_cov += border_info.delta;
        af_product /= border_info.af_old_interval;

        last_pos = Some(locus.to_owned());
    }
    cnv_intervals
}

#[derive(Builder, Getters)]
pub struct BreakendIterval {
    #[getset(get = "pub")]
    interval: Interval,
    #[getset(get = "pub")]
    af: f64,
    #[getset(get = "pub")]
    event_probs: BTreeMap<String, f64>,
}

#[derive(Builder, Getters, Debug)]
pub struct CnvCandidate {
    #[getset(get = "pub")]
    breakend: Breakend,
    #[getset(get = "pub")]
    af: f64,
    #[getset(get = "pub")]
    event_probs: BTreeMap<String, f64>,
}

/// Reads an input BCF file, processes each record, and writes out modified CNV records
pub fn find_candidates(breakends_bcf: PathBuf, outbcf: Option<PathBuf>) -> Result<()> {
    let mut bcf_reader = Reader::from_path(&breakends_bcf)
        .with_context(|| format!("Error opening input file: {:?}", breakends_bcf))?;

    let header = create_header_from_existing(bcf_reader.header())?;

    let mut breakends = Vec::new();

    let header_orig = bcf_reader.header().clone();
    let aux_info = AuxInfoCollector::new(&[], &bcf_reader)
        .with_context(|| "Failed to create AuxInfoCollector")?;
    // Get all breakends from the input BCF file
    for record_result in bcf_reader.records() {
        let record = record_result.with_context(|| "Error reading record")?;
        let breakend = match parse_breakend(&record, &header_orig, &aux_info)? {
            Some(b) if b.is_left_to_right() && b.join().is_some() => b,
            _ => continue,
        };
        let af: f64 = {
            let af_values = record.format(b"AF").float()?; // BufferBacked<Vec<&[f32]>, _>
            if af_values.is_empty() || af_values[0].is_empty() {
                f64::NAN
            } else {
                af_values[0][0] as f64
            }
        };
        let candidate = CnvCandidateBuilder::default()
            .breakend(breakend)
            .af(af)
            .event_probs(BTreeMap::new())
            .build()?;

        breakends.push(candidate);
    }
    let intervals: Vec<BreakendIterval> = breakends_to_intervals(breakends);
    // Write bcf file
    write_cnv_records(intervals, &header, outbcf)?;

    Ok(())
}
