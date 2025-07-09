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
    cnv_intervals: Vec<Interval>,
    header: &Header,
    outbcf: Option<PathBuf>,
) -> Result<()> {
    let mut bcf_writer = match outbcf {
        Some(path) => Writer::from_path(path, header, true, Format::Bcf)
            .with_context(|| "Error opening BCF writer".to_string())?,
        None => Writer::from_stdout(header, true, Format::Bcf)
            .with_context(|| "Error opening BCF writer".to_string())?,
    };

    for interval in cnv_intervals {
        let mut cnv_record = bcf_writer.empty_record();

        cnv_record.set_pos(interval.range().start as i64);
        cnv_record.set_qual(f32::missing());
        // TODO: What is the REF allele?
        cnv_record.set_alleles(&[b"N", b"<CNV>"])?;
        cnv_record
            .push_info_integer(b"END", &[interval.range().end as i32])
            .with_context(|| "Failed to push END info string")?;
        cnv_record
            .push_info_string(b"SVTYPE", &[b"CNV"])
            .with_context(|| "Failed to push SVTYPE info string")?;

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

fn breakends_to_intervals(breakends: Vec<Breakend>) -> Vec<Interval> {
    let mut all_intervals = BTreeMap::new();
    for b in breakends {
        if let Some(join) = b.join() {
            let (breakpoint_start, breakpoint_end) = (b.locus(), join.locus());

            let (interval_start, interval_end) = if breakpoint_start.pos() < breakpoint_end.pos() {
                (breakpoint_start.to_owned(), breakpoint_end.to_owned())
            } else {
                (breakpoint_end.to_owned(), breakpoint_start.to_owned())
            };

            *all_intervals.entry(interval_start).or_insert(0) += 1;
            *all_intervals.entry(interval_end).or_insert(0) -= 1;
        }
    }
    // Compute intervals from borders
    let mut cnv_intervals = Vec::new();
    let mut last_pos: Option<Locus> = None;
    let mut current_cov = 0;

    for (locus, &delta) in all_intervals.iter() {
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
                cnv_intervals.push(interval);
            }
        }
        current_cov += delta;
        last_pos = Some(locus.to_owned());
    }
    cnv_intervals
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
        breakends.push(breakend);
    }
    let intervals: Vec<Interval> = breakends_to_intervals(breakends);
    // Write bcf file
    write_cnv_records(intervals, &header, outbcf)?;

    Ok(())
}
