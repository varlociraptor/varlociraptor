use anyhow::{Context, Result};
use bio_types::genome::{AbstractLocus, Locus};
use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::{Format, Header, Read, Reader, Record, Writer};
use std::collections::BTreeMap;
use std::path::PathBuf;

use crate::utils::aux_info::AuxInfoCollector;
use crate::variants::model::VariantPrecision;
use crate::variants::types::breakends::Breakend;

#[derive(Debug)]
pub struct Interval {
    pub start: Locus,
    pub end: Locus,
    pub dir: i32,
}

impl Interval {
    pub fn new(start: Locus, end: Locus, dir: i32) -> Self {
        Self { start, end, dir }
    }
}

fn write_records(intervals: Vec<Interval>, header: &Header, outbcf: Option<PathBuf>) -> Result<()> {
    let mut bcf_writer = match outbcf {
        Some(path) => Writer::from_path(path, header, true, Format::Bcf)
            .with_context(|| "Error opening BCF writer".to_string())?,
        None => Writer::from_stdout(header, true, Format::Bcf)
            .with_context(|| "Error opening BCF writer".to_string())?,
    };

    for interval in intervals {
        let mut cnv_type = b"NAN";
        if interval.dir < 0 {
            cnv_type = b"INS";
        }
        if interval.dir > 0 {
            cnv_type = b"DEL";
        }
        let mut cnv_record = bcf_writer.empty_record();
        cnv_record.set_pos(interval.start.pos() as i64);
        cnv_record.set_qual(f32::NAN);
        cnv_record.set_alleles(&[b"<CNV>"])?;
        // TODO: Get the chromosome dynamically
        let end_value = format!("{}:{}", interval.end.contig(), interval.end.pos());
        let end = end_value.as_str();
        cnv_record
            .push_info_string(b"ENDING", &[end.as_bytes()])
            .with_context(|| "Failed to push END info string")?;

        cnv_record.push_info_string(b"CNVTYPE", &[cnv_type])?;

        bcf_writer
            .write(&cnv_record)
            .with_context(|| "Failed to write BCF record".to_string())?;
    }
    Ok(())
}

/// Builds a new BCF header based on the input header
fn create_header_from_existing(existing: &HeaderView) -> Result<Header> {
    let mut header = Header::new();

    for rid in 0..existing.contig_count() {
        let name =
            std::str::from_utf8(existing.rid2name(rid)?).context("Invalid UTF-8 in contig name")?;
        let header_contig_line = format!(r#"##contig=<ID={}>"#, name);
        header.push_record(header_contig_line.as_bytes());
    }

    header.push_record(b"##INFO=<ID=CNVTYPE,Number=1,Type=String,Description=\"Type of CNV\">");
    header.push_record(
        b"##INFO=<ID=ENDING,Number=1,Type=String,Description=\"Ending position of breakend\">",
    );
    Ok(header)
}

fn create_breakend_from_record(
    record: &Record,
    header: &HeaderView,
    aux_info_collector: &AuxInfoCollector,
) -> Option<Breakend> {
    let contig =
        String::from_utf8(header.rid2name(record.rid().unwrap()).unwrap().to_vec()).unwrap();

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
        aux_info_collector.collect(record).unwrap(),
    )
    .unwrap();
    breakend
}

fn breakends_to_intervals(breakends: Vec<Breakend>) -> Vec<Interval> {
    // For every Cnv interval set the start of the interval to +1 and the end of the interval to -1
    let mut deltas = BTreeMap::new();
    for b in breakends {
        if let Some(join) = b.join() {
            let start = b.locus().to_owned();
            let end = join.locus().to_owned();
            dbg!(&start, &end);
            *deltas.entry(start).or_insert(0) += 1;
            *deltas.entry(end).or_insert(0) -= 1;
        }
    }
    // Compute intervals from borders
    let mut intervals_deltas = Vec::new();
    let mut last_pos: Option<Locus> = None;
    let mut current_cov = 0;

    for (locus, &delta) in deltas.iter() {
        if let Some(prev_locus) = last_pos {
            if prev_locus.pos() < locus.pos() && current_cov != 0 {
                intervals_deltas.push(Interval::new(prev_locus, locus.to_owned(), current_cov));
            }
        }
        current_cov += delta;
        last_pos = Some(locus.to_owned());
    }
    intervals_deltas
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
        let breakend = match create_breakend_from_record(&record, &header_orig, &aux_info) {
            Some(b) if b.is_left_to_right() && b.join().is_some() => b,
            _ => continue,
        };
        breakends.push(breakend);
    }
    let intervals: Vec<Interval> = breakends_to_intervals(breakends);
    // Write bcf file
    write_records(intervals, &header, outbcf)?;

    Ok(())
}
