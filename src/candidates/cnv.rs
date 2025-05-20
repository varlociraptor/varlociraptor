use anyhow::{Context, Result};
use bio_types::genome::{AbstractLocus, Locus};
use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::{Format, Header, Read, Reader, Record, Writer};
use std::collections::BTreeMap;
use std::path::PathBuf;

use crate::utils::aux_info::AuxInfoCollector;
use crate::variants::model::VariantPrecision;
use crate::variants::types::breakends::Breakend;

#[derive(Default)]
pub struct Interval {
    pub start: i32,
    pub end: i32,
    pub dir: i32,
}

impl Interval {
    pub fn new(start: i32, end: i32, dir: i32) -> Self {
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
        cnv_record.set_pos(interval.start as i64);
        cnv_record.set_qual(f32::NAN);
        cnv_record.set_alleles(&[b"<CNV>"])?;
        cnv_record.push_info_integer(b"END", &[interval.end])?;
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
        b"##INFO=<ID=END,Number=1,Type=Integer,Description=\"Ending position of breakend\">",
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
            let start = b.locus().pos() as usize;
            let end = join.locus().pos() as usize;

            if start < end {
                *deltas.entry(start).or_insert(0) += 1;
                *deltas.entry(end).or_insert(0) -= 1;
            } else if start > end {
                *deltas.entry(end).or_insert(0) -= 1;
                *deltas.entry(start).or_insert(0) += 1;
            }
        }
    }
    // Compute intervals from borders
    let mut intervals_deltas = Vec::new();
    let mut last_pos = 0;
    let mut current_cov = 0;

    for (&pos, &delta) in deltas.iter() {
        if last_pos != 0 && last_pos < pos && current_cov != 0 {
            intervals_deltas.push(Interval::new(last_pos as i32, pos as i32, current_cov));
        }
        current_cov += delta;
        last_pos = pos;
    }
    intervals_deltas
}

/// Reads an input BCF file, processes each record, and writes out modified Cnv records
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
