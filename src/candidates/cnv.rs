use anyhow::{Context, Result};
use bio_types::genome::{AbstractLocus, Locus};
use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::{Format, Header, Read, Reader, Record, Writer};
use std::path::PathBuf;

use crate::utils::aux_info::AuxInfoCollector;
use crate::variants::model::VariantPrecision;
use crate::variants::types::breakends::Breakend;

/// GRIDSS calculates quality scores according to the model outlined in the paper.
/// As GRIDSS does not yet perform multiple test correction or score recalibration, QUAL scores are vastly overestimated for all variants.
/// As a rule of thumb, variants that have QUAL >= 1000 and have assemblies from both sides of the breakpoint (AS > 0 & RAS > 0) are considered of high quality,
/// variants with QUAL >= 500 but that can only be assembled from one breakend (AS > 0 | RAS > 0) are considered of intermediate quality,
/// and variants with low QUAL score or lack any supporting assemblies are considered to be of low quality.
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
        // Locus::new(header.rid2name(record.rid())?, record.pos() as u64),
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

/// Helper function to extract the first integer value from an INFO field
fn get_first_info_int(record: &Record, tag: &[u8]) -> Option<i32> {
    record
        .info(tag)
        .integer()
        .ok()
        .flatten()
        .and_then(|v| v.first().copied())
}

/// Creates a new CNV record based on an existing BCF record
fn create_record(record: &Record, bcf_writer: &Writer, breakend: Breakend) -> Result<Record> {
    let mut cnv_record = bcf_writer.empty_record();

    cnv_record.set_rid(record.rid());
    cnv_record.set_pos(record.pos());
    cnv_record.set_qual(f32::NAN); // Sets quality to missing

    let old_allele = record.alleles()[0].to_vec();
    cnv_record.set_alleles(&[old_allele.as_slice(), b"<CNV>"])?;

    let as_val = get_first_info_int(record, b"AS");
    let ras_val = get_first_info_int(record, b"RAS");

    let gridss_qual = match (as_val, ras_val) {
        (Some(1), Some(1)) if record.qual() > 1000.0 => "HIGH",
        (Some(1), _) | (_, Some(1)) if record.qual() > 500.0 => "MEDIUM",
        _ => "LOW",
    };
    cnv_record
        .push_info_string(b"QUAL", &[gridss_qual.as_bytes()])
        .with_context(|| "Failed to push QUAL info string")?;
    let mate_locus = breakend.join().as_ref().unwrap().locus();
    let end_value = format!("{}:{}", mate_locus.contig(), mate_locus.pos());
    let end = end_value.as_str();
    cnv_record
        .push_info_string(b"END", &[end.as_bytes()])
        .with_context(|| "Failed to push END info string")?;

    Ok(cnv_record)
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

    header.push_record(
        b"##INFO=<ID=QUAL,Number=1,Type=String,Description=\"Variant classification\">",
    );
    header.push_record(
        b"##INFO=<ID=END,Number=1,Type=String,Description=\"Ending position of breakend\">",
    );

    Ok(header)
}

/// Reads an input BCF file, processes each record, and writes out modified CNV records
pub fn find_candidates(breakends_bcf: PathBuf, outbcf: Option<PathBuf>) -> Result<()> {
    let mut bcf_reader = Reader::from_path(&breakends_bcf)
        .with_context(|| format!("Error opening input file: {:?}", breakends_bcf))?;

    let header = create_header_from_existing(bcf_reader.header())?;

    let mut bcf_writer = match outbcf {
        Some(path) => Writer::from_path(path, &header, true, Format::Bcf)
            .with_context(|| "Error opening BCF writer".to_string())?,
        None => Writer::from_stdout(&header, true, Format::Bcf)
            .with_context(|| "Error opening BCF writer".to_string())?,
    };

    let header_orig = bcf_reader.header().clone();
    let aux_info = AuxInfoCollector::new(&[], &bcf_reader)
        .with_context(|| "Failed to create AuxInfoCollector")?;
    for record_result in bcf_reader.records() {
        let record = record_result.with_context(|| "Error reading record")?;
        let breakend = create_breakend_from_record(&record, &header_orig, &aux_info).unwrap();
        if breakend.is_left_to_right() {
            let cnv_record = create_record(&record, &bcf_writer, breakend)?;

            bcf_writer
                .write(&cnv_record)
                .with_context(|| "Failed to write BCF record".to_string())?;
        }
    }
    Ok(())
}
