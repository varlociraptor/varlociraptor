use anyhow::{Context, Result};
use bio_types::genome::Locus;
use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::record::{Buffer, Numeric};
use rust_htslib::bcf::{Format, Header, Read, Reader, Record, Writer};
use std::path::PathBuf;

use crate::utils::aux_info::{self, AuxInfo, AuxInfoCollector};
use crate::variants::model::VariantPrecision;
use crate::variants::types::breakends::Breakend;

fn create_breakend_from_record(
    record: Record,
    header: &HeaderView,
    aux_info_collector: &AuxInfoCollector,
) -> Option<Breakend> {
    let contig =
        String::from_utf8(header.rid2name(record.rid().unwrap()).unwrap().to_vec()).unwrap();

    let mut mateid = None;
    if let Ok(Some(values)) = record.info(b"MATEID").string() {
        // Nimm den ersten Wert (wenn vorhanden) und kopiere ihn in einen Vec<u8>
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
        aux_info_collector.collect(&record).unwrap(),
    )
    .unwrap();
    dbg!(&breakend);
    breakend
}

pub fn find_candidates(breakends_bcf: PathBuf, outbcf: Option<PathBuf>) -> Result<()> {
    let bcf_reader = Reader::from_path(breakends_bcf.clone()).expect("Error opening input file.");
    let mut header = Header::new();
    let aux_info_collector = AuxInfoCollector::new(&[], &bcf_reader).unwrap();
    let headerview = bcf_reader.header();
    // TODO: I create it 2 time because of ownership issues
    let mut bcf_reader = Reader::from_path(breakends_bcf).expect("Error opening input file.");

    // TODO: Find contig from old file
    let header_contig_line = format!(r#"##contig=<ID={}>"#, "J02459");
    header.push_record(header_contig_line.as_bytes());

    let mut bcf_writer: Writer;
    match outbcf {
        Some(path) => {
            bcf_writer = Writer::from_path(path, &header, true, Format::Bcf)
                .with_context(|| "Error opening BCF writer".to_string())?;
        }
        None => {
            bcf_writer = Writer::from_stdout(&header, true, Format::Bcf)
                .with_context(|| "Error opening BCF writer".to_string())?;
        }
    }

    // Wir gehen direkt durch alle Datensätze und verändern sie on-the-fly
    for record_result in bcf_reader.records() {
        let record = record_result.expect("Error reading record");
        let mut cnv_record = bcf_writer.empty_record();

        cnv_record.set_rid(record.rid());
        cnv_record.set_pos(record.pos());
        cnv_record.set_qual(f32::missing());

        let old_allele = record.alleles()[0].to_vec();
        cnv_record.set_alleles(&[old_allele.as_slice(), b"<CNV>"])?;

        let breakend = create_breakend_from_record(record, headerview, &aux_info_collector);
        // TODO: Filter on breakends to use only csv positions
        bcf_writer
            .write(&cnv_record)
            .with_context(|| "Failed to write modified BCF record".to_string())?;
    }

    Ok(())
}
