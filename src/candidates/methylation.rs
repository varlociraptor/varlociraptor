use anyhow::{Context, Ok, Result};
use bio::io::fasta;
use itertools::Itertools;
use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{Format, Header, Writer};
use std::collections::HashSet;
use std::fs::File;
use std::path::PathBuf;

/// Find all methylation candidates "CG" in a FASTA File
///
/// # Arguments
///
/// * `infasta` - path to FASTA with genome
/// * `outbcf` - path to BCF with found methylation candidates (None for stdout)
// TODO: add implementation for other methylation types (CHH, ..., given via a pattern arg)
pub fn find_candidates(infasta: PathBuf, outbcf: Option<PathBuf>) -> Result<()> {
    // Open FASTA File
    let reader =
        fasta::Reader::from_file(infasta).with_context(|| format!("error reading FASTA file"))?;
    let mut data: Vec<(String, i64)> = vec![];

    // Collect all chromosomes and positions of candidates
    // TODO: consider using an IndexedReader to write the BCF header first and then write the methylation candidate records on the fly while reading the FASTA.
    for result in reader.records() {
        let fasta_record = result.with_context(|| format!("error parsing FASTA record"))?;
        let sequence = String::from_utf8_lossy(fasta_record.seq()).to_string();
        let candidates: Vec<_> = sequence.match_indices("CG").map(|(idx, _)| idx).collect();
        // For every candidate collect the information
        for position in candidates {
            let contig = fasta_record.id().to_owned();
            let pos = position as i64;
            data.push((contig, pos));
        }
    }
    //Write the BCF header (every contig appears once)
    let mut bcf_header = Header::new();
    for contig_id in data.clone().into_iter().map(|(contig, _)| contig).unique() {
        let header_contig_line = format!(r#"##contig=<ID={}>"#, contig_id);
        bcf_header.push_record(header_contig_line.as_bytes());
    }

    //Create a BCF writer depending on the output (to file or to stdout)
    let mut bcf_writer;
    match outbcf {
        Some(path) => {
            bcf_writer = Writer::from_path(path, &bcf_header, true, Format::Bcf)
                .with_context(|| format!("error opening BCF writer"))?;
        }
        None => {
            bcf_writer = Writer::from_stdout(&bcf_header, true, Format::Bcf)
                .with_context(|| format!("error opening BCF writer"))?;
        }
    }

    //Prepare the records
    let mut record = bcf_writer.empty_record();
    for (contig, pos) in data {
        let rid = bcf_writer
            .header()
            .name2rid(contig.as_bytes())
            .with_context(|| format!("error finding contig {contig} in header."))?;
        record.set_rid(Some(rid));
        record.set_pos(pos);
        let new_alleles: &[&[u8]] = &[b"CG", b"<METH>"];
        record.set_alleles(new_alleles);
        record.set_qual(f32::missing());

        // Write record
        bcf_writer
            .write(&record)
            .with_context(|| format!("failed to write BCF record with methylation candidate"))?;
    }
    Ok(())
}
