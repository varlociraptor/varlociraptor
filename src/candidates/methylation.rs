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
/// * `outbcf` - path to VCF with found methylation candidates (None for stdout)
// TODO: add implementation for other methylation types (CHH, ..., given via a pattern arg)
pub fn find_candidates(infasta: PathBuf, outvcf: Option<PathBuf>) -> Result<()> {
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
    //Write the VCF header (every contig appears once)
    let mut vcf_header = Header::new();
    for contig_id in data.clone().into_iter().map(|(contig, _)| contig).unique() {
        let header_contig_line = format!(r#"##contig=<ID={}>"#, contig_id);
        vcf_header.push_record(header_contig_line.as_bytes());
    }

    //Create a VCF writer depending on the output (to file or to stdout)
    let mut vcf_writer;
    match outvcf {
        Some(path) => {
            vcf_writer = Writer::from_path(path, &vcf_header, true, Format::Vcf)
                .with_context(|| format!("error opening BCF writer"))?;
        }
        None => {
            vcf_writer = Writer::from_stdout(&vcf_header, true, Format::Vcf)
                .with_context(|| format!("error opening BCF writer"))?;
        }
    }

    //Prepare the records
    let mut record = vcf_writer.empty_record();
    for (contig, pos) in data {
        let rid = vcf_writer
            .header()
            .name2rid(contig.as_bytes())
            .with_context(|| format!("error finding contig {contig} in header."))?;
        record.set_rid(Some(rid));
        record.set_pos(pos);
        let new_alleles: &[&[u8]] = &[b"CG", b"<METH>"];
        record.set_alleles(new_alleles);
        record.set_qual(f32::missing());

        // Write record
        vcf_writer
            .write(&record)
            .with_context(|| format!("failed to write BCF record with methylation candidate"))?;
    }
    Ok(())
}
