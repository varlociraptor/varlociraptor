use aho_corasick::AhoCorasick;
use anyhow::{Context, Ok, Result};
use bio::io::fasta::Reader;
use itertools::Itertools;
use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{Format, Header, Writer};
use std::path::PathBuf;

/// Find all methylation candidates in a FASTA File
///
/// # Arguments
///
/// * `infasta` - path to FASTA with genome
/// * `outbcf` - path to BCF with found methylation candidates (None for stdout)
// TODO: add implementation for other methylation types (CHH, ..., given via a pattern arg)
pub fn find_candidates(
    infasta: PathBuf,
    motifs: Vec<String>,
    outbcf: Option<PathBuf>,
) -> Result<()> {
    // Open FASTA File
    let reader =
        Reader::from_file(infasta).with_context(|| "error reading FASTA file".to_string())?;
    let mut data = Vec::new();
    let ac_motifs = AhoCorasick::new(&motifs)?;
    // Collect all chromosomes and positions of candidates
    // TODO: consider using an IndexedReader to write the BCF header first and then write the methylation candidate records on the fly while reading the FASTA.
    for result in reader.records() {
        let fasta_record = result.with_context(|| "error parsing FASTA record".to_string())?;
        let sequence = String::from_utf8_lossy(fasta_record.seq()).to_string();
        let contig = fasta_record.id().to_owned();

        for mat in ac_motifs.find_iter(&sequence) {
            let pos = mat.start() as i64;
            let motif = motifs[mat.pattern()].to_string();
            data.push((contig.clone(), pos, motif));
        }
    }
    //Write the BCF header (every contig appears once)
    let mut bcf_header = Header::new();
    for contig_id in data
        .clone()
        .into_iter()
        .map(|(contig, _, _)| contig)
        .unique()
    {
        let header_contig_line = format!(r#"##contig=<ID={}>"#, contig_id);
        bcf_header.push_record(header_contig_line.as_bytes());
    }

    //Create a BCF writer depending on the output (to file or to stdout)
    let mut bcf_writer;
    match outbcf {
        Some(path) => {
            bcf_writer = Writer::from_path(path, &bcf_header, true, Format::Bcf)
                .with_context(|| "error opening BCF writer".to_string())?;
        }
        None => {
            bcf_writer = Writer::from_stdout(&bcf_header, true, Format::Bcf)
                .with_context(|| "error opening BCF writer".to_string())?;
        }
    }

    //Prepare the records
    let mut record = bcf_writer.empty_record();
    for (contig, pos, motif) in data {
        let rid = bcf_writer
            .header()
            .name2rid(contig.as_bytes())
            .with_context(|| format!("error finding contig {contig} in header."))?;
        record.set_rid(Some(rid));
        record.set_pos(pos);
        let new_alleles: &[&[u8]] = &[motif.as_bytes(), b"<METH>"];
        record
            .set_alleles(new_alleles)
            .with_context(|| "error setting alleles".to_string())?;
        record.set_qual(f32::missing());

        // Write record
        bcf_writer
            .write(&record)
            .with_context(|| "failed to write BCF record with methylation candidate".to_string())?;
    }
    Ok(())
}
