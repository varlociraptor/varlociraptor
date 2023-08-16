use bio::io::fasta;
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

pub fn find_candidates(infasta: PathBuf, outvcf: Option<PathBuf>) {
    // Open FASTA File
    let fasta_file: File = File::open(infasta).expect("Unable to open");
    let reader = fasta::Reader::new(fasta_file);
    let mut data: Vec<(String, i64)> = vec![];

    // Collect all chromosomes and positions of candidates
    for result in reader.records() {
        let fasta_record = result.expect("Error during fasta record parsing");
        let sequence: String = String::from_utf8_lossy(fasta_record.seq()).to_string();
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
    let unique_strs: HashSet<_> = data.clone().into_iter().map(|(s, _)| s).collect();
    let unique_contig_list: Vec<String> = unique_strs.into_iter().collect();
    for contig_id in unique_contig_list {
        let header_contig_line: String = format!(r#"##contig=<ID={}>"#, contig_id);
        vcf_header.push_record(header_contig_line.as_bytes());
    }

    //Create a VCF writer depending on the output (to file or to stdout)
    let mut vcf_writer;
    match outvcf {
        Some(path) => {
            vcf_writer = Writer::from_path(path, &vcf_header, true, Format::Vcf).unwrap();
        }
        None => {
            vcf_writer = Writer::from_stdout(&vcf_header, true, Format::Vcf).unwrap();
        }
    }

    //Prepare the records
    let mut record = vcf_writer.empty_record();
    for rec in data {
        let rid = vcf_writer
            .header()
            .name2rid(rec.0.to_string().as_bytes())
            .unwrap();
        record.set_rid(Some(rid));
        record.set_pos(rec.1);
        let new_alleles: &[&[u8]] = &[b"CG", b"<METH>"];
        record.set_alleles(new_alleles);
        record.set_qual(f32::missing());

        // Write record
        match vcf_writer.write(&record) {
            Ok(_) => {}
            Err(e) => {
                eprintln!("Error: {}", e)
            }
        }
    }
}
