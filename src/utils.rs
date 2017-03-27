use std::ops::Range;
use std::error::Error;
use std::fs;
use std::str;

use itertools::Itertools;
use rust_htslib::bcf;
use bio::io::fasta;

use model;
use BCFError;


/// Collect variants from a given ´bcf::Record`.
pub fn collect_variants(
    record: &mut bcf::Record,
    omit_snvs: bool,
    omit_indels: bool,
    indel_len_range: Option<Range<u32>>,
    exclusive_end: bool
) -> Result<Vec<Option<model::Variant>>, Box<Error>> {
    let pos = record.pos();
    let svlen = match record.info(b"SVLEN").integer() {
        Ok(Some(svlen)) => Some(svlen[0].abs() as u32),
        _ => None
    };
    let end = match record.info(b"END").integer() {
        Ok(Some(end)) => {
            let mut end = end[0] as u32;
            if exclusive_end {
                // this happens with DELLY
                end -= 1;
            }
            Some(end)
        },
        _ => None
    };
    let inslen = match record.info(b"INSLEN").integer() {
        Ok(Some(inslen)) => Some(inslen[0] as u32),
        _ => None
    };
    // TODO avoid cloning svtype
    let svtype = match record.info(b"SVTYPE").string() {
        Ok(Some(svtype)) => Some(svtype[0].to_owned()),
        _ => None
    };

    // check if len is within the given range
    let is_valid_len = |svlen| {
        if let Some(ref len_range) = indel_len_range {
            // TODO replace with Range::contains once stabilized
            if svlen < len_range.start || svlen >= len_range.end {
                return false;
            }
        }
        true
    };

    let variants = if let Some(svtype) = svtype {
        vec![
            if omit_indels {
                None
            } else if svtype == b"INS" {
                let svlen = match (svlen, inslen) {
                    (Some(svlen), _)     => svlen,
                    (None, Some(inslen)) => inslen,
                    _ => {
                        return Err(Box::new(BCFError::MissingTag("SVLEN or INSLEN".to_owned())));
                    }
                };
                if is_valid_len(svlen) {
                    Some(model::Variant::Insertion(svlen))
                } else {
                    None
                }
            } else if svtype == b"DEL" {
                let svlen = match(svlen, end) {
                    (Some(svlen), _)  => svlen,
                    (None, Some(end)) => end - 1 - pos,
                    _ => {
                        return Err(Box::new(BCFError::MissingTag("SVLEN or END".to_owned())));
                    }
                };
                if is_valid_len(svlen) {
                    Some(model::Variant::Deletion(svlen))
                } else {
                    None
                }
            } else {
                None
            }
        ]
    } else {
        let alleles = record.alleles();
        let ref_allele = alleles[0];

        alleles.iter().skip(1).map(|alt_allele| {
            if alt_allele.len() == 1 && ref_allele.len() == 1 {
                // SNV
                if omit_snvs {
                    None
                } else {
                    Some(model::Variant::SNV(alt_allele[0]))
                }
            } else if alt_allele.len() == ref_allele.len() {
                // neither indel nor SNV
                None
            } else {
                let indel_len = (alt_allele.len() as i32 - ref_allele.len() as i32).abs() as u32;

                if omit_indels {
                    None
                } else if !is_valid_len(indel_len) {
                    None
                } else if alt_allele.len() < ref_allele.len() {
                    Some(model::Variant::Deletion((ref_allele.len() - alt_allele.len()) as u32))
                } else {
                    Some(model::Variant::Insertion((alt_allele.len() - ref_allele.len()) as u32))
                }
            }
        }).collect_vec()
    };

    Ok(variants)
}


/// A lazy buffer for reference sequences.
pub struct ReferenceBuffer {
    reader: fasta::IndexedReader<fs::File>,
    chrom: Option<Vec<u8>>,
    sequence: Vec<u8>
}


impl ReferenceBuffer {
    pub fn new(fasta: fasta::IndexedReader<fs::File>) -> Self {
        ReferenceBuffer {
            reader: fasta,
            chrom: None,
            sequence: Vec::new()
        }
    }

    /// Load given chromosome and return it as a slice. This is O(1) if chromosome was loaded before.
    pub fn seq(&mut self, chrom: &[u8]) -> Result<&[u8], Box<Error>> {
        if let Some(ref last_chrom) = self.chrom {
            if last_chrom == &chrom {
                return Ok(&self.sequence);
            }
        }

        try!(self.reader.read_all(try!(str::from_utf8(chrom)), &mut self.sequence));
        self.chrom = Some(chrom.to_owned());

        Ok(&self.sequence)
    }
}
