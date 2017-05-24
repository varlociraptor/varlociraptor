use std::ops::Range;
use std::error::Error;
use std::fs;
use std::str;

use itertools::Itertools;
use rust_htslib::bcf;
use bio::io::fasta;
use bio::stats::{LogProb, PHREDProb};
use ordered_float::NotNaN;

use model;
use BCFError;
use utils;
use Event;


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
            let mut end = end[0] as u32 - 1;
            if exclusive_end {
                // this happens with DELLY
                debug!("fixing END tag");
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
                // get sequence
                let alleles = record.alleles();
                if alleles.len() > 2 {
                    return Err(Box::new(
                        BCFError::InvalidRecord("SVTYPE=INS but more than one ALT allele".to_owned())
                    ))
                }
                let ref_allele = alleles[0];
                let alt_allele = alleles[1];

                if alt_allele == b"<INS>" {
                    // don't support insertions without exact sequence
                    None
                } else {
                    let len = alt_allele.len() - ref_allele.len();

                    if is_valid_len(len as u32) {
                        Some(model::Variant::Insertion(alt_allele[ref_allele.len()..].to_owned()))
                    } else {
                        None
                    }
                }
            } else if svtype == b"DEL" {
                let svlen = match(svlen, end) {
                    (Some(svlen), _)  => svlen,
                    (None, Some(end)) => end - pos,
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
            if alt_allele[0] == b'<' {
                // skip allele if it is a special tag (such alleles have been handled above)
                None
            } else if alt_allele.len() == 1 && ref_allele.len() == 1 {
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
                    Some(model::Variant::Insertion(alt_allele[ref_allele.len()..].to_owned()))
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


/// Collect distribution of posterior probabilities from a VCF file that has been written by
/// libprosic.
pub fn collect_prob_dist<E: Event>(
    calls: &bcf::Reader,
    event: &E,
    vartype: &model::VariantType) -> Result<Vec<NotNaN<f64>>, Box<Error>> {
    let mut record = bcf::Record::new();
    let mut prob_dist = Vec::new();
    loop {
        if let Err(e) = calls.read(&mut record) {
            if e.is_eof() {
                break;
            } else {
                return Err(Box::new(e));
            }
        }

        let variants = try!(utils::collect_variants(&mut record, false, false, None, false));
        let tag = event.tag_name("PROB");
        let event_probs = try!(record.info(tag.as_bytes()).float());
        if let Some(event_probs) = event_probs {
            // tag present
            for (variant, event_prob) in variants.into_iter().zip(event_probs.into_iter()) {
                if let Some(variant) = variant {
                    if !variant.is_type(vartype) || event_prob.is_nan() {
                        continue;
                    }
                    let event_prob = LogProb::from(PHREDProb(*event_prob as f64));
                    prob_dist.push(try!(NotNaN::new(*event_prob)));
                }
            }
        }
    }
    prob_dist.sort();
    Ok(prob_dist)
}
