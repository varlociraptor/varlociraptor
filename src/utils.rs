use std::ops::Range;
use std::error::Error;

use itertools::Itertools;
use rust_htslib::bcf;

use model;
use BCFError;


/// Collect variants from a given ´bcf::Record`.
pub fn collect_variants(
    record: &mut bcf::Record,
    omit_snvs: bool,
    omit_indels: bool,
    indel_len_range: Option<Range<u32>>
) -> Result<Vec<Option<model::Variant>>, Box<Error>> {
    let pos = record.pos();
    let svlen = match record.info(b"SVLEN").integer() {
        Ok(Some(svlen)) => Some(svlen[0].abs() as u32),
        _ => None
    };
    let end = match record.info(b"END").integer() {
        Ok(Some(end)) => Some(end[0] as u32),
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
            if svlen < len_range.start && svlen >= len_range.end {
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
                } else if is_valid_len(indel_len) {
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
