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
    // TODO avoid cloning svtype
    let svtypes = record.info(b"SVTYPE").string().map(|values| {
        values.map(|values| values.into_iter().map(|svtype| svtype.to_owned()).collect_vec())
    });

    let variants = if let Ok(Some(svtypes)) = svtypes {
        // obtain svlen if it is present or raise error
        let svlens = record.info(b"SVLEN").integer().map(|values| values.map(|values| values.to_owned()));
        let pos = record.pos() as i32;
        let svlens = match svlens {
            Ok(Some(svlens)) => svlens,
            _ => {
                match record.info(b"END").integer() {
                    Ok(Some(ends)) => ends.into_iter().map(|end| end - 1 - pos).collect_vec(),
                    Err(e) => return Err(Box::new(e)),
                    Ok(None) => return Err(Box::new(BCFError::MissingTag("SVLEN".to_owned())))
                }
            }
        };
        svtypes.iter().zip(svlens).map(|(svtype, svlen)| {
            let svlen = svlen.abs() as u32;
            if let Some(ref len_range) = indel_len_range {
                // TODO replace with Range::contains once stabilized
                if svlen < len_range.start && svlen >= len_range.end {
                    return None;
                }
            }

            if omit_indels {
                None
            } else if svtype == b"INS" {
                Some(model::Variant::Insertion(svlen))
            } else if svtype == b"DEL" {
                Some(model::Variant::Deletion(svlen))
            } else {
                None
            }
        }).collect_vec()
    } else {
        let alleles = record.alleles();
        let ref_allele = alleles[0];

        alleles.iter().skip(1).map(|alt_allele| {
            if alt_allele.len() == 1 && ref_allele.len() == 1 {
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
                if let Some(ref len_range) = indel_len_range {
                    // TODO replace with Range::contains once stabilized
                    if indel_len < len_range.start && indel_len >= len_range.end {
                        return None;
                    }
                }

                if omit_indels {
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
