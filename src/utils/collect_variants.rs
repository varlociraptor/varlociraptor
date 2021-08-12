use anyhow::Result;
use itertools::Itertools;
use rust_htslib::{bcf, bcf::record::Numeric};

use crate::errors;
use crate::utils::SimpleCounter;
use crate::variants::model;

#[derive(
    Hash, PartialEq, Eq, EnumString, EnumIter, IntoStaticStr, EnumVariantNames, Display, Debug,
)]
pub(crate) enum SkipReason {
    #[strum(serialize = "imprecise variants (will be supported in a future release)")]
    Imprecise,
    #[strum(serialize = "inversions with missing END tag")]
    InversionMissingEndTag,
    #[strum(serialize = "duplications with missing END tag")]
    DuplicationMissingEndTag,
    #[strum(serialize = "inversion with more than a single <INV> allele")]
    InversionInvalidAlt,
    #[strum(serialize = "duplication with more than a single <DUP> allele")]
    DuplicationInvalidAlt,
    #[strum(serialize = "breakend without EVENT tag (will be supported in a future release)")]
    BreakendNoEvent,
}

/// Collect variants from a given ´bcf::Record`.
pub(crate) fn collect_variants(
    record: &mut bcf::Record,
    skip_imprecise: bool,
    skips: Option<&mut SimpleCounter<SkipReason>>,
) -> Result<Vec<model::Variant>> {
    let imprecise = record.info(b"IMPRECISE").flag().ok().unwrap_or(false);

    let skip_incr = |reason| {
        if let Some(skips) = skips {
            skips.incr(reason);
        }
    };

    if skip_imprecise && imprecise {
        skip_incr(SkipReason::Imprecise);
        return Ok(Vec::with_capacity(1));
    }

    let pos = record.pos() as u64;
    let svlens = match record.info(b"SVLEN").integer() {
        Ok(Some(svlens)) => Some(
            svlens
                .iter()
                .map(|l| {
                    if !l.is_missing() {
                        Some(l.abs() as u64)
                    } else {
                        None
                    }
                })
                .collect_vec(),
        ),
        _ => None,
    };
    let end = match record.info(b"END").integer() {
        Ok(Some(end)) => {
            let end = end[0] as u64 - 1;
            Some(end)
        }
        _ => None,
    };
    // TODO avoid cloning svtype
    let svtype = match record.info(b"SVTYPE").string() {
        Ok(Some(svtype)) => Some(svtype[0].to_owned()),
        _ => None,
    };

    let event = match record.info(b"EVENT").string() {
        Ok(Some(event)) => Some(event[0].to_owned()),
        _ => None,
    };

    let is_valid_insertion_alleles = |ref_allele: &[u8], alt_allele: &[u8]| {
        alt_allele == b"<INS>"
            || (ref_allele.len() < alt_allele.len()
                && ref_allele == &alt_allele[..ref_allele.len()])
    };

    let is_valid_deletion_alleles = |ref_allele: &[u8], alt_allele: &[u8]| {
        alt_allele == b"<DEL>"
            || (ref_allele.len() > alt_allele.len()
                && &ref_allele[..alt_allele.len()] == alt_allele)
    };

    let mut variants = Vec::new();

    if let Some(svtype) = svtype {
        if svtype == b"INV" {
            let alleles = record.alleles();
            if alleles.len() != 2 {
                skip_incr(SkipReason::InversionInvalidAlt);
            } else if let Some(end) = end {
                let len = end + 1 - pos; // end is inclusive, pos as well.
                variants.push(model::Variant::Inversion(len));
            } else {
                skip_incr(SkipReason::InversionMissingEndTag);
            }
        } else if svtype == b"DUP" {
            let alleles = record.alleles();
            if alleles.len() != 2 {
                skip_incr(SkipReason::DuplicationInvalidAlt);
            } else if let Some(end) = end {
                let len = end + 1 - pos; // end is inclusive, pos as well.
                variants.push(model::Variant::Duplication(len));
            } else {
                skip_incr(SkipReason::DuplicationMissingEndTag)
            }
        } else if svtype == b"BND" {
            let alleles = record.alleles();
            if let Some(ref event) = event {
                for spec in &alleles[1..] {
                    variants.push(model::Variant::Breakend {
                        event: event.clone(),
                        ref_allele: alleles[0].to_owned(),
                        spec: spec.to_vec(),
                    })
                }
            } else {
                skip_incr(SkipReason::BreakendNoEvent);
            }
        } else if svtype == b"INS" {
            // get sequence
            let alleles = record.alleles();
            if alleles.len() > 2 {
                return Err(errors::Error::InvalidBCFRecord {
                    msg: "SVTYPE=INS but more than one ALT allele".to_owned(),
                }
                .into());
            }
            let ref_allele = alleles[0];
            let alt_allele = alleles[1];

            if alt_allele != b"<INS>" {
                // don't support insertions without exact sequence
                if is_valid_insertion_alleles(ref_allele, alt_allele) {
                    variants.push(model::Variant::Insertion(
                        alt_allele[ref_allele.len()..].to_owned(),
                    ));
                }
            }
        } else if svtype == b"DEL" {
            let svlen = match (svlens, end) {
                (Some(ref svlens), _) if svlens[0].is_some() => svlens[0].unwrap(),
                (None, Some(end)) => end - (pos + 1), // pos is pointing to the allele before the DEL
                _ => {
                    return Err(errors::Error::MissingBCFTag {
                        name: "SVLEN or END".to_owned(),
                    }
                    .into());
                }
            };
            if svlen == 0 {
                return Err(errors::Error::InvalidBCFRecord {
                    msg: "Absolute value of SVLEN or END - POS must be greater than zero."
                        .to_owned(),
                }
                .into());
            }
            let alleles = record.alleles();
            if alleles.len() > 2 {
                return Err(errors::Error::InvalidBCFRecord {
                    msg: "SVTYPE=DEL but more than one ALT allele".to_owned(),
                }
                .into());
            }
            let ref_allele = alleles[0];
            let alt_allele = alleles[1];

            if alt_allele == b"<DEL>" || is_valid_deletion_alleles(ref_allele, alt_allele) {
                variants.push(model::Variant::Deletion(svlen));
            }
        }
    } else {
        let alleles = record.alleles();
        let ref_allele = alleles[0];

        for (i, alt_allele) in alleles.iter().skip(1).enumerate() {
            if alt_allele == b"<*>" {
                // dummy non-ref allele, signifying potential homozygous reference site
                variants.push(model::Variant::None);
            } else if alt_allele == b"<DEL>" {
                if let Some(ref svlens) = svlens {
                    if let Some(svlen) = svlens[i] {
                        variants.push(model::Variant::Deletion(svlen));
                    }
                    // TODO fail with an error in else case
                }
            } else if alt_allele[0] == b'<' {
                // skip any other special alleles
            } else if alt_allele.len() == 1 && ref_allele.len() == 1 {
                // SNV
                variants.push(model::Variant::Snv(alt_allele[0]));
            } else if alt_allele.len() == ref_allele.len() {
                // MNV
                variants.push(model::Variant::Mnv(alt_allele.to_vec()));
            } else {
                // TODO fix position if variant is like this: cttt -> ct
                if is_valid_deletion_alleles(ref_allele, alt_allele) {
                    variants.push(model::Variant::Deletion(
                        (ref_allele.len() - alt_allele.len()) as u64,
                    ));
                } else if is_valid_insertion_alleles(ref_allele, alt_allele) {
                    variants.push(model::Variant::Insertion(
                        alt_allele[ref_allele.len()..].to_owned(),
                    ));
                } else {
                    // arbitrary replacement
                    variants.push(model::Variant::Replacement {
                        ref_allele: ref_allele.to_owned(),
                        alt_allele: alt_allele.to_vec(),
                    });
                }
            }
        }
    }

    Ok(variants)
}
