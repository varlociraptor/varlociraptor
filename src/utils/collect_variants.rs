use std::convert::TryFrom;

use anyhow::Result;
use bio_types::genome::AbstractLocus;
use itertools::Itertools;
use rust_htslib::{bcf, bcf::record::Numeric};

use crate::errors;
use crate::utils::SimpleCounter;
use crate::variants::model::{self, HaplotypeIdentifier, VariantPrecision};

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

#[derive(Debug, Getters, Clone)]
#[getset(get = "pub(crate)")]
pub(crate) struct VariantInfo {
    pub(crate) variant: model::Variant,
    pub(crate) haplotype: Option<HaplotypeIdentifier>,
}

/// Collect variants from a given ´bcf::Record`.
pub fn collect_variants(
    record: &mut bcf::Record,
    skip_imprecise: bool,
    skips: Option<&mut SimpleCounter<SkipReason>>,
) -> Result<Vec<VariantInfo>> {
    let imprecise = record.info(b"IMPRECISE").flag().ok().unwrap_or(false);

    let skip_incr = |reason| {
        if let Some(skips) = skips {
            skips.incr(reason);
        }
    };

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

    let is_valid_insertion_alleles = |ref_allele: &[u8], alt_allele: &[u8]| {
        alt_allele == b"<INS>"
            || (ref_allele.len() < alt_allele.len()
                && ref_allele == &alt_allele[..ref_allele.len()]
                && ref_allele.len() == 1)
    };

    let is_valid_deletion_alleles = |ref_allele: &[u8], alt_allele: &[u8]| {
        alt_allele == b"<DEL>"
            || (ref_allele.len() > alt_allele.len()
                && &ref_allele[..alt_allele.len()] == alt_allele
                && alt_allele.len() == 1)
    };

    let haplotype = HaplotypeIdentifier::from(record)?;

    let mut variants = Vec::new();

    let mut push_variant = |variant| {
        variants.push(VariantInfo {
            variant,
            haplotype: haplotype.clone(),
        })
    };

    if skip_imprecise && imprecise && svtype.as_ref().map_or(false, |svtype| svtype != b"BND") {
        skip_incr(SkipReason::Imprecise);
        return Ok(Vec::with_capacity(1));
    }

    if let Some(svtype) = svtype {
        if svtype == b"INV" {
            let alleles = record.alleles();
            if alleles.len() != 2 {
                skip_incr(SkipReason::InversionInvalidAlt);
            } else if let Some(end) = end {
                let len = end + 1 - pos; // end is inclusive, pos as well.
                push_variant(model::Variant::Inversion(len));
            } else {
                skip_incr(SkipReason::InversionMissingEndTag);
            }
        } else if svtype == b"DUP" {
            let alleles = record.alleles();
            if alleles.len() != 2 {
                skip_incr(SkipReason::DuplicationInvalidAlt);
            } else if let Some(end) = end {
                let len = end + 1 - pos; // end is inclusive, pos as well.
                push_variant(model::Variant::Duplication(len));
            } else {
                skip_incr(SkipReason::DuplicationMissingEndTag)
            }
        } else if svtype == b"BND" {
            let alleles = record.alleles();
            if haplotype.is_some() {
                for spec in &alleles[1..] {
                    let rec: &bcf::Record = &*record;
                    push_variant(model::Variant::Breakend {
                        ref_allele: alleles[0].to_owned(),
                        spec: spec.to_vec(),
                        precision: VariantPrecision::try_from(rec)?,
                    })
                }
            } else {
                skip_incr(SkipReason::BreakendNoEvent);
            }
        } else if svtype == b"INS" {
            // get sequence
            let alleles = record.alleles();
            if alleles.len() > 2 {
                return Err(errors::invalid_bcf_record(
                    record.contig(),
                    record.pos(),
                    "SVTYPE=INS but more than one ALT allele",
                )
                .into());
            }
            let ref_allele = alleles[0];
            let alt_allele = alleles[1];

            if alt_allele != b"<INS>" {
                // don't support insertions without exact sequence
                if is_valid_insertion_alleles(ref_allele, alt_allele) {
                    push_variant(model::Variant::Insertion(
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
                return Err(errors::invalid_bcf_record(
                    record.contig(),
                    record.pos(),
                    "Absolute value of SVLEN or END - POS must be greater than zero.",
                )
                .into());
            }
            let alleles = record.alleles();
            if alleles.len() > 2 {
                return Err(errors::invalid_bcf_record(
                    record.contig(),
                    record.pos(),
                    "SVTYPE=DEL but more than one ALT allele",
                )
                .into());
            }
            let ref_allele = alleles[0];
            let alt_allele = alleles[1];

            if alt_allele == b"<DEL>" || is_valid_deletion_alleles(ref_allele, alt_allele) {
                push_variant(model::Variant::Deletion(svlen));
            }
        }
    } else {
        let alleles = record.alleles();
        let ref_allele = alleles[0];

        for (i, alt_allele) in alleles.iter().skip(1).enumerate() {
            if alt_allele == b"<*>" {
                // dummy non-ref allele, signifying potential homozygous reference site
                push_variant(model::Variant::None);
            } else if alt_allele == b"<DEL>" {
                if let Some(ref svlens) = svlens {
                    if let Some(svlen) = svlens[i] {
                        push_variant(model::Variant::Deletion(svlen));
                    }
                    // TODO fail with an error in else case
                }
            } else if alt_allele[0] == b'<' {
                // skip any other special alleles
            } else if alt_allele.len() == 1 && ref_allele.len() == 1 {
                // SNV
                push_variant(model::Variant::Snv(alt_allele[0]));
            } else if alt_allele.len() == ref_allele.len() {
                // MNV
                push_variant(model::Variant::Mnv(alt_allele.to_vec()));
            } else if is_valid_deletion_alleles(ref_allele, alt_allele) {
                push_variant(model::Variant::Deletion(
                    (ref_allele.len() - alt_allele.len()) as u64,
                ));
            } else if is_valid_insertion_alleles(ref_allele, alt_allele) {
                push_variant(model::Variant::Insertion(
                    alt_allele[ref_allele.len()..].to_owned(),
                ));
            } else {
                // arbitrary replacement
                push_variant(model::Variant::Replacement {
                    ref_allele: ref_allele.to_owned(),
                    alt_allele: alt_allele.to_vec(),
                });
            }
        }
    }

    Ok(variants)
}
