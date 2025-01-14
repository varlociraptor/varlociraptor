use std::convert::TryFrom;

use anyhow::{bail, Result};
use bio::stats::{LogProb, Prob};
use bio_types::genome::AbstractLocus;
use itertools::Itertools;
use rust_htslib::{bcf, bcf::record::Numeric};

use crate::errors;
use crate::utils::SimpleCounter;
use crate::variants::model::{self, HaplotypeIdentifier, VariantPrecision};

#[derive(
    Hash, PartialEq, Eq, EnumString, EnumIter, IntoStaticStr, VariantNames, Display, Debug,
)]
pub enum SkipReason {
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

#[derive(Debug, Getters, CopyGetters, Clone)]
pub struct VariantInfo {
    #[getset(get = "pub(crate)")]
    pub(crate) variant: model::Variant,
    #[getset(get = "pub(crate)")]
    pub(crate) haplotype: Option<HaplotypeIdentifier>,
    #[getset(get_copy = "pub(crate)")]
    heterozygosity: Option<LogProb>,
    #[getset(get_copy = "pub(crate)")]
    somatic_effective_mutation_rate: Option<LogProb>,
}

/// Collect variants from a given ´bcf::Record`.
pub fn collect_variants(
    record: &mut bcf::Record,
    skip_imprecise: bool,
    skips: Option<&mut SimpleCounter<SkipReason>>,
    variant_heterozygosity_field: Option<&[u8]>,
    variant_somatic_effective_mutation_rate_field: Option<&[u8]>,
) -> Result<Vec<VariantInfo>> {
    // TODO ignore imprecise variants for now?
    // let nonzero_bounds = |tag| {
    //     record.info(tag).integer().ok().map_or(false, |ci| {
    //         ci.map_or(false, |ci| ci.iter().any(|bound| *bound != 0))
    //     })
    // };

    // let imprecise = record.info(b"IMPRECISE").flag().ok().unwrap_or(false)
    //     || (!record.info(b"PRECISE").flag().ok().unwrap_or(false)
    //         && (nonzero_bounds(b"CIPOS") || nonzero_bounds(b"CIEND")));

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
                        Some(l.unsigned_abs() as u64)
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

    let get_plain_prob_field = |key: Option<&[u8]>| -> Result<Vec<Option<LogProb>>> {
        if let Some(key) = key {
            let field = record.info(key);
            if let Some(values) = field.float()? {
                if values.len() != record.alleles().len() - 1 {
                    bail!(errors::Error::InvalidVariantPrior);
                }
                return Ok(values
                    .iter()
                    .map(|value| {
                        if value.is_missing() {
                            None
                        } else {
                            Some(LogProb::from(Prob((*value) as f64)))
                        }
                    })
                    .collect());
            }
        }
        Ok(vec![None; record.alleles().len() - 1])
    };

    let variant_heterozygosities = get_plain_prob_field(variant_heterozygosity_field)?;
    let variant_somatic_effective_mutation_rates =
        get_plain_prob_field(variant_somatic_effective_mutation_rate_field)?;

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

    let mut push_variant = |variant, i: usize| {
        variants.push(VariantInfo {
            variant,
            haplotype: haplotype.clone(),
            heterozygosity: variant_heterozygosities[i],
            somatic_effective_mutation_rate: variant_somatic_effective_mutation_rates[i],
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
                push_variant(model::Variant::Inversion(len), 0);
            } else {
                skip_incr(SkipReason::InversionMissingEndTag);
            }
        } else if svtype == b"DUP" {
            let alleles = record.alleles();
            if alleles.len() != 2 {
                skip_incr(SkipReason::DuplicationInvalidAlt);
            } else if let Some(end) = end {
                let len = end + 1 - pos; // end is inclusive, pos as well.
                push_variant(model::Variant::Duplication(len), 0);
            } else {
                skip_incr(SkipReason::DuplicationMissingEndTag)
            }
        } else if svtype == b"BND" {
            let alleles = record.alleles();
            if haplotype.is_some() {
                for (i, spec) in alleles[1..].iter().enumerate() {
                    let rec: &bcf::Record = &*record;

                    push_variant(
                        model::Variant::Breakend {
                            ref_allele: alleles[0].to_owned(),
                            spec: spec.to_vec(),
                            precision: VariantPrecision::try_from(rec)?,
                        },
                        i,
                    )
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
                    push_variant(
                        model::Variant::Insertion(alt_allele[ref_allele.len()..].to_owned()),
                        0,
                    );
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
                push_variant(model::Variant::Deletion(svlen), 0);
            }
        }
    } else {
        let alleles = record.alleles();
        let ref_allele = alleles[0];

        for (i, alt_allele) in alleles.iter().skip(1).enumerate() {
            if alt_allele == b"<*>" {
                // dummy non-ref allele, signifying potential homozygous reference site
                push_variant(model::Variant::None, i);
            } else if alt_allele == b"<DEL>" {
                if let Some(ref svlens) = svlens {
                    if let Some(svlen) = svlens[i] {
                        push_variant(model::Variant::Deletion(svlen), i);
                    }
                    // TODO fail with an error in else case
                }
            } else if alt_allele[0] == b'<' {
                // skip any other special alleles
            } else if alt_allele.len() == 1 && ref_allele.len() == 1 {
                // SNV
                push_variant(model::Variant::Snv(alt_allele[0]), i);
            } else if alt_allele.len() == ref_allele.len() {
                // MNV
                push_variant(model::Variant::Mnv(alt_allele.to_vec()), i);
            } else if is_valid_deletion_alleles(ref_allele, alt_allele) {
                push_variant(
                    model::Variant::Deletion((ref_allele.len() - alt_allele.len()) as u64),
                    i,
                );
            } else if is_valid_insertion_alleles(ref_allele, alt_allele) {
                push_variant(
                    model::Variant::Insertion(alt_allele[ref_allele.len()..].to_owned()),
                    i,
                );
            } else {
                // arbitrary replacement
                push_variant(
                    model::Variant::Replacement {
                        ref_allele: ref_allele.to_owned(),
                        alt_allele: alt_allele.to_vec(),
                    },
                    i,
                );
            }
        }
    }

    Ok(variants)
}
