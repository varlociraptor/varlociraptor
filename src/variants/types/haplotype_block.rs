// Copyright 2021 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::mem;

use anyhow::Result;

use bio::stats::LogProb;
use itertools::Itertools;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::evidence::observations::read_observation::ReadObservable;
use crate::variants::evidence::realignment::edit_distance::EditDistance;
use crate::variants::evidence::realignment::Realignable;
use crate::variants::types::{
    AlleleSupport, AlleleSupportBuilder, Evidence, MultiLocus, ReadVariant,
};

use super::ToVariantRepresentation;

pub(crate) trait HaplotypeVariant:
    ReadVariant + ReadObservable + ToVariantRepresentation
{
}

impl<V> HaplotypeVariant for V where V: ReadVariant + ReadObservable + ToVariantRepresentation {}

#[derive(Default, Getters)]
pub(crate) struct HaplotypeBlock {
    #[getset(get = "pub(crate)")]
    variants: Vec<Box<dyn HaplotypeVariant>>,
    loci: MultiLocus,
}

impl HaplotypeBlock {
    pub(crate) fn new() -> Self {
        HaplotypeBlock {
            variants: Vec::new(),
            loci: MultiLocus::new(Vec::new()),
        }
    }

    pub(crate) fn push_variant(&mut self, variant: Box<dyn HaplotypeVariant>) {
        self.loci.extend(variant.loci().iter().cloned());
        self.variants.push(variant);
    }
}

impl ReadVariant for HaplotypeBlock {
    fn is_imprecise(&self) -> bool {
        false
    }

    fn is_valid_evidence(
        &self,
        evidence: &Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> Option<Vec<usize>> {
        let mut locus_offset = 0;
        let valid_indices: Vec<usize> = self
            .variants
            .iter()
            .enumerate()
            .filter_map(|(i, variant)| {
                let ret = if variant
                    .is_valid_evidence(evidence, alignment_properties)
                    .is_some()
                {
                    Some(i + locus_offset)
                } else {
                    None
                };
                locus_offset += variant.loci().len();
                ret
            })
            .collect();

        if !valid_indices.is_empty() {
            Some(valid_indices)
        } else {
            None
        }
    }

    fn loci(&self) -> &MultiLocus {
        &self.loci
    }

    fn allele_support(
        &self,
        evidence: &Evidence,
        alignment_properties: &AlignmentProperties,
        _alt_variants: &[Box<dyn Realignable>],
    ) -> Result<Option<AlleleSupport>> {
        let support: Vec<Option<_>> = self
            .variants
            .iter()
            .map(|variant| variant.allele_support(evidence, alignment_properties, &[]))
            .collect::<Result<_>>()?;
        let support = support.into_iter().flatten().collect_vec();
        if support.is_empty() {
            Ok(None)
        } else {
            Ok(Some(haplotype_support(&support)))
        }
    }

    fn prob_sample_alt(
        &self,
        _evidence: &Evidence,
        _alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        // TODO combine sampling probs of all involved variants, reuse is_valid_evidence information for that
        LogProb::ln_one()
    }
}

fn haplotype_support(variant_supports: &[AlleleSupport]) -> AlleleSupport {
    let prob_alt_allele = variant_supports
        .iter()
        .map(|support| support.prob_alt_allele)
        .sum();

    let mut third_allele_evidence: Option<EditDistance> = None;
    for variant_support in variant_supports {
        if let Some(ref other_dist) = variant_support.third_allele_evidence {
            if let Some(ref mut dist) = third_allele_evidence {
                dist.update(other_dist);
            } else {
                third_allele_evidence = variant_support.third_allele_evidence;
            }
        }
    }

    // METHOD: calculate prob_ref via dynamic programming:
    // any combination of ref supports should suffice

    // init values
    let n = variant_supports.len();
    let mut current = vec![LogProb::ln_zero(); n + 1];
    let mut last = vec![LogProb::ln_zero(); n + 1];
    // initially, there is no ref allele on the haplotype yet
    last[0] = LogProb::ln_one();

    for variant_support in variant_supports {
        // still no ref allele
        current[0] = last[0] + variant_support.prob_alt_allele();
        for n_ref_allele in 1..=n {
            current[n_ref_allele] = (
                // one more ref allele
                last[n_ref_allele - 1] + variant_support.prob_ref_allele()
            )
            .ln_add_exp(
                // same number of ref alleles as before (this variant must be alt allele hence)
                last[n_ref_allele] + variant_support.prob_alt_allele(),
            );
        }
        mem::swap(&mut current, &mut last);
    }

    let prob_ref_allele = LogProb::ln_sum_exp(&last[1..]);

    AlleleSupportBuilder::default()
        .prob_alt_allele(prob_alt_allele)
        .prob_ref_allele(prob_ref_allele)
        .strand(variant_supports[0].strand())
        .third_allele_evidence(third_allele_evidence)
        .build()
        .unwrap()
}
