// Copyright 2021 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::mem;

use anyhow::Result;

use bio::stats::LogProb;
use itertools::Itertools;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::evidence::observations::read_observation::{
    Evidence, Observable, SingleEndEvidence,
};
use crate::variants::evidence::realignment::Realignable;
use crate::variants::types::{
    AlleleSupport, AlleleSupportBuilder, MultiLocus, PairedEndEvidence, SingleLocus, Variant,
};

use super::ToVariantRepresentation;

pub(crate) trait SingleLocusSingleEndVariant:
    Variant<Loci = SingleLocus, Evidence = SingleEndEvidence>
    + Observable<SingleEndEvidence>
    + ToVariantRepresentation
{
}

impl<V> SingleLocusSingleEndVariant for V where
    V: Variant<Loci = SingleLocus, Evidence = SingleEndEvidence>
        + ToVariantRepresentation
        + Observable<SingleEndEvidence>
{
}

pub(crate) trait SingleLocusPairedEndVariant:
    Variant<Loci = SingleLocus, Evidence = PairedEndEvidence>
    + Observable<PairedEndEvidence>
    + ToVariantRepresentation
{
}

impl<V> SingleLocusPairedEndVariant for V where
    V: Variant<Loci = SingleLocus, Evidence = PairedEndEvidence>
        + ToVariantRepresentation
        + Observable<PairedEndEvidence>
{
}

pub(crate) trait MultiLocusSingleEndVariant:
    Variant<Loci = MultiLocus, Evidence = SingleEndEvidence>
    + Observable<SingleEndEvidence>
    + ToVariantRepresentation
{
}

impl<V> MultiLocusSingleEndVariant for V where
    V: Variant<Loci = MultiLocus, Evidence = SingleEndEvidence>
        + ToVariantRepresentation
        + Observable<SingleEndEvidence>
{
}

pub(crate) trait MultiLocusPairedEndVariant:
    Variant<Loci = MultiLocus, Evidence = PairedEndEvidence>
    + Observable<PairedEndEvidence>
    + ToVariantRepresentation
{
}

impl<V> MultiLocusPairedEndVariant for V where
    V: Variant<Loci = MultiLocus, Evidence = PairedEndEvidence>
        + ToVariantRepresentation
        + Observable<PairedEndEvidence>
{
}

pub(crate) trait HaplotypeVariant: Variant + ToVariantRepresentation {}

#[derive(Default, Getters)]
pub(crate) struct HaplotypeBlock {
    #[getset(get = "pub(crate)")]
    single_locus_single_end_evidence_variants: Vec<Box<dyn SingleLocusSingleEndVariant>>,
    #[getset(get = "pub(crate)")]
    single_locus_paired_end_evidence_variants: Vec<Box<dyn SingleLocusPairedEndVariant>>,
    #[getset(get = "pub(crate)")]
    multi_locus_single_end_evidence_variants: Vec<Box<dyn MultiLocusSingleEndVariant>>,
    #[getset(get = "pub(crate)")]
    multi_locus_paired_end_evidence_variants: Vec<Box<dyn MultiLocusPairedEndVariant>>,
    loci: MultiLocus,
}

impl HaplotypeBlock {
    pub(crate) fn new() -> Self {
        HaplotypeBlock {
            single_locus_single_end_evidence_variants: Vec::new(),
            single_locus_paired_end_evidence_variants: Vec::new(),
            multi_locus_single_end_evidence_variants: Vec::new(),
            multi_locus_paired_end_evidence_variants: Vec::new(),
            loci: MultiLocus::new(Vec::new()),
        }
    }

    pub(crate) fn push_single_locus_single_end_evidence_variant(
        &mut self,
        variant: Box<dyn SingleLocusSingleEndVariant>,
    ) {
        self.loci.push(variant.loci().to_owned());
        self.single_locus_single_end_evidence_variants.push(variant);
    }

    pub(crate) fn push_single_locus_paired_end_evidence_variant(
        &mut self,
        variant: Box<dyn SingleLocusPairedEndVariant>,
    ) {
        self.loci.push(variant.loci().to_owned());
        self.single_locus_paired_end_evidence_variants.push(variant);
    }

    pub(crate) fn push_multi_locus_single_end_evidence_variant(
        &mut self,
        variant: Box<dyn MultiLocusSingleEndVariant>,
    ) {
        self.loci.extend(variant.loci().iter().cloned());
        self.multi_locus_single_end_evidence_variants.push(variant);
    }

    pub(crate) fn push_multi_locus_paired_end_evidence_variant(
        &mut self,
        variant: Box<dyn MultiLocusPairedEndVariant>,
    ) {
        self.loci.extend(variant.loci().iter().cloned());
        self.multi_locus_paired_end_evidence_variants.push(variant);
    }
}

impl Variant for HaplotypeBlock {
    type Evidence = PairedEndEvidence;
    type Loci = MultiLocus;

    fn is_valid_evidence(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> Option<Vec<usize>> {
        let mut valid_indices: Vec<usize> = self
            .single_locus_single_end_evidence_variants
            .iter()
            .enumerate()
            .filter_map(|(i, variant)| {
                if evidence.into_single_end_evidence().iter().any(|evidence| {
                    variant
                        .is_valid_evidence(evidence, alignment_properties)
                        .is_some()
                }) {
                    Some(i)
                } else {
                    None
                }
            })
            .collect();
        let mut locus_offset = self.single_locus_single_end_evidence_variants.len();

        valid_indices.extend(
            self.single_locus_paired_end_evidence_variants
                .iter()
                .enumerate()
                .filter_map(|(i, variant)| {
                    if variant
                        .is_valid_evidence(evidence, alignment_properties)
                        .is_some()
                    {
                        Some(i + locus_offset)
                    } else {
                        None
                    }
                }),
        );
        locus_offset += self.single_locus_paired_end_evidence_variants.len();

        valid_indices.extend(
            self.multi_locus_single_end_evidence_variants
                .iter()
                .enumerate()
                .filter_map(|(i, variant)| {
                    let valid = evidence
                        .into_single_end_evidence()
                        .iter()
                        .filter_map(|evidence| {
                            variant.is_valid_evidence(evidence, alignment_properties)
                        })
                        .flatten()
                        .unique()
                        .collect_vec();
                    let ret = if valid.is_empty() {
                        None
                    } else {
                        Some(
                            valid
                                .into_iter()
                                .map(|idx| idx + locus_offset)
                                .collect_vec(),
                        )
                    };
                    locus_offset += variant.loci().len();
                    ret
                })
                .flatten(),
        );

        valid_indices.extend(
            self.multi_locus_paired_end_evidence_variants
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
                }),
        );

        if !valid_indices.is_empty() {
            Some(valid_indices)
        } else {
            None
        }
    }

    fn loci(&self) -> &Self::Loci {
        &self.loci
    }

    fn allele_support(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
        _alt_variants: &[Box<dyn Realignable>],
    ) -> Result<Option<AlleleSupport>> {
        let support: Vec<Option<AlleleSupport>> = self
            .single_locus_single_end_evidence_variants
            .iter()
            .map(|variant| {
                evidence
                    .into_single_end_evidence()
                    .iter()
                    .map(|evidence| variant.allele_support(evidence, alignment_properties, &[]))
                    .collect_vec()
            })
            .flatten()
            .chain(
                self.single_locus_paired_end_evidence_variants
                    .iter()
                    .map(|variant| variant.allele_support(evidence, alignment_properties, &[])),
            )
            .chain(
                self.multi_locus_single_end_evidence_variants
                    .iter()
                    .map(|variant| {
                        evidence
                            .into_single_end_evidence()
                            .iter()
                            .map(|evidence| {
                                variant.allele_support(evidence, alignment_properties, &[])
                            })
                            .collect_vec()
                    })
                    .flatten(),
            )
            .chain(
                self.multi_locus_paired_end_evidence_variants
                    .iter()
                    .map(|variant| variant.allele_support(evidence, alignment_properties, &[])),
            )
            .collect::<Result<Vec<Option<AlleleSupport>>>>()?;
        let mut support = support
            .into_iter()
            .filter_map(|support| support)
            .collect_vec();
        if support.is_empty() {
            Ok(None)
        } else {
            // if evidence.name() == b"simulated.105" {
            //     dbg!(haplotype_support(&support, true));
            // }
            Ok(Some(haplotype_support(&support)))
        }
    }

    fn prob_sample_alt(
        &self,
        _evidence: &Self::Evidence,
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
        .build()
        .unwrap()
}
