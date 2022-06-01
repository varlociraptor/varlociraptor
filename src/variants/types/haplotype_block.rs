// Copyright 2021 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::fmt::Debug;

use anyhow::Result;

use bio::stats::LogProb;
use itertools::Itertools;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::evidence::observations::read_observation::{Observable, SingleEndEvidence};
use crate::variants::evidence::realignment::Realignable;
use crate::variants::model;
use crate::variants::types::{AlleleSupport, MultiLocus, PairedEndEvidence, SingleLocus, Variant};

pub(crate) trait SingleLocusSingleEndVariant:
    Variant<Loci = SingleLocus, Evidence = SingleEndEvidence> + Observable<SingleEndEvidence> + Debug
{
}

pub(crate) trait SingleLocusPairedEndVariant:
    Variant<Loci = SingleLocus, Evidence = PairedEndEvidence> + Observable<PairedEndEvidence> + Debug
{
}

pub(crate) trait MultiLocusSingleEndVariant:
    Variant<Loci = MultiLocus, Evidence = SingleEndEvidence> + Observable<SingleEndEvidence> + Debug
{
}

pub(crate) trait MultiLocusPairedEndVariant:
    Variant<Loci = MultiLocus, Evidence = PairedEndEvidence> + Observable<PairedEndEvidence> + Debug
{
}

#[derive(Debug, Default, Getters)]
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
            let mut haplotype_support = support.pop().unwrap();
            for locus_support in &support {
                haplotype_support.merge(locus_support);
            }
            Ok(Some(haplotype_support))
        }
    }

    fn prob_sample_alt(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        // TODO combine sampling probs of all involved variants, reuse is_valid_evidence information for that
        LogProb::ln_one()
    }

    fn to_variant_representation<'a>(&'a self) -> Box<dyn Iterator<Item = model::Variant> + 'a> {
        Box::new(
            self.single_locus_single_end_evidence_variants
                .iter()
                .map(|variant| variant.to_variant_representation())
                .flatten()
                .chain(
                    self.single_locus_paired_end_evidence_variants
                        .iter()
                        .map(|variant| variant.to_variant_representation())
                        .flatten(),
                )
                .chain(
                    self.multi_locus_single_end_evidence_variants
                        .iter()
                        .map(|variant| variant.to_variant_representation())
                        .flatten(),
                )
                .chain(
                    self.multi_locus_paired_end_evidence_variants
                        .iter()
                        .map(|variant| variant.to_variant_representation())
                        .flatten(),
                ),
        )
    }
}
