// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cell::RefCell;
use std::ops::Range;
use std::rc::Rc;
use std::sync::Arc;
use std::{cmp, iter};

use anyhow::Result;

use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval, AbstractLocus};

use crate::default_ref_base_emission;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::reference;
use crate::utils::homopolymers::{extend_homopolymer_stretch, is_homopolymer_seq};
use crate::variants::evidence::realignment::pairhmm::{
    RefBaseEmission, RefBaseVariantEmission, VariantEmission,
};
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::model;
use crate::variants::sampling_bias::{ReadSamplingBias, SamplingBias};
use crate::variants::types::{AlleleSupport, MultiLocus, PairedEndEvidence, SingleLocus, Variant};

pub(crate) struct Insertion<R: Realigner> {
    locus: MultiLocus,
    ins_seq: Rc<Vec<u8>>,
    realigner: RefCell<R>,
    homopolymer: Option<Range<u64>>,
}

impl<R: Realigner> Insertion<R> {
    pub(crate) fn new(locus: genome::Locus, ins_seq: Vec<u8>, realigner: R) -> Result<Self> {
        let start = locus.pos() as usize;
        let ref_seq = &realigner.ref_buffer().seq(locus.contig())?;

        let homopolymer = if is_homopolymer_seq(&ins_seq) {
            let end = (start
                + ins_seq.len()
                + extend_homopolymer_stretch(ins_seq[0], &mut ref_seq[start + 1..].iter()))
                as u64;
            let start = start as u64 + 1
                - extend_homopolymer_stretch(ins_seq[0], &mut ref_seq[..start + 1].iter().rev())
                    as u64;
            Some(start..end)
        } else {
            None
        };

        Ok(Insertion {
            locus: MultiLocus::new(vec![SingleLocus::new(genome::Interval::new(
                locus.contig().to_owned(),
                locus.pos()..locus.pos() + 1,
            ))]),
            ins_seq: Rc::new(ins_seq),
            realigner: RefCell::new(realigner),
            homopolymer,
        })
    }

    pub(crate) fn locus(&self) -> &SingleLocus {
        &self.locus[0]
    }
}

impl<R: Realigner> Realignable for Insertion<R> {
    fn alt_emission_params(
        &self,
        ref_buffer: Arc<reference::Buffer>,
        _: &genome::Interval,
        ref_window: usize,
    ) -> Result<Vec<Box<dyn RefBaseVariantEmission>>> {
        let l = self.ins_seq.len() as usize;
        let start = self.locus().range().start as usize;

        let ref_seq = ref_buffer.seq(self.locus().contig())?;

        let ref_seq_len = ref_seq.len();
        Ok(vec![Box::new(InsertionEmissionParams {
            ref_seq,
            ref_offset: start.saturating_sub(ref_window),
            ref_end: cmp::min(start + l + ref_window, ref_seq_len),
            ins_start: start,
            ins_len: l,
            ins_end: start + l,
            ins_seq: Rc::clone(&self.ins_seq),
            homopolymer: self.homopolymer.clone(),
        })])
    }
}

impl<R: Realigner> SamplingBias for Insertion<R> {
    fn feasible_bases(
        &self,
        read_len: u64,
        alignment_properties: &AlignmentProperties,
    ) -> Option<u64> {
        if let Some(len) = self.enclosable_len() {
            if let Some(maxlen) = alignment_properties.max_ins_cigar_len {
                if len <= (maxlen as u64) {
                    return Some(read_len);
                }
            }
        }
        alignment_properties
            .frac_max_softclip
            .map(|maxfrac| (read_len as f64 * maxfrac) as u64)
    }

    fn enclosable_len(&self) -> Option<u64> {
        Some(self.ins_seq.len() as u64)
    }
}

impl<R: Realigner> ReadSamplingBias for Insertion<R> {}

impl<R: Realigner> Variant for Insertion<R> {
    type Evidence = PairedEndEvidence;
    type Loci = MultiLocus;

    fn homopolymer_indel_len(&self) -> Option<i8> {
        if self.homopolymer.is_some() {
            Some(self.ins_seq.len() as i8)
        } else {
            None
        }
    }

    fn is_valid_evidence(
        &self,
        evidence: &Self::Evidence,
        _: &AlignmentProperties,
    ) -> Option<Vec<usize>> {
        if match evidence {
            PairedEndEvidence::SingleEnd(read) => !self.locus().overlap(read, true).is_none(),
            PairedEndEvidence::PairedEnd { left, right } => {
                !self.locus().overlap(left, true).is_none()
                    || !self.locus().overlap(right, true).is_none()
            }
        } {
            Some(vec![0])
        } else {
            None
        }
    }

    /// Return variant loci.
    fn loci(&self) -> &Self::Loci {
        &self.locus
    }

    /// Calculate probability for alt and reference allele.
    fn allele_support(
        &self,
        evidence: &Self::Evidence,
        _alignment_properties: &AlignmentProperties,
        alt_variants: &[Box<dyn Realignable>],
    ) -> Result<Option<AlleleSupport>> {
        match evidence {
            PairedEndEvidence::SingleEnd(record) => {
                Ok(Some(self.realigner.borrow_mut().allele_support(
                    record,
                    self.locus.iter(),
                    self,
                    alt_variants,
                )?))
            }
            PairedEndEvidence::PairedEnd { left, right } => {
                let left_support = self.realigner.borrow_mut().allele_support(
                    left,
                    self.locus.iter(),
                    self,
                    alt_variants,
                )?;
                let right_support = self.realigner.borrow_mut().allele_support(
                    right,
                    self.locus.iter(),
                    self,
                    alt_variants,
                )?;

                let mut support = left_support;

                support.merge(&right_support);

                Ok(Some(support))
            }
        }
    }

    fn prob_sample_alt(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        match evidence {
            PairedEndEvidence::PairedEnd { left, right } => {
                // METHOD: we do not require the fragment to enclose the variant.
                // Hence, we treat both reads independently.
                (self
                    .prob_sample_alt_read(left.seq().len() as u64, alignment_properties)
                    .ln_one_minus_exp()
                    + self
                        .prob_sample_alt_read(right.seq().len() as u64, alignment_properties)
                        .ln_one_minus_exp())
                .ln_one_minus_exp()
            }
            PairedEndEvidence::SingleEnd(read) => {
                self.prob_sample_alt_read(read.seq().len() as u64, alignment_properties)
            }
        }
    }

    fn to_variant_representation(&self) -> Box<dyn Iterator<Item = model::Variant>> {
        Box::new(iter::once(model::Variant::Insertion(self.ins_seq.to_vec())))
    }
}

/// Emission parameters for PairHMM over insertion allele.
pub(crate) struct InsertionEmissionParams {
    ref_seq: Arc<Vec<u8>>,
    ref_offset: usize,
    ref_end: usize,
    ins_start: usize,
    ins_end: usize,
    ins_len: usize,
    ins_seq: Rc<Vec<u8>>,
    homopolymer: Option<Range<u64>>,
}

impl RefBaseEmission for InsertionEmissionParams {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        let i_ = i + self.ref_offset;
        if i_ <= self.ins_start {
            self.ref_seq[i_]
        } else if i_ > self.ins_end {
            self.ref_seq[i_ - self.ins_len]
        } else {
            self.ins_seq[i_ - (self.ins_start + 1)]
        }
    }

    fn variant_homopolymer_ref_range(&self) -> Option<Range<u64>> {
        self.homopolymer.clone()
    }

    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset + self.ins_len
    }

    default_ref_base_emission!();
}

impl VariantEmission for InsertionEmissionParams {
    fn is_homopolymer_indel(&self) -> bool {
        self.homopolymer.is_some()
    }
}
