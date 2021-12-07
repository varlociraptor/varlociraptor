// Copyright 2021 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cell::RefCell;
use std::cmp;
use std::ops::Range;
use std::rc::Rc;
use std::sync::Arc;

use anyhow::Result;
use bio::stats::pairhmm::EmissionParameters;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval};

use crate::default_ref_base_emission;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::reference;
use crate::utils::homopolymers::HomopolymerIndelOperation;
use crate::variants::evidence::realignment::pairhmm::{
    ReadEmission, RefBaseEmission, RefBaseVariantEmission, VariantEmission,
};
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::sampling_bias::{ReadSamplingBias, SamplingBias};
use crate::variants::types::{AlleleSupport, MultiLocus, PairedEndEvidence, SingleLocus, Variant};

pub(crate) struct Replacement<R: Realigner> {
    locus: MultiLocus,
    replacement: Rc<Vec<u8>>,
    realigner: RefCell<R>,
    homopolymer_indel_len: Option<i8>,
}

impl<R: Realigner> Replacement<R> {
    pub(crate) fn new(locus: genome::Interval, replacement: Vec<u8>, realigner: R) -> Result<Self> {
        let ref_seq = &realigner.ref_buffer().seq(locus.contig())?;
        let homopolymer_indel_len = HomopolymerIndelOperation::from_text_and_pattern_global(
            &ref_seq[locus.range().start as usize..locus.range().end as usize],
            &replacement,
        )
        .map(|op| op.len());

        Ok(Replacement {
            locus: MultiLocus::new(vec![SingleLocus::new(locus)]),
            replacement: Rc::new(replacement),
            realigner: RefCell::new(realigner),
            homopolymer_indel_len,
        })
    }

    pub(crate) fn locus(&self) -> &SingleLocus {
        &self.locus[0]
    }

    fn ref_len(&self) -> usize {
        (self.locus().range().end - self.locus().range().start) as usize
    }

    fn is_deletion(&self) -> bool {
        self.replacement.len() < self.ref_len()
    }

    fn is_insertion(&self) -> bool {
        self.replacement.len() > self.ref_len()
    }
}

impl<R: Realigner> Realignable for Replacement<R> {
    fn alt_emission_params(
        &self,
        ref_buffer: Arc<reference::Buffer>,
        _: &genome::Interval,
        ref_window: usize,
    ) -> Result<Vec<Box<dyn RefBaseVariantEmission>>> {
        let repl_alt_len = self.replacement.len() as usize;
        let repl_ref_len = self.ref_len();

        let start = self.locus().range().start as usize;

        let ref_seq = ref_buffer.seq(self.locus().contig())?;

        let ref_seq_len = ref_seq.len();

        Ok(vec![Box::new(ReplacementEmissionParams {
            ref_seq,
            ref_offset: start.saturating_sub(ref_window),
            ref_end: cmp::min(start + repl_ref_len + ref_window, ref_seq_len),
            repl_start: start,
            repl_alt_end: start + repl_alt_len,
            repl_alt_len,
            repl_ref_len,
            repl_seq: Rc::clone(&self.replacement),
            is_homopolymer_indel: self.homopolymer_indel_len.is_some(),
        })])
    }
}

impl<R: Realigner> SamplingBias for Replacement<R> {
    fn feasible_bases(
        &self,
        read_len: u64,
        alignment_properties: &AlignmentProperties,
    ) -> Option<u64> {
        let len = self.enclosable_len().unwrap();
        if self.is_insertion() {
            if let Some(maxlen) = alignment_properties.max_ins_cigar_len {
                if len <= (maxlen as u64) {
                    return Some(read_len);
                }
            }
        } else if self.is_deletion() {
            if let Some(maxlen) = alignment_properties.max_del_cigar_len {
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
        let len_diff = self.replacement.len() as isize - self.ref_len() as isize;
        Some(if self.is_insertion() {
            len_diff as u64
        } else if self.is_deletion() {
            (-len_diff) as u64
        } else {
            unreachable!("bug: replacements have to either delete or insert something")
        })
    }
}

impl<R: Realigner> ReadSamplingBias for Replacement<R> {}

impl<R: Realigner> Variant for Replacement<R> {
    type Evidence = PairedEndEvidence;
    type Loci = MultiLocus;

    fn homopolymer_indel_len(&self) -> Option<i8> {
        self.homopolymer_indel_len
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
}

/// Emission parameters for PairHMM over replacement allele.
pub(crate) struct ReplacementEmissionParams {
    ref_seq: Arc<Vec<u8>>,
    ref_offset: usize,
    ref_end: usize,
    repl_start: usize,
    repl_alt_end: usize,
    repl_alt_len: usize,
    repl_ref_len: usize,
    repl_seq: Rc<Vec<u8>>,
    is_homopolymer_indel: bool,
}

impl RefBaseEmission for ReplacementEmissionParams {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        let i_ = i + self.ref_offset;
        if i_ < self.repl_start {
            self.ref_seq[i_]
        } else if i_ >= self.repl_alt_end {
            self.ref_seq[i_ - self.repl_alt_len + self.repl_ref_len]
        } else {
            self.repl_seq[i_ - (self.repl_start)]
        }
    }

    fn variant_homopolymer_ref_range(&self) -> Option<Range<u64>> {
        if self.is_homopolymer_indel {
            Some(self.repl_start as u64..(self.repl_start + self.repl_ref_len) as u64)
        } else {
            None
        }
    }

    #[inline]
    fn len_x(&self) -> usize {
        // The window is shrunken later on with shrink to hit.
        // Therefore it can happen that the window does not cover the entire variant anymore.
        // Hence, we have to consider that the first computation yields 0, in which case we
        // will simply take the second solution.
        let altered_len =
            (self.ref_end - self.ref_offset + self.repl_alt_len).saturating_sub(self.repl_ref_len);
        if altered_len == 0 {
            self.ref_end - self.ref_offset
        } else {
            altered_len
        }
    }

    default_ref_base_emission!();
}

impl VariantEmission for ReplacementEmissionParams {
    fn is_homopolymer_indel(&self) -> bool {
        self.is_homopolymer_indel
    }
}
