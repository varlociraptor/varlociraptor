// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cell::RefCell;
use std::cmp;
use std::rc::Rc;
use std::sync::Arc;

use anyhow::Result;
use bio::stats::pairhmm::EmissionParameters;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval, AbstractLocus};

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::reference;
use crate::variants::evidence::realignment::pairhmm::{ReadEmission, RefBaseEmission};
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::sampling_bias::{ReadSamplingBias, SamplingBias};
use crate::variants::types::{AlleleSupport, MultiLocus, PairedEndEvidence, SingleLocus, Variant};
use crate::{default_emission, default_ref_base_emission};

pub(crate) struct Insertion {
    locus: MultiLocus,
    ins_seq: Rc<Vec<u8>>,
    realigner: RefCell<Realigner>,
}

impl Insertion {
    pub(crate) fn new(locus: genome::Locus, ins_seq: Vec<u8>, realigner: Realigner) -> Self {
        Insertion {
            locus: MultiLocus::new(vec![SingleLocus::new(genome::Interval::new(
                locus.contig().to_owned(),
                locus.pos()..locus.pos() + 1,
            ))]),
            ins_seq: Rc::new(ins_seq),
            realigner: RefCell::new(realigner),
        }
    }

    pub(crate) fn locus(&self) -> &SingleLocus {
        &self.locus[0]
    }
}

impl<'a> Realignable<'a> for Insertion {
    type EmissionParams = InsertionEmissionParams<'a>;

    fn alt_emission_params(
        &self,
        read_emission_params: Rc<ReadEmission<'a>>,
        ref_buffer: Arc<reference::Buffer>,
        _: &genome::Interval,
        ref_window: usize,
    ) -> Result<InsertionEmissionParams<'a>> {
        let l = self.ins_seq.len() as usize;
        let start = self.locus().range().start as usize;

        let ref_seq = ref_buffer.seq(self.locus().contig())?;

        let ref_seq_len = ref_seq.len();
        Ok(InsertionEmissionParams {
            ref_seq,
            ref_offset: start.saturating_sub(ref_window),
            ref_end: cmp::min(start + l + ref_window, ref_seq_len),
            ins_start: start,
            ins_len: l,
            ins_end: start + l,
            ins_seq: Rc::clone(&self.ins_seq),
            read_emission: read_emission_params,
        })
    }
}

impl SamplingBias for Insertion {
    fn enclosable_len(&self) -> Option<u64> {
        Some(self.ins_seq.len() as u64)
    }
}

impl ReadSamplingBias for Insertion {}

impl Variant for Insertion {
    type Evidence = PairedEndEvidence;
    type Loci = MultiLocus;

    fn is_valid_evidence(&self, evidence: &Self::Evidence) -> Option<Vec<usize>> {
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
    ) -> Result<Option<AlleleSupport>> {
        match evidence {
            PairedEndEvidence::SingleEnd(record) => Ok(Some(
                self.realigner
                    .borrow_mut()
                    .allele_support(record, self.locus.iter(), self)?,
            )),
            PairedEndEvidence::PairedEnd { left, right } => {
                let left_support =
                    self.realigner
                        .borrow_mut()
                        .allele_support(left, self.locus.iter(), self)?;
                let right_support =
                    self.realigner
                        .borrow_mut()
                        .allele_support(right, self.locus.iter(), self)?;

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

/// Emission parameters for PairHMM over insertion allele.
#[derive(Debug)]
pub(crate) struct InsertionEmissionParams<'a> {
    ref_seq: Arc<Vec<u8>>,
    ref_offset: usize,
    ref_end: usize,
    ins_start: usize,
    ins_end: usize,
    ins_len: usize,
    ins_seq: Rc<Vec<u8>>,
    read_emission: Rc<ReadEmission<'a>>,
}

impl<'a> RefBaseEmission for InsertionEmissionParams<'a> {
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

    default_ref_base_emission!();
}

impl<'a> EmissionParameters for InsertionEmissionParams<'a> {
    default_emission!();

    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset + self.ins_len
    }
}
