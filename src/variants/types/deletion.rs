// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cell::RefCell;
use std::cmp;
use std::ops::Range;

use std::sync::Arc;

use anyhow::Result;

use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval};

use crate::default_ref_base_emission;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::reference;
use crate::utils::homopolymers::{extend_homopolymer_stretch, is_homopolymer_seq};
use crate::variants::evidence::realignment::pairhmm::{
    RefBaseEmission, RefBaseVariantEmission, VariantEmission,
};
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::model;
use crate::variants::sampling_bias::{FragmentSamplingBias, ReadSamplingBias, SamplingBias};
use crate::variants::types::{AlleleSupport, MultiLocus, PairedEndEvidence, SingleLocus, Variant};

use super::{IsizeObservable, ToVariantRepresentation};

#[derive(Debug)]
pub(crate) struct Deletion<R: Realigner> {
    locus: SingleLocus,
    fetch_loci: MultiLocus,
    realigner: RefCell<R>,
    homopolymer: Option<Range<u64>>,
    len: u64,
}

impl<R: Realigner> Deletion<R> {
    pub(crate) fn new(locus: genome::Interval, realigner: R) -> Result<Self> {
        let start = locus.range().start;
        let end = locus.range().end;
        let len = end - start;
        let centerpoint = start + (len as f64 / 2.0).round() as u64;
        let contig = locus.contig().to_owned();
        let ref_seq = realigner.ref_buffer().seq(&contig)?;
        let del_seq = &ref_seq[start as usize + 1..end as usize + 1];

        let homopolymer = if is_homopolymer_seq(del_seq) {
            let start = start + 1
                - extend_homopolymer_stretch(
                    del_seq[0],
                    &mut ref_seq[..start as usize + 1].iter().rev(),
                ) as u64;
            let end = end
                + 1
                + extend_homopolymer_stretch(del_seq[0], &mut ref_seq[end as usize + 1..].iter())
                    as u64;
            if end - start > 1 {
                Some(start..end)
            } else {
                None
            }
        } else {
            None
        };

        let fetch_loci = MultiLocus::new(vec![
            SingleLocus::new(genome::Interval::new(contig.clone(), start..start + 1)),
            SingleLocus::new(genome::Interval::new(
                contig.clone(),
                centerpoint..centerpoint + 1,
            )),
            SingleLocus::new(genome::Interval::new(contig, end - 1..end)),
        ]);

        Ok(Deletion {
            locus: SingleLocus::new(locus),
            fetch_loci,
            realigner: RefCell::new(realigner),
            homopolymer,
            len,
        })
    }

    pub(crate) fn centerpoint(&self) -> u64 {
        self.fetch_loci[1].range().start
    }
}

impl<R: Realigner> IsizeObservable for Deletion<R> {}

impl<R: Realigner> SamplingBias for Deletion<R> {
    fn feasible_bases(
        &self,
        read_len: u64,
        alignment_properties: &AlignmentProperties,
    ) -> Option<u64> {
        if let Some(len) = self.enclosable_len() {
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
        Some(self.locus.range().end - self.locus.range().start)
    }
}

impl<R: Realigner> FragmentSamplingBias for Deletion<R> {}
impl<R: Realigner> ReadSamplingBias for Deletion<R> {}

impl<R: Realigner> Realignable for Deletion<R> {
    fn alt_emission_params(
        &self,
        ref_buffer: Arc<reference::Buffer>,
        _: &genome::Interval,
        ref_window: usize,
    ) -> Result<Vec<Box<dyn RefBaseVariantEmission>>> {
        let start = self.locus.range().start as usize;
        let end = self.locus.range().end as usize;
        let ref_seq = ref_buffer.seq(self.locus.contig())?;

        Ok(vec![Box::new(DeletionEmissionParams {
            del_start: start,
            del_len: end - start,
            ref_offset: start.saturating_sub(ref_window),
            ref_end: cmp::min(start + ref_window, ref_seq.len() - self.len() as usize),
            ref_seq,
            homopolymer: self.homopolymer.clone(),
            ref_offset_override: None,
            ref_end_override: None,
        })])
    }
}

impl<R: Realigner> Variant for Deletion<R> {
    fn is_imprecise(&self) -> bool {
        false
    }

    fn homopolymer_indel_len(&self) -> Option<i8> {
        // METHOD: enable HomopolymerArtifact to detect e.g. homopolymer errors due to PCR
        if self.homopolymer.is_some() {
            Some(-(self.len() as i8))
        } else {
            None
        }
    }

    fn is_valid_evidence(
        &self,
        evidence: &PairedEndEvidence,
        alignment_properties: &AlignmentProperties,
    ) -> Option<Vec<usize>> {
        match evidence {
            PairedEndEvidence::SingleEnd(read) => {
                if !self.locus.overlap(read.record(), true).is_none() {
                    Some(vec![0])
                } else {
                    None
                }
            }
            PairedEndEvidence::PairedEnd { left, right } => {
                if alignment_properties.insert_size.is_some() {
                    let right_cigar = right.record().cigar_cached().unwrap();
                    let encloses_centerpoint = (left.record().pos() as u64) < self.centerpoint()
                        && right_cigar.end_pos() as u64 > self.centerpoint();
                    // METHOD: only keep fragments that enclose the centerpoint and have at least one overlapping read.
                    // Importantly, enclosed reads have to be allowed as well. Otherwise, we bias against reference
                    // reads, since they are more unlikely to overlap a breakend and span the centerpoint at the same time,
                    // in particular for large deletions.
                    if encloses_centerpoint
                        && (!self.locus.overlap(left.record(), true).is_none()
                            || !self.locus.overlap(right.record(), true).is_none())
                    {
                        Some(vec![0])
                    } else {
                        None
                    }
                } else if !self.locus.overlap(left.record(), true).is_none()
                    || !self.locus.overlap(right.record(), true).is_none()
                {
                    Some(vec![0])
                } else {
                    None
                }
            }
        }
    }

    /// Return variant loci.
    fn loci(&self) -> &MultiLocus {
        &self.fetch_loci
    }

    fn allele_support(
        &self,
        evidence: &PairedEndEvidence,
        alignment_properties: &AlignmentProperties,
        alt_variants: &[Box<dyn Realignable>],
    ) -> Result<Option<AlleleSupport>> {
        match evidence {
            PairedEndEvidence::SingleEnd(record) => {
                Ok(Some(self.realigner.borrow_mut().allele_support(
                    record.record(),
                    &[&self.locus],
                    self,
                    alt_variants,
                    alignment_properties,
                )?))
            }
            PairedEndEvidence::PairedEnd { left, right } => {
                // METHOD: Extract insert size information for fragments (e.g. read pairs) spanning an indel of interest
                // Here we calculate the product of insert size based and alignment based probabilities.
                // This has the benefit that the calculation automatically checks for consistence between
                // insert size and overlapping alignmnments.
                // This sports the following desirable behavior:
                //
                // * If there is no clear evidence from either the insert size or the alignment, the factors
                //   simply scale because the probabilities of the corresponding type of evidence will be equal.
                // * If reads and fragments agree, 1 stays 1 and 0 stays 0.
                // * If reads and fragments disagree (the promising part!), the other type of evidence will
                //   scale potential false positive probabilities down.
                // * Since there is only one observation per fragment, there is no double counting when
                //   estimating allele frequencies. Before, we had one observation for an overlapping read
                //   and potentially another observation for the corresponding fragment.
                let left_support = self.realigner.borrow_mut().allele_support(
                    left.record(),
                    &[&self.locus],
                    self,
                    alt_variants,
                    alignment_properties,
                )?;
                let right_support = self.realigner.borrow_mut().allele_support(
                    right.record(),
                    &[&self.locus],
                    self,
                    alt_variants,
                    alignment_properties,
                )?;

                let mut support = left_support;
                support.merge(&right_support);

                if alignment_properties.insert_size.is_some() {
                    let isize_support = self.allele_support_isize(
                        left.record(),
                        right.record(),
                        alignment_properties,
                        self.len(),
                    )?;
                    support.merge(&isize_support);
                }

                Ok(Some(support))
            }
        }
    }

    fn prob_sample_alt(
        &self,
        evidence: &PairedEndEvidence,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        match evidence {
            PairedEndEvidence::PairedEnd { left, right } => {
                if alignment_properties.insert_size.is_some() {
                    self.prob_sample_alt_fragment(
                        left.record().seq().len() as u64,
                        right.record().seq().len() as u64,
                        alignment_properties,
                    )
                } else {
                    // METHOD: we do not require the fragment to enclose the variant.
                    // Hence, we treat both reads independently.
                    (self
                        .prob_sample_alt_read(
                            left.record().seq().len() as u64,
                            alignment_properties,
                        )
                        .ln_one_minus_exp()
                        + self
                            .prob_sample_alt_read(
                                right.record().seq().len() as u64,
                                alignment_properties,
                            )
                            .ln_one_minus_exp())
                    .ln_one_minus_exp()
                }
            }
            PairedEndEvidence::SingleEnd(read) => {
                self.prob_sample_alt_read(read.record().seq().len() as u64, alignment_properties)
            }
        }
    }
}

impl<R: Realigner> ToVariantRepresentation for Deletion<R> {
    fn to_variant_representation(&self) -> model::Variant {
        model::Variant::Deletion(self.len)
    }
}

/// Emission parameters for PairHMM over deletion allele.
pub(crate) struct DeletionEmissionParams {
    ref_seq: Arc<Vec<u8>>,
    ref_offset: usize,
    ref_end: usize,
    del_start: usize,
    del_len: usize,
    homopolymer: Option<Range<u64>>,
    ref_offset_override: Option<usize>,
    ref_end_override: Option<usize>,
}

impl RefBaseEmission for DeletionEmissionParams {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        let i_ = i + self.ref_offset;
        if i_ <= self.del_start {
            self.ref_seq[i_]
        } else {
            self.ref_seq[i_ + self.del_len]
        }
    }

    fn variant_homopolymer_ref_range(&self) -> Option<Range<u64>> {
        self.homopolymer.clone()
    }

    fn variant_ref_range(&self) -> Option<Range<u64>> {
        Some(self.del_start as u64..(self.del_start as u64 + self.del_len as u64))
    }

    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset
    }

    default_ref_base_emission!();
}

impl VariantEmission for DeletionEmissionParams {
    fn is_homopolymer_indel(&self) -> bool {
        self.homopolymer.is_some()
    }

    fn alt_vs_ref_len_diff(&self) -> isize {
        -(self.del_len as isize)
    }
}
