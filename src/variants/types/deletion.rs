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
use bio_types::genome::{self, AbstractInterval};
use rust_htslib::bam;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::evidence::insert_size::estimate_insert_size;
use crate::variants::evidence::realignment::pairhmm::{ReadEmission, RefBaseEmission};
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::sampling_bias::{FragmentSamplingBias, ReadSamplingBias, SamplingBias};
use crate::variants::types::{
    AlleleSupport, AlleleSupportBuilder, MultiLocus, PairedEndEvidence, SingleLocus, Variant,
};
use crate::{default_emission, default_ref_base_emission};

pub(crate) struct Deletion {
    locus: SingleLocus,
    fetch_loci: MultiLocus,
    realigner: RefCell<Realigner>,
}

impl Deletion {
    pub(crate) fn new(locus: genome::Interval, realigner: Realigner) -> Self {
        let start = locus.range().start;
        let end = locus.range().end;
        let len = end - start;
        let centerpoint = start + (len as f64 / 2.0).round() as u64;
        let contig = locus.contig().to_owned();

        let fetch_loci = MultiLocus::new(vec![
            SingleLocus(genome::Interval::new(contig.clone(), start..start + 1)),
            SingleLocus(genome::Interval::new(
                contig.clone(),
                centerpoint..centerpoint + 1,
            )),
            SingleLocus(genome::Interval::new(contig, end - 1..end)),
        ]);

        Deletion {
            locus: SingleLocus(locus),
            fetch_loci,
            realigner: RefCell::new(realigner),
        }
    }

    pub(crate) fn centerpoint(&self) -> u64 {
        self.fetch_loci[1].range().start
    }

    pub(crate) fn start(&self) -> u64 {
        self.fetch_loci[0].range().start
    }

    pub(crate) fn end(&self) -> u64 {
        self.fetch_loci[2].range().end
    }

    pub(crate) fn allele_support_isize(
        &self,
        left_record: &bam::Record,
        right_record: &bam::Record,
        alignment_properties: &AlignmentProperties,
    ) -> Result<AlleleSupport> {
        let insert_size = estimate_insert_size(left_record, right_record)?;

        let p_ref = self.isize_pmf(insert_size, 0.0, alignment_properties);
        let p_alt = self.isize_pmf(insert_size, self.len() as f64, alignment_properties);

        if (p_ref == LogProb::ln_zero()
            && !self.is_within_sd(insert_size, self.len() as f64, alignment_properties))
            || (p_alt == LogProb::ln_zero()
                && !self.is_within_sd(insert_size, 0.0, alignment_properties))
        {
            // METHOD: We cannot consider insert size as a reliable estimate here, because it is
            // outside of the numerical resolution for one of the alleles, and not within a
            // standard deviation away from the mean for the other allele.
            Ok(AlleleSupportBuilder::default()
                .prob_ref_allele(LogProb::ln_one())
                .prob_alt_allele(LogProb::ln_one())
                .no_strand_info()
                .build()
                .unwrap())
        } else {
            Ok(AlleleSupportBuilder::default()
                .prob_ref_allele(p_ref)
                .prob_alt_allele(p_alt)
                .no_strand_info()
                .build()
                .unwrap())
        }
    }
}

impl SamplingBias for Deletion {
    fn len(&self) -> u64 {
        self.locus.range().end - self.locus.range().start
    }
}

impl FragmentSamplingBias for Deletion {}
impl ReadSamplingBias for Deletion {}

impl<'a> Realignable<'a> for Deletion {
    type EmissionParams = DeletionEmissionParams<'a>;

    fn alt_emission_params(
        &self,
        read_emission_params: Rc<ReadEmission<'a>>,
        ref_seq: Arc<Vec<u8>>,
        ref_window: usize,
    ) -> DeletionEmissionParams<'a> {
        let start = self.locus.range().start as usize;
        let end = self.locus.range().end as usize;
        DeletionEmissionParams {
            del_start: start,
            del_len: end - start,
            ref_offset: start.saturating_sub(ref_window),
            ref_end: cmp::min(start + ref_window, ref_seq.len() - self.len() as usize),
            ref_seq,
            read_emission: read_emission_params,
        }
    }
}

impl Variant for Deletion {
    type Evidence = PairedEndEvidence;
    type Loci = MultiLocus;

    fn is_valid_evidence(&self, evidence: &Self::Evidence) -> Option<Vec<usize>> {
        match evidence {
            PairedEndEvidence::SingleEnd(read) => {
                if !self.locus.overlap(read, true).is_none() {
                    Some(vec![0])
                } else {
                    None
                }
            }
            PairedEndEvidence::PairedEnd { left, right } => {
                let right_cigar = right.cigar_cached().unwrap();
                let encloses_centerpoint = (left.pos() as u64) < self.centerpoint()
                    && right_cigar.end_pos() as u64 > self.centerpoint();
                // METHOD: only keep fragments that enclose the centerpoint and have at least one overlapping read.
                // Importantly, enclosed reads have to be allowed as well. Otherwise, we bias against reference
                // reads, since they are more unlikely to overlap a breakend and span the centerpoint at the same time,
                // in particular for large deletions.
                if encloses_centerpoint
                    && (!self.locus.overlap(left, true).is_none()
                        || !self.locus.overlap(right, true).is_none())
                {
                    Some(vec![0])
                } else {
                    None
                }
            }
        }
    }

    /// Return variant loci.
    fn loci(&self) -> &Self::Loci {
        &self.fetch_loci
    }

    fn allele_support(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> Result<Option<AlleleSupport>> {
        match evidence {
            PairedEndEvidence::SingleEnd(record) => Ok(Some(
                self.realigner
                    .borrow_mut()
                    .allele_support(record, &self.locus, self)?,
            )),
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
                let left_support =
                    self.realigner
                        .borrow_mut()
                        .allele_support(left, &self.locus, self)?;
                let right_support =
                    self.realigner
                        .borrow_mut()
                        .allele_support(right, &self.locus, self)?;
                let isize_support = self.allele_support_isize(left, right, alignment_properties)?;

                let mut support = left_support;
                support.merge(&right_support).merge(&isize_support);

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
            PairedEndEvidence::PairedEnd { left, right } => self.prob_sample_alt_fragment(
                left.seq().len() as u64,
                right.seq().len() as u64,
                alignment_properties,
            ),
            PairedEndEvidence::SingleEnd(read) => {
                self.prob_sample_alt_read(read.seq().len() as u64, alignment_properties)
            }
        }
    }
}

/// Emission parameters for PairHMM over deletion allele.
#[derive(Debug)]
pub(crate) struct DeletionEmissionParams<'a> {
    ref_seq: Arc<Vec<u8>>,
    ref_offset: usize,
    ref_end: usize,
    del_start: usize,
    del_len: usize,
    read_emission: Rc<ReadEmission<'a>>,
}

impl<'a> RefBaseEmission for DeletionEmissionParams<'a> {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        let i_ = i + self.ref_offset;
        if i_ <= self.del_start {
            self.ref_seq[i_]
        } else {
            self.ref_seq[i_ + self.del_len]
        }
    }

    default_ref_base_emission!();
}

impl<'a> EmissionParameters for DeletionEmissionParams<'a> {
    default_emission!();

    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset
    }
}
