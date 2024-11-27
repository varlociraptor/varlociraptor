// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;
use std::collections::BTreeMap;
use std::rc::Rc;

use anyhow::Result;
use bio::stats::{LogProb, PHREDProb};
use bio_types::genome::{self, AbstractInterval};
use rust_htslib::bam;
use vec_map::VecMap;

use crate::errors::Error;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::utils::homopolymers::HomopolymerErrorModel;
use crate::utils::PROB_05;
use crate::variants::evidence::observations::read_observation::{
    Evidence, ExtendedRecord, Observable, PairedEndEvidence, ReadObservation, SingleEndEvidence,
    Strand,
};
use crate::variants::sample;

pub(crate) mod breakends;
pub(crate) mod deletion;
pub(crate) mod duplication;
pub(crate) mod haplotype_block;
pub(crate) mod insertion;
pub(crate) mod inversion;
pub(crate) mod methylation;
pub(crate) mod mnv;
pub(crate) mod none;
pub(crate) mod replacement;
pub(crate) mod snv;

pub(crate) use deletion::Deletion;
pub(crate) use duplication::Duplication;
pub(crate) use insertion::Insertion;
pub(crate) use inversion::Inversion;
pub(crate) use methylation::Methylation;
pub(crate) use mnv::Mnv;
pub(crate) use none::None;
pub(crate) use replacement::Replacement;
pub(crate) use snv::Snv;

use super::evidence::insert_size::estimate_insert_size;
use super::evidence::observations::fragment_id_factory::FragmentIdFactory;
use super::evidence::realignment::edit_distance::EditDistance;
use super::evidence::realignment::Realignable;
use super::model;
use super::sampling_bias::FragmentSamplingBias;

#[derive(Debug, CopyGetters, Getters, Builder)]
pub(crate) struct AlleleSupport {
    prob_ref_allele: LogProb,
    prob_alt_allele: LogProb,
    #[getset(get_copy = "pub")]
    strand: Strand,
    #[builder(default)]
    #[getset(get_copy = "pub")]
    read_position: Option<u32>,
    #[builder(default)]
    #[getset(get_copy = "pub")]
    homopolymer_indel_len: Option<i8>,
    #[getset(get_copy = "pub")]
    third_allele_evidence: Option<EditDistance>,
}

impl AlleleSupport {
    fn both_alleles_impossible(&self) -> bool {
        *self.prob_ref_allele == f64::NEG_INFINITY && *self.prob_alt_allele == f64::NEG_INFINITY
    }

    pub(crate) fn prob_ref_allele(&self) -> LogProb {
        if self.both_alleles_impossible() {
            *PROB_05
        } else {
            self.prob_ref_allele
        }
    }

    pub(crate) fn prob_alt_allele(&self) -> LogProb {
        if self.both_alleles_impossible() {
            *PROB_05
        } else {
            self.prob_alt_allele
        }
    }

    /// METHOD: This is an estimate of the allele likelihood at the true location in case
    /// the read is mismapped. The value has to be approximately in the range of prob_alt
    /// and prob_ref. Otherwise it could cause numerical problems, by dominating the
    /// likelihood such that subtle differences in allele frequencies become numercically
    /// invisible in the resulting likelihood.
    pub(crate) fn prob_missed_allele(&self) -> LogProb {
        self.prob_ref_allele().ln_add_exp(self.prob_alt_allele()) - LogProb(2.0_f64.ln())
    }

    pub(crate) fn merge(&mut self, other: &AlleleSupport) -> &mut Self {
        // TODO set read position to None if both allele supports have one
        self.prob_ref_allele += other.prob_ref_allele;
        self.prob_alt_allele += other.prob_alt_allele;

        if self.strand == Strand::None {
            self.strand = other.strand;
            self.homopolymer_indel_len = other.homopolymer_indel_len;
        } else if other.strand != Strand::None {
            if self.strand != other.strand {
                self.strand = Strand::Both;
            }
            self.homopolymer_indel_len =
                match (self.homopolymer_indel_len, other.homopolymer_indel_len) {
                    (Some(indel_len), Some(_other_indel_len)) => Some(indel_len), // just keep one of them
                    (Some(indel_len), None) => Some(indel_len),
                    (None, Some(indel_len)) => Some(indel_len),
                    (None, None) => None,
                }
        }

        match (
            &mut self.third_allele_evidence,
            &other.third_allele_evidence,
        ) {
            (Some(edit_dist), Some(other_edit_dist)) => edit_dist.update(other_edit_dist),
            (None, Some(other_edit_dist)) => self.third_allele_evidence = Some(*other_edit_dist),
            (Some(_), None) => (),
            (None, None) => (),
        }

        self
    }
}

pub(crate) trait Variant {
    type Evidence: Evidence;
    type Loci: Loci;

    fn is_imprecise(&self) -> bool;

    /// Determine whether the evidence is suitable to assessing probabilities
    /// (i.e. overlaps the variant in the right way).
    ///
    /// # Returns
    ///
    /// The index of the loci for which this evidence is valid, `None` if invalid.
    fn is_valid_evidence(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> Option<Vec<usize>>;

    /// Return variant loci.
    fn loci(&self) -> &Self::Loci;

    /// Calculate probability for alt and reference allele.
    fn allele_support(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
        alt_variants: &[Box<dyn Realignable>],
    ) -> Result<Option<AlleleSupport>>;

    /// Calculate probability to sample a record length like the given one from the alt allele.
    fn prob_sample_alt(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb;

    /// Return the homopolymer indel len of the variant, if any.
    fn homopolymer_indel_len(&self) -> Option<i8> {
        None
    }
}

pub(crate) trait IsizeObservable: Variant + FragmentSamplingBias {
    fn allele_support_isize(
        &self,
        left_record: &bam::Record,
        right_record: &bam::Record,
        alignment_properties: &AlignmentProperties,
        alt_del_len: u64,
    ) -> Result<AlleleSupport> {
        if alignment_properties.insert_size.unwrap().sd == 0.0 {
            return Err(Error::UnrealisticIsizeSd.into());
        }

        let insert_size = estimate_insert_size(left_record, right_record)?;

        // TODO: sum over all possible lens if there are multiple ones and use some reasonable prior
        // or take the best p_alt.
        let p_ref = self.isize_pmf(insert_size, 0.0, alignment_properties);
        let p_alt = self.isize_pmf(insert_size, alt_del_len as f64, alignment_properties);

        if (p_ref == LogProb::ln_zero()
            && !self.is_within_sd(insert_size, alt_del_len as f64, alignment_properties))
            || (p_alt == LogProb::ln_zero()
                && !self.is_within_sd(insert_size, 0.0, alignment_properties))
        {
            // METHOD: We cannot consider insert size as a reliable estimate here, because it is
            // outside of the numerical resolution for one of the alleles, and not within a
            // standard deviation away from the mean for the other allele.
            Ok(AlleleSupportBuilder::default()
                .prob_ref_allele(LogProb::ln_one())
                .prob_alt_allele(LogProb::ln_one())
                .strand(Strand::None)
                .third_allele_evidence(None)
                .build()
                .unwrap())
        } else {
            Ok(AlleleSupportBuilder::default()
                .prob_ref_allele(p_ref)
                .prob_alt_allele(p_alt)
                .strand(Strand::None)
                .third_allele_evidence(None)
                .build()
                .unwrap())
        }
    }

    fn len(&self) -> u64 {
        self.enclosable_len().unwrap()
    }
}

pub(crate) trait ToVariantRepresentation {
    fn to_variant_representation(&self) -> model::Variant;
}

impl<V> Observable<SingleEndEvidence, SingleLocus> for V
where
    V: Variant<Evidence = SingleEndEvidence, Loci = SingleLocus>,
{
    fn prob_mapping(&self, evidence: &SingleEndEvidence) -> LogProb {
        let prob_mismapping = LogProb::from(PHREDProb(evidence.mapq() as f64));
        prob_mismapping.ln_one_minus_exp()
    }

    fn min_mapq(&self, evidence: &SingleEndEvidence) -> u8 {
        evidence.mapq()
    }

    fn extract_observations(
        &self,
        buffer: &mut sample::RecordBuffer,
        alignment_properties: &mut AlignmentProperties,
        max_depth: usize,
        alt_variants: &[Box<dyn Realignable>],
        observation_id_factory: &mut Option<&mut FragmentIdFactory>,
    ) -> Result<Vec<ReadObservation>> {
        let locus = self.loci();
        buffer.fetch(locus, false)?;

        let homopolymer_error_model = HomopolymerErrorModel::new(self, alignment_properties);

        let candidates: Vec<_> = buffer
            .iter()
            .filter_map(|record| {
                // METHOD: First, we check whether the record contains an indel in the cigar.
                // We store the maximum indel size to update the global estimates, in case
                // it is larger in this region.
                alignment_properties.update_max_cigar_ops_len(record.as_ref(), false);

                let evidence = SingleEndEvidence::new(record);
                if self
                    .is_valid_evidence(&evidence, alignment_properties)
                    .is_some()
                {
                    Some(evidence)
                } else {
                    None
                }
            })
            .collect();

        let mut subsampler = sample::SubsampleCandidates::new(max_depth, candidates.len());

        let mut observations = Vec::new();
        for evidence in candidates {
            if subsampler.keep() {
                if let Some(obs) = self.evidence_to_observation(
                    &evidence,
                    alignment_properties,
                    &homopolymer_error_model,
                    alt_variants,
                    observation_id_factory,
                )? {
                    observations.push(obs);
                }
            }
        }
        Ok(observations)
    }
}

impl<V> Observable<PairedEndEvidence, MultiLocus> for V
where
    V: Variant<Evidence = PairedEndEvidence, Loci = MultiLocus>,
{
    fn prob_mapping(&self, evidence: &PairedEndEvidence) -> LogProb {
        let prob = |record: &bam::Record| LogProb::from(PHREDProb(record.mapq() as f64));
        match evidence {
            PairedEndEvidence::SingleEnd(record) => prob(record.record()).ln_one_minus_exp(),
            PairedEndEvidence::PairedEnd { left, right } => {
                // METHOD: take maximum of the (log-spaced) mapping quality of the left and the right read.
                // In BWA, MAPQ is influenced by the mate, hence they are not independent
                // and we can therefore not multiply them. By taking the maximum, we
                // make a conservative choice (since 1-mapq is the mapping probability).
                let mut p = prob(left.record());
                let mut q = prob(right.record());
                if p < q {
                    std::mem::swap(&mut p, &mut q);
                }
                p.ln_one_minus_exp()
            }
        }
    }

    fn min_mapq(&self, evidence: &PairedEndEvidence) -> u8 {
        match evidence {
            PairedEndEvidence::SingleEnd(record) => record.record().mapq(),
            PairedEndEvidence::PairedEnd { left, right } => {
                cmp::min(left.record().mapq(), right.record().mapq())
            }
        }
    }

    fn extract_observations(
        &self,
        buffer: &mut sample::RecordBuffer,
        alignment_properties: &mut AlignmentProperties,
        max_depth: usize,
        alt_variants: &[Box<dyn Realignable>],
        observation_id_factory: &mut Option<&mut FragmentIdFactory>,
    ) -> Result<Vec<ReadObservation>> {
        // We cannot use a hash function here because candidates have to be considered
        // in a deterministic order. Otherwise, subsampling high-depth regions will result
        // in slightly different probabilities each time.
        let mut candidate_records = BTreeMap::new();

        let homopolymer_error_model = HomopolymerErrorModel::new(self, alignment_properties);

        let mut fetches = buffer.build_fetches(true);
        for locus in self.loci().iter() {
            fetches.push(locus);
        }
        for interval in fetches.iter() {
            // Fetch intervals cannot overlap. This is ensured by the way they are built.
            buffer.fetch(interval, true)?;

            for record in buffer.iter() {
                // METHOD: First, we check whether the record contains an indel in the cigar.
                // We store the maximum indel size to update the global estimates, in case
                // it is larger in this region.
                alignment_properties.update_max_cigar_ops_len(record.as_ref(), false);

                // We look at the whole fragment at once.

                // TODO move this text to the right place:
                // We ensure fair sampling by checking if the whole fragment overlaps the
                // centerpoint. Only taking the internal segment would not be fair,
                // because then the second read of reference fragments tends to cross
                // the centerpoint and the fragment would be discarded.
                // The latter would not happen for alt (deletion) fragments, because the second
                // read would map right of the variant in that case.

                // We always choose the leftmost and the rightmost alignment, thereby also
                // considering supplementary alignments.

                if !candidate_records.contains_key(record.qname()) {
                    // this is the first (primary or supplementary alignment in the pair
                    candidate_records.insert(record.qname().to_owned(), Candidate::new(record));
                } else if let Some(candidate) = candidate_records.get_mut(record.qname()) {
                    // this is either the last alignment or one in the middle
                    if (candidate.left.is_first_in_template() && record.is_first_in_template())
                        && (candidate.left.is_last_in_template() && record.is_last_in_template())
                    {
                        // Ignore another partial alignment right of the first.
                        continue;
                    }
                    // replace right record (we seek for the rightmost (partial) alignment)
                    candidate.right = Some(record);
                }
            }
        }

        let mut candidates = Vec::new();
        let mut locus_depth = VecMap::new();
        let mut push_evidence = |evidence: PairedEndEvidence, idx| {
            candidates.push(evidence);
            for i in idx {
                let count = locus_depth.entry(i).or_insert(0);
                *count += 1;
            }
        };

        for candidate in candidate_records.values() {
            if let Some(ref right) = candidate.right {
                if candidate.left.mapq() == 0 || right.mapq() == 0 {
                    // Ignore pairs with ambiguous alignments.
                    // The statistical model does not consider them anyway.
                    continue;
                }
                let evidence = PairedEndEvidence::PairedEnd {
                    // buffer.get_methylation_probs returns None if we do not deal with PacBio or Nanopore methylation
                    left: ExtendedRecord::new(
                        candidate.left.to_owned(),
                        buffer.get_methylation_probs(candidate.left.to_owned()),
                    ),
                    right: ExtendedRecord::new(
                        right.to_owned(),
                        buffer.get_methylation_probs(right.to_owned()),
                    ),
                };
                if let Some(idx) = self.is_valid_evidence(&evidence, alignment_properties) {
                    push_evidence(evidence, idx);
                }
            } else {
                // this is a single alignment with unmapped mate or mate outside of the
                // region of interest
                let evidence = PairedEndEvidence::SingleEnd(ExtendedRecord::new(
                    candidate.left.to_owned(),
                    buffer.get_methylation_probs(candidate.left.to_owned()),
                ));
                if let Some(idx) = self.is_valid_evidence(&evidence, alignment_properties) {
                    push_evidence(evidence, idx);
                }
            }
        }

        // METHOD: if all loci exceed the maximum depth, we subsample the evidence.
        // We cannot decide this per locus, because we risk adding more biases if loci have different alt allele sampling biases.
        let subsample = locus_depth.values().all(|depth| *depth > max_depth);
        let mut subsampler = sample::SubsampleCandidates::new(max_depth, candidates.len());

        let mut observations = Vec::new();
        for evidence in &candidates {
            if !subsample || subsampler.keep() {
                if let Some(obs) = self.evidence_to_observation(
                    evidence,
                    alignment_properties,
                    &homopolymer_error_model,
                    alt_variants,
                    observation_id_factory,
                )? {
                    observations.push(obs);
                }
            }
        }

        Ok(observations)
    }
}

impl<V> Observable<PairedEndEvidence, SingleLocus> for V
where
    V: Variant<Evidence = PairedEndEvidence, Loci = SingleLocus>,
{
    fn prob_mapping(&self, evidence: &PairedEndEvidence) -> LogProb {
        let prob = |record: &bam::Record| LogProb::from(PHREDProb(record.mapq() as f64));
        match evidence {
            PairedEndEvidence::SingleEnd(record) => prob(record.record()).ln_one_minus_exp(),
            PairedEndEvidence::PairedEnd { left, right } => {
                // METHOD: take maximum of the (log-spaced) mapping quality of the left and the right read.
                // In BWA, MAPQ is influenced by the mate, hence they are not independent
                // and we can therefore not multiply them. By taking the maximum, we
                // make a conservative choice (since 1-mapq is the mapping probability).
                let mut p = prob(left.record());
                let mut q = prob(right.record());
                if p < q {
                    std::mem::swap(&mut p, &mut q);
                }
                p.ln_one_minus_exp()
            }
        }
    }

    fn min_mapq(&self, evidence: &PairedEndEvidence) -> u8 {
        match evidence {
            PairedEndEvidence::SingleEnd(record) => record.record().mapq(),
            PairedEndEvidence::PairedEnd { left, right } => {
                cmp::min(left.record().mapq(), right.record().mapq())
            }
        }
    }

    fn extract_observations(
        &self,
        buffer: &mut sample::RecordBuffer,
        alignment_properties: &mut AlignmentProperties,
        max_depth: usize,
        alt_variants: &[Box<dyn Realignable>],
        observation_id_factory: &mut Option<&mut FragmentIdFactory>,
    ) -> Result<Vec<ReadObservation>> {
        // We cannot use a hash function here because candidates have to be considered
        // in a deterministic order. Otherwise, subsampling high-depth regions will result
        // in slightly different probabilities each time.
        let mut candidate_records = BTreeMap::new();

        let locus = self.loci();
        buffer.fetch(locus, true)?;
        let homopolymer_error_model = HomopolymerErrorModel::new(self, alignment_properties);
        for record in buffer.iter() {
            // METHOD: First, we check whether the record contains an indel in the cigar.
            // We store the maximum indel size to update the global estimates, in case
            // it is larger in this region.
            alignment_properties.update_max_cigar_ops_len(record.as_ref(), false);

            // We look at the whole fragment at once.

            // TODO move this text to the right place:
            // We ensure fair sampling by checking if the whole fragment overlaps the
            // centerpoint. Only taking the internal segment would not be fair,
            // because then the second read of reference fragments tends to cross
            // the centerpoint and the fragment would be discarded.
            // The latter would not happen for alt (deletion) fragments, because the second
            // read would map right of the variant in that case.

            // We always choose the leftmost and the rightmost alignment, thereby also
            // considering supplementary alignments.

            if !candidate_records.contains_key(record.qname()) {
                // this is the first (primary or supplementary alignment in the pair
                candidate_records.insert(record.qname().to_owned(), Candidate::new(record));
            } else if let Some(candidate) = candidate_records.get_mut(record.qname()) {
                // this is either the last alignment or one in the middle
                if (candidate.left.is_first_in_template() && record.is_first_in_template())
                    && (candidate.left.is_last_in_template() && record.is_last_in_template())
                {
                    // Ignore another partial alignment right of the first.
                    continue;
                }
                // replace right record (we seek for the rightmost (partial) alignment)
                candidate.right = Some(record);
            }
        }

        let mut candidates = Vec::new();
        let mut locus_depth = VecMap::new();
        let mut push_evidence = |evidence: PairedEndEvidence, idx| {
            candidates.push(evidence);
            for i in idx {
                let count = locus_depth.entry(i).or_insert(0);
                *count += 1;
            }
        };

        for candidate in candidate_records.values() {
            if let Some(ref right) = candidate.right {
                if candidate.left.mapq() == 0 || right.mapq() == 0 {
                    // Ignore pairs with ambiguous alignments.
                    // The statistical model does not consider them anyway.
                    continue;
                }
                let evidence = PairedEndEvidence::PairedEnd {
                    left: ExtendedRecord::new(
                        candidate.left.to_owned(),
                        buffer.get_methylation_probs(candidate.left.to_owned()),
                    ),
                    right: ExtendedRecord::new(
                        right.to_owned(),
                        buffer.get_methylation_probs(right.to_owned()),
                    ),
                };
                if let Some(idx) = self.is_valid_evidence(&evidence, alignment_properties) {
                    push_evidence(evidence, idx);
                }
            } else {
                // this is a single alignment with unmapped mate or mate outside of the
                // region of interest
                let evidence = PairedEndEvidence::SingleEnd(ExtendedRecord::new(
                    candidate.left.to_owned(),
                    buffer.get_methylation_probs(candidate.left.to_owned()),
                ));
                if let Some(idx) = self.is_valid_evidence(&evidence, alignment_properties) {
                    push_evidence(evidence, idx);
                }
            }
        }

        // METHOD: if all loci exceed the maximum depth, we subsample the evidence.
        // We cannot decide this per locus, because we risk adding more biases if loci have different alt allele sampling biases.
        let subsample = locus_depth.values().all(|depth| *depth > max_depth);
        let mut subsampler = sample::SubsampleCandidates::new(max_depth, candidates.len());
        let mut observations = Vec::new();
        for evidence in &candidates {
            if !subsample || subsampler.keep() {
                if let Some(obs) = self.evidence_to_observation(
                    evidence,
                    alignment_properties,
                    &homopolymer_error_model,
                    alt_variants,
                    observation_id_factory,
                )? {
                    observations.push(obs);
                }
            }
        }

        Ok(observations)
    }
}

pub(crate) trait Loci {
    fn contig(&self) -> Option<&str>;
    fn first_pos(&self) -> u64;
}

#[derive(Debug, Derefable, Builder, new, Clone)]
pub(crate) struct SingleLocus {
    #[deref]
    interval: genome::Interval,
}

impl AsRef<SingleLocus> for SingleLocus {
    fn as_ref(&self) -> &SingleLocus {
        self
    }
}

impl std::fmt::Display for SingleLocus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let start = self.interval.range().start + 1;
        let end = self.interval.range().end;
        if end > start {
            write!(f, "{}:{}-{}", self.interval.contig(), start, end)
        } else {
            write!(f, "{}:{}", self.interval.contig(), start)
        }
    }
}

impl SingleLocus {
    pub(crate) fn overlap(&self, record: &bam::Record, consider_clips: bool) -> Overlap {
        let mut pos = record.pos() as u64;
        let cigar = record.cigar_cached().unwrap();
        let mut end_pos = record.cigar_cached().unwrap().end_pos() as u64;

        if consider_clips {
            // consider soft clips for overlap detection
            pos = pos.saturating_sub(cigar.leading_softclips() as u64);
            end_pos += cigar.trailing_softclips() as u64;
        }

        if pos <= self.range().start {
            if end_pos >= self.range().end {
                return Overlap::Enclosing;
            } else if end_pos > self.range().start {
                return Overlap::Left;
            }
        } else if end_pos >= self.range().end && pos < self.range().end {
            return Overlap::Right;
        } else if pos >= self.range().start && end_pos <= self.range().end {
            return Overlap::Enclosed;
        }

        Overlap::None
    }

    fn outside_overlap(&self, record: &bam::Record) -> bool {
        let reverse_read = Self::read_reverse_strand(record.inner.core.flag);
        let pos = record.pos() as u64;
        if pos == self.range().start + 1 && reverse_read {
            return true;
        }
        false
    }

    /// Finds out whether the given string is a forward or reverse string.
    ///
    /// # Returns
    ///
    /// True if read given read is a reverse read, false if it is a forward read
    pub fn read_reverse_strand(flag: u16) -> bool {
        let read_paired = 0b1;
        let read_mapped_porper_pair = 0b01;
        let read_reverse = 0b10000;
        let mate_reverse = 0b100000;
        let first_in_pair = 0b1000000;
        let second_in_pair = 0b10000000;
        if (flag & read_paired != 0
            && flag & read_mapped_porper_pair != 0
            && flag & read_reverse != 0
            && flag & first_in_pair != 0)
            || (flag & read_paired != 0
                && flag & read_mapped_porper_pair != 0
                && flag & mate_reverse != 0
                && flag & second_in_pair != 0)
            || (flag & read_reverse != 0
                && (flag & read_paired == 0 || flag & read_mapped_porper_pair == 0))
        {
            return true;
        }
        false
    }
}

impl Loci for SingleLocus {
    fn first_pos(&self) -> u64 {
        self.range().start
    }
    fn contig(&self) -> Option<&str> {
        Some(self.interval.contig())
    }
}

#[derive(new, Default, Debug, Derefable, Clone)]
pub(crate) struct MultiLocus {
    #[deref(mutable)]
    loci: Vec<SingleLocus>,
}

impl MultiLocus {
    pub(crate) fn from_single_locus(locus: SingleLocus) -> Self {
        MultiLocus::new(vec![locus])
    }
}

impl Loci for MultiLocus {
    fn first_pos(&self) -> u64 {
        self[0].first_pos()
    }

    fn contig(&self) -> Option<&str> {
        let contig = self.loci[0].interval.contig();
        let is_single_contig = self.loci[1..]
            .iter()
            .all(|locus| locus.interval.contig() == contig);
        if is_single_contig {
            Some(contig)
        } else {
            None
        }
    }
}

#[derive(Debug)]
struct Candidate {
    left: Rc<bam::Record>,
    right: Option<Rc<bam::Record>>,
}

impl Candidate {
    fn new(record: Rc<bam::Record>) -> Self {
        Candidate {
            left: record,
            right: None,
        }
    }
}

/// Describes whether read overlaps a variant in a valid or invalid (too large overlap) way.
#[derive(Debug, PartialEq, Eq)]
pub(crate) enum Overlap {
    Enclosing,
    Left,
    Right,
    Enclosed,
    None,
}

impl Overlap {
    pub(crate) fn is_none(&self) -> bool {
        matches!(self, Overlap::None)
    }
}
