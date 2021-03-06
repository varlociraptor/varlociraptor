// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::BTreeMap;
use std::rc::Rc;

use anyhow::Result;
use bio::stats::{LogProb, PHREDProb};
use bio_types::genome::{self, AbstractInterval};
use rust_htslib::bam;
use vec_map::VecMap;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::evidence::observation::{
    Evidence, Observable, Observation, PairedEndEvidence, SingleEndEvidence, Strand,
};
use crate::variants::sample;

pub(crate) mod breakends;
pub(crate) mod deletion;
pub(crate) mod duplication;
pub(crate) mod insertion;
pub(crate) mod inversion;
pub(crate) mod mnv;
pub(crate) mod none;
pub(crate) mod replacement;
pub(crate) mod snv;

pub(crate) use deletion::Deletion;
pub(crate) use duplication::Duplication;
pub(crate) use insertion::Insertion;
pub(crate) use inversion::Inversion;
pub(crate) use mnv::MNV;
pub(crate) use none::None;
pub(crate) use replacement::Replacement;
pub(crate) use snv::SNV;

#[derive(Debug, CopyGetters, Builder)]
#[getset(get_copy = "pub")]
pub(crate) struct AlleleSupport {
    prob_ref_allele: LogProb,
    prob_alt_allele: LogProb,
    strand: Strand,
    #[builder(default)]
    read_position: Option<u32>,
}

impl AlleleSupport {
    /// METHOD: This is an estimate of the allele likelihood at the true location in case
    /// the read is mismapped. The value has to be approximately in the range of prob_alt
    /// and prob_ref. Otherwise it could cause numerical problems, by dominating the
    /// likelihood such that subtle differences in allele frequencies become numercically
    /// invisible in the resulting likelihood.
    pub(crate) fn prob_missed_allele(&self) -> LogProb {
        self.prob_ref_allele.ln_add_exp(self.prob_alt_allele) - LogProb(2.0_f64.ln())
    }

    pub(crate) fn merge(&mut self, other: &AlleleSupport) -> &mut Self {
        self.prob_ref_allele += other.prob_ref_allele;
        self.prob_alt_allele += other.prob_alt_allele;
        if self.strand == Strand::None {
            self.strand = other.strand;
        } else if other.strand != Strand::None && self.strand != other.strand {
            self.strand = Strand::Both;
        }

        self
    }
}

pub(crate) trait Variant {
    type Evidence: Evidence;
    type Loci: Loci;

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
    ) -> Result<Option<AlleleSupport>>;

    /// Calculate probability to sample a record length like the given one from the alt allele.
    fn prob_sample_alt(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb;
}

impl<V> Observable<SingleEndEvidence> for V
where
    V: Variant<Evidence = SingleEndEvidence, Loci = SingleLocus>,
{
    fn prob_mapping(&self, evidence: &SingleEndEvidence) -> LogProb {
        let prob_mismapping = LogProb::from(PHREDProb(evidence.mapq() as f64));
        prob_mismapping.ln_one_minus_exp()
    }

    fn extract_observations(
        &self,
        buffer: &mut sample::RecordBuffer,
        alignment_properties: &mut AlignmentProperties,
        max_depth: usize,
    ) -> Result<Vec<Observation>> {
        let locus = self.loci();
        buffer.fetch(locus, false)?;

        let candidates: Vec<_> = buffer
            .iter()
            .filter_map(|record| {
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
                if let Some(obs) = self.evidence_to_observation(&evidence, alignment_properties)? {
                    observations.push(obs);
                }
            }
        }
        Ok(observations)
    }
}

impl<V> Observable<PairedEndEvidence> for V
where
    V: Variant<Evidence = PairedEndEvidence, Loci = MultiLocus>,
{
    fn prob_mapping(&self, evidence: &PairedEndEvidence) -> LogProb {
        let prob = |record: &bam::Record| LogProb::from(PHREDProb(record.mapq() as f64));
        match evidence {
            PairedEndEvidence::SingleEnd(record) => prob(record).ln_one_minus_exp(),
            PairedEndEvidence::PairedEnd { left, right } => {
                // METHOD: take maximum of the (log-spaced) mapping quality of the left and the right read.
                // In BWA, MAPQ is influenced by the mate, hence they are not independent
                // and we can therefore not multiply them. By taking the maximum, we
                // make a conservative choice (since 1-mapq is the mapping probability).
                let mut p = prob(left);
                let mut q = prob(right);
                if p < q {
                    std::mem::swap(&mut p, &mut q);
                }
                p.ln_one_minus_exp()
            }
        }
    }

    fn extract_observations(
        &self,
        buffer: &mut sample::RecordBuffer,
        alignment_properties: &mut AlignmentProperties,
        max_depth: usize,
    ) -> Result<Vec<Observation>> {
        // We cannot use a hash function here because candidates have to be considered
        // in a deterministic order. Otherwise, subsampling high-depth regions will result
        // in slightly different probabilities each time.
        let mut candidate_records = BTreeMap::new();

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
                    left: Rc::clone(&candidate.left),
                    right: Rc::clone(right),
                };
                if let Some(idx) = self.is_valid_evidence(&evidence, alignment_properties) {
                    push_evidence(evidence, idx);
                }
            } else {
                // this is a single alignment with unmapped mate or mate outside of the
                // region of interest
                let evidence = PairedEndEvidence::SingleEnd(Rc::clone(&candidate.left));
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
                if let Some(obs) = self.evidence_to_observation(evidence, alignment_properties)? {
                    observations.push(obs);
                }
            }
        }

        Ok(observations)
    }
}

pub(crate) trait Loci {}

#[derive(Debug, Derefable, Builder, new, Clone)]
pub(crate) struct SingleLocus {
    #[deref]
    interval: genome::Interval,
    #[builder(default = "true")]
    #[new(value = "true")]
    from_left: bool,
    #[builder(default = "true")]
    #[new(value = "true")]
    from_right: bool,
}

impl AsRef<SingleLocus> for SingleLocus {
    fn as_ref(&self) -> &SingleLocus {
        self
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
}

impl Loci for SingleLocus {}

#[derive(new, Default, Debug, Derefable, Clone)]
pub(crate) struct MultiLocus {
    #[deref(mutable)]
    loci: Vec<SingleLocus>,
}

impl Loci for MultiLocus {}

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
#[derive(Debug)]
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
