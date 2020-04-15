use std::collections::{BTreeMap, HashSet};

use anyhow::Result;
use bio::stats::{LogProb, PHREDProb};
use bio_types::genome::{self, AbstractInterval};
use rust_htslib::bam;
use vec_map::VecMap;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::model::evidence::observation::{
    Evidence, Observable, Observation, PairedEndEvidence, SingleEndEvidence, Strand, StrandBuilder,
};
use crate::model::sample;

pub mod deletion;
pub mod fragment_enclosable;
pub mod mnv;
pub mod realignable;
pub mod snv;

pub use fragment_enclosable::FragmentEnclosable;

#[derive(Debug, CopyGetters, new)]
#[getset(get_copy = "pub")]
pub struct AlleleProb {
    ref_allele: LogProb,
    alt_allele: LogProb,
}

impl AlleleProb {
    /// METHOD: This is an estimate of the allele likelihood at the true location in case
    /// the read is mismapped. The value has to be approximately in the range of prob_alt
    /// and prob_ref. Otherwise it could cause numerical problems, by dominating the
    /// likelihood such that subtle differences in allele frequencies become numercically
    /// invisible in the resulting likelihood.
    pub fn missed_allele(&self) -> LogProb {
        self.ref_allele.ln_add_exp(self.alt_allele) - LogProb(2.0_f64.ln())
    }
}

pub trait Variant<'a> {
    type Evidence: Evidence<'a>;
    type Loci: Loci;

    /// Determine whether the evidence is suitable to assessing probabilities
    /// (i.e. overlaps the variant in the right way).
    ///
    /// # Returns
    ///
    /// The index of the loci for which this evidence is valid, `None` if invalid.
    fn is_valid_evidence(&self, evidence: &Self::Evidence) -> Option<Vec<usize>>;

    /// Return variant loci.
    fn loci(&self) -> &Self::Loci;

    /// Calculate probability for alt and reference allele.
    fn prob_alleles(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> Result<Option<AlleleProb>>;

    /// Calculate probability to sample a record length like the given one from the alt allele.
    fn prob_sample_alt(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb;
}

impl<'a, V> Observable<'a, SingleEndEvidence<'a>> for V
where
    V: Variant<'a, Evidence = SingleEndEvidence<'a>, Loci = SingleLocus>,
{
    fn prob_mapping(&self, evidence: &SingleEndEvidence<'a>) -> LogProb {
        let prob_mismapping = LogProb::from(PHREDProb(evidence.mapq() as f64));
        let prob_mapping = prob_mismapping.ln_one_minus_exp();
        prob_mapping
    }

    fn strand(&self, evidence: &SingleEndEvidence<'a>) -> Strand {
        let reverse = evidence.flags() & 0x10 != 0;

        StrandBuilder::default()
            .reverse(reverse)
            .forward(!reverse)
            .build()
            .unwrap()
    }

    fn extract_observations(
        &self,
        buffer: &'a mut sample::RecordBuffer,
        alignment_properties: &mut AlignmentProperties,
        max_depth: usize,
    ) -> Result<Vec<Observation>> {
        let locus = self.loci();
        buffer.fetch(locus, false)?;

        let candidates: Vec<_> = buffer
            .iter()
            .filter_map(|record| {
                let evidence = SingleEndEvidence::new(record);
                if self.is_valid_evidence(&evidence).is_some() {
                    Some(evidence)
                } else {
                    None
                }
            })
            .collect();

        let subsampler = sample::SubsampleCandidates::new(max_depth, candidates.len());

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

impl<'a, V> Observable<'a, PairedEndEvidence<'a>> for V
where
    V: Variant<'a, Evidence = PairedEndEvidence<'a>, Loci = MultiLocus>,
{
    fn prob_mapping(&self, evidence: &PairedEndEvidence<'a>) -> LogProb {
        let prob = |record: &bam::Record| LogProb::from(PHREDProb(record.mapq() as f64));
        match evidence {
            PairedEndEvidence::SingleEnd(record) => prob(record),
            PairedEndEvidence::PairedEnd { left, right } => prob(left) + prob(right),
        }
    }

    fn strand(&self, evidence: &PairedEndEvidence<'a>) -> Strand {
        let is_reverse = |record: &bam::Record| record.flags() & 0x10 != 0;

        let (forward, reverse) = match evidence {
            PairedEndEvidence::SingleEnd(record) => {
                let reverse = is_reverse(record);
                (!reverse, reverse)
            }
            PairedEndEvidence::PairedEnd { left, right } => {
                let left_reverse = is_reverse(left);
                let right_reverse = is_reverse(right);

                let reverse = left_reverse || right_reverse;
                (!reverse, reverse)
            }
        };

        StrandBuilder::default()
            .reverse(reverse)
            .forward(forward)
            .build()
            .unwrap()
    }

    /// TODO move this into the prob_allele implementation
    /// Here we calculate the product of insert size based and alignment based probabilities.
    /// This has the benefit that the calculation automatically checks for consistence between
    /// insert size and overlapping alignmnments.
    /// This sports the following desirable behavior:
    ///
    /// * If there is no clear evidence from either the insert size or the alignment, the factors
    ///   simply scale because the probabilities of the corresponding type of evidence will be equal.
    /// * If reads and fragments agree, 1 stays 1 and 0 stays 0.
    /// * If reads and fragments disagree (the promising part!), the other type of evidence will
    ///   scale potential false positive probabilities down.
    /// * Since there is only one observation per fragment, there is no double counting when
    ///   estimating allele frequencies. Before, we had one observation for an overlapping read
    ///   and potentially another observation for the corresponding fragment.
    fn extract_observations(
        &self,
        buffer: &'a mut sample::RecordBuffer,
        alignment_properties: &mut AlignmentProperties,
        max_depth: usize,
    ) -> Result<Vec<Observation>> {
        // We cannot use a hash function here because candidates have to be considered
        // in a deterministic order. Otherwise, subsampling high-depth regions will result
        // in slightly different probabilities each time.
        let mut candidate_records = BTreeMap::new();

        // TODO centerpoint locus for deletions (in trait impl of self.loci())
        let mut fetches = buffer.build_fetches(true);
        for locus in self.loci().iter() {
            fetches.push(locus);
        }
        for interval in fetches.iter() {
            buffer.fetch(interval, true)?;

            for record in buffer.iter() {
                // METHOD: First, we check whether the record contains an indel in the cigar.
                // We store the maximum indel size to update the global estimates, in case
                // it is larger in this region.
                alignment_properties.update_max_cigar_ops_len(record);

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
                        // ignore another partial alignment right of the first
                        continue;
                    }
                    // replace right record (we seek for the rightmost (partial) alignment)
                    candidate.right = Some(record);
                }
            }
        }

        let mut candidates = VecMap::new();
        let push_evidence = |evidence, idx| {
            for i in idx {
                let entry = candidates.entry(i).or_insert_with(|| Vec::new());
                entry.push(evidence);
            }
        };

        for candidate in candidate_records.values() {
            if let Some(right) = candidate.right {
                if candidate.left.mapq() == 0 || right.mapq() == 0 {
                    // Ignore pairs with ambiguous alignments.
                    // The statistical model does not consider them anyway.
                    continue;
                }
                let evidence = PairedEndEvidence::PairedEnd {
                    left: candidate.left,
                    right: right,
                };
                if let Some(idx) = self.is_valid_evidence(&evidence) {
                    push_evidence(evidence, idx);
                }
            } else {
                // this is a single alignment with unmapped mate or mate outside of the
                // region of interest
                let evidence = PairedEndEvidence::SingleEnd(candidate.left);
                if let Some(idx) = self.is_valid_evidence(&evidence) {
                    push_evidence(evidence, idx);
                }
            }
        }

        // METHOD: if all loci exceed the maximum depth, we subsample the evidence.
        // We cannot decide this per locus, because we risk adding more biases if loci have different alt allele sampling biases.
        let subsample = candidates
            .values()
            .map(|locus_candidates| locus_candidates.len())
            .all(|l| l > max_depth);
        let candidates: HashSet<_> = candidates.values().flatten().collect();

        let subsampler = sample::SubsampleCandidates::new(max_depth, candidates.len());

        let mut observations = Vec::new();
        for evidence in candidates {
            if !subsample || subsampler.keep() {
                if let Some(obs) = self.evidence_to_observation(evidence, alignment_properties)? {
                    observations.push(obs);
                }
            }
        }

        Ok(observations)
    }
}

pub trait Loci {}

#[derive(Debug, Derefable, new)]
pub struct SingleLocus(#[deref] genome::Interval);

impl SingleLocus {
    pub fn overlap(&self, record: &bam::Record, consider_clips: bool) -> Overlap {
        let mut pos = record.pos() as u64;
        let cigar = record.cigar_cached().unwrap();
        let mut end_pos = record.cigar_cached().unwrap().end_pos() as u64;

        if consider_clips {
            // consider soft clips for overlap detection
            pos = pos.saturating_sub(cigar.leading_softclips() as u64);
            end_pos = end_pos + cigar.trailing_softclips() as u64;
        }

        if pos <= self.range().start {
            if end_pos > self.range().end {
                return Overlap::Enclosing;
            } else if end_pos > self.range().start {
                return Overlap::Left;
            }
        } else if end_pos > self.range().end && pos < self.range().end {
            return Overlap::Right;
        }

        Overlap::None
    }
}

impl Loci for SingleLocus {}

#[derive(new, Default, Debug, Derefable)]
pub struct MultiLocus {
    #[deref]
    loci: Vec<SingleLocus>,
}

impl Loci for MultiLocus {}

struct Candidate<'a> {
    left: &'a bam::Record,
    right: Option<&'a bam::Record>,
}

impl<'a> Candidate<'a> {
    fn new(record: &'a bam::Record) -> Self {
        Candidate {
            left: record,
            right: None,
        }
    }
}

/// Describes whether read overlaps a variant in a valid or invalid (too large overlap) way.
#[derive(Debug)]
pub enum Overlap {
    Enclosing,
    Left,
    Right,
    None,
}

impl Overlap {
    pub fn is_none(&self) -> bool {
        if let Overlap::None = self {
            true
        } else {
            false
        }
    }
}
