use std::collections::BTreeMap;

use anyhow::Result;
use bio::stats::{LogProb, PHREDProb, Prob};
use rust_htslib::bam;
use bio_types::{genome, genome::AbstractInterval};

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::model::evidence::observation::{
    Evidence, Observable, Observation, ObservationBuilder, PairedEndEvidence, SingleEndEvidence,
    Strand, StrandBuilder,
};
use crate::model::evidence;
use crate::model::sample;
use crate::utils::Overlap;

pub mod mnv;
pub mod snv;

#[derive(Debug, CopyGetters, new)]
#[getset(get_copy = "pub")]
pub struct AlleleProb {
    ref_allele: LogProb,
    alt_allele: LogProb,
}

impl AlleleProb {
    pub fn missed_allele(&self) -> LogProb {
        self.ref_allele.ln_add_exp(self.alt_allele) - LogProb(2.0_f64.ln())
    }
}

pub trait Variant<'a> {
    type Evidence: Evidence<'a>;
    type Loci: Loci;

    /// Determine whether the evidence is suitable to assessing probabilities 
    /// (i.e. overlaps the variant in the right way).
    fn is_valid_evidence(&self, evidence: &Self::Evidence) -> bool;

    /// Return variant loci.
    fn loci(&self) -> &Self::Loci;

    /// Calculate probability for alt and reference allele.
    fn prob_alleles(&self, evidence: &Self::Evidence) -> Result<Option<AlleleProb>>;

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
        // TODO smaller window, constrained by read length only!
        buffer.fetch(locus)?;

        let candidates: Vec<_> = buffer
            .iter()
            .filter_map(|record| {
                let evidence = SingleEndEvidence::new(record);
                if self.is_valid_evidence(&evidence) {
                    Some(evidence)
                } else {
                    None
                }
            })
            .collect();

        self.evidence_to_observations(&candidates, alignment_properties, max_depth)
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

        // TODO centerpoint locus for deletions (in trait impl)
        let mut fetches = buffer.build_fetches();
        for locus in self.loci().iter() {
            fetches.push(locus);
        }
        for interval in fetches.iter() {
            buffer.fetch(interval)?;

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
                        && (candidate.left.is_last_in_template()
                            && record.is_last_in_template())
                    {
                        // ignore another partial alignment right of the first
                        continue;
                    }
                    // replace right record (we seek for the rightmost (partial) alignment)
                    candidate.right = Some(record);
                }
            }
        }

        let candidates: Vec<_> = candidate_records.values().filter_map(|candidate| {
            if let Some(right) = candidate.right {
                if candidate.left.mapq() == 0 || right.mapq() == 0 {
                    // Ignore pairs with ambiguous alignments.
                    // The statistical model does not consider them anyway.
                    return None;
                }
                let evidence = PairedEndEvidence::PairedEnd { left: candidate.left, right: right };
                if self.is_valid_evidence(&evidence) {
                    Some(evidence)
                } else {
                    None
                }
            } else {
                // this is a single alignment with unmapped mate or mate outside of the
                // region of interest
                let evidence = PairedEndEvidence::SingleEnd(candidate.left);
                if self.is_valid_evidence(&evidence) {
                    Some(evidence)
                } else {
                    None
                }
            }
        }).collect();

        self.evidence_to_observations(&candidates, alignment_properties, max_depth)
    }
}

pub trait Loci {}

#[derive(Debug, Derefable, new)]
pub struct SingleLocus(#[deref] genome::Interval);

impl Loci for SingleLocus {}

#[derive(Default, Debug, Derefable)]
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