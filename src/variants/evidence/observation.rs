// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::hash::{Hash, Hasher};
use std::ops::Deref;
use std::rc::Rc;

use anyhow::Result;
use bio::stats::{LogProb, Prob};
use rust_htslib::bam;
use serde::ser::{SerializeStruct, Serializer};
use serde::Serialize;
// use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::bayesian::BayesFactor;
use itertools::Itertools;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::sample;
use crate::variants::types::Variant;

/// Calculate expected value of sequencing depth, considering mapping quality.
pub(crate) fn expected_depth(obs: &[Observation]) -> u32 {
    LogProb::ln_sum_exp(&obs.iter().map(|o| o.prob_mapping).collect_vec())
        .exp()
        .round() as u32
}

/// An observation for or against a variant.
#[derive(Clone, Debug, Builder, Default)]
pub(crate) struct Observation {
    /// Posterior probability that the read/read-pair has been mapped correctly (1 - MAPQ).
    #[builder(private)]
    pub(crate) prob_mapping: LogProb,
    /// Posterior probability that the read/read-pair has been mapped incorrectly (MAPQ).
    #[builder(private)]
    pub(crate) prob_mismapping: LogProb,
    /// Probability that the read/read-pair comes from the alternative allele.
    pub(crate) prob_alt: LogProb,
    /// Probability that the read/read-pair comes from the reference allele.
    pub(crate) prob_ref: LogProb,
    /// Probability that the read/read-pair comes from an unknown allele at an unknown true
    /// locus (in case it is mismapped). This should usually be set as the product of the maxima
    /// of prob_ref and prob_alt per read.
    pub(crate) prob_missed_allele: LogProb,
    /// Probability to sample the alt allele
    pub(crate) prob_sample_alt: LogProb,
    /// Probability to overlap with both strands
    #[builder(private)]
    pub(crate) prob_double_overlap: LogProb,
    /// Probability to overlap with one strand only (1-prob_double_overlap)
    #[builder(private)]
    pub(crate) prob_single_overlap: LogProb,
    /// Probability to observe any overlapping strand combination (when not associated with alt allele)
    pub(crate) prob_any_strand: LogProb,
    /// Observation relies on forward strand evidence
    pub(crate) forward_strand: bool,
    /// Observation relies on reverse strand evidence
    pub(crate) reverse_strand: bool,
}

impl ObservationBuilder {
    pub(crate) fn prob_mapping_mismapping(&mut self, prob_mapping: LogProb) -> &mut Self {
        self.prob_mapping(prob_mapping)
            .prob_mismapping(prob_mapping.ln_one_minus_exp())
    }

    pub(crate) fn prob_overlap(&mut self, prob_double_overlap: LogProb) -> &mut Self {
        self.prob_double_overlap(prob_double_overlap)
            .prob_single_overlap(prob_double_overlap.ln_one_minus_exp())
    }
}

impl Observation {
    pub(crate) fn bayes_factor_alt(&self) -> BayesFactor {
        BayesFactor::new(self.prob_alt, self.prob_ref)
    }
}

impl Serialize for Observation {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut s = serializer.serialize_struct("Observation", 3)?;
        s.serialize_field("prob_mapping", &self.prob_mapping)?;
        s.serialize_field("prob_mismapping", &self.prob_mismapping)?;
        s.serialize_field("prob_alt", &self.prob_alt)?;
        s.serialize_field("prob_ref", &self.prob_ref)?;
        s.serialize_field("prob_sample_alt", &self.prob_sample_alt)?;
        s.end()
    }
}

/// Something that can be converted into observations.
pub(crate) trait Observable<E>: Variant<Evidence = E>
where
    E: Evidence + Eq + Hash,
{
    fn extract_observations(
        &self,
        buffer: &mut sample::RecordBuffer,
        alignment_properties: &mut AlignmentProperties,
        max_depth: usize,
    ) -> Result<Vec<Observation>>;

    /// Convert MAPQ (from read mapper) to LogProb for the event that the read maps
    /// correctly.
    fn prob_mapping(&self, evidence: &E) -> LogProb;

    /// Calculate an observation from the given evidence.
    fn evidence_to_observation(
        &self,
        evidence: &E,
        alignment_properties: &AlignmentProperties,
    ) -> Result<Option<Observation>> {
        Ok(match self.allele_support(evidence, alignment_properties)? {
            // METHOD: only consider allele support if it comes either from forward or reverse strand.
            // Unstranded observations (e.g. only insert size), are too unreliable, or do not contain
            // any information (e.g. no overlap).
            Some(allele_support)
                if allele_support.forward_strand() || allele_support.reverse_strand() =>
            {
                Some(
                    ObservationBuilder::default()
                        .prob_mapping_mismapping(self.prob_mapping(evidence))
                        .prob_alt(allele_support.prob_alt_allele())
                        .prob_ref(allele_support.prob_ref_allele())
                        .prob_sample_alt(self.prob_sample_alt(evidence, alignment_properties))
                        .prob_missed_allele(allele_support.prob_missed_allele())
                        .prob_overlap(LogProb::ln_zero()) // no double overlap possible
                        .prob_any_strand(LogProb::from(Prob(0.5)))
                        .forward_strand(allele_support.forward_strand())
                        .reverse_strand(allele_support.reverse_strand())
                        .build()
                        .unwrap(),
                )
            }
            _ => None,
        })
    }
}

pub(crate) trait Evidence {}

#[derive(new, Clone, Eq, Debug)]
pub(crate) struct SingleEndEvidence {
    inner: Rc<bam::Record>,
}

impl Deref for SingleEndEvidence {
    type Target = bam::Record;

    fn deref(&self) -> &bam::Record {
        self.inner.as_ref()
    }
}

impl Evidence for SingleEndEvidence {}

impl PartialEq for SingleEndEvidence {
    fn eq(&self, other: &Self) -> bool {
        self.qname() == other.qname()
    }
}

impl Hash for SingleEndEvidence {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.qname().hash(state);
    }
}

#[derive(Clone, Eq, Debug)]
pub(crate) enum PairedEndEvidence {
    SingleEnd(Rc<bam::Record>),
    PairedEnd {
        left: Rc<bam::Record>,
        right: Rc<bam::Record>,
    },
}

impl PairedEndEvidence {
    pub(crate) fn qname(&self) -> &[u8] {
        match self {
            PairedEndEvidence::SingleEnd(rec) => rec.qname(),
            PairedEndEvidence::PairedEnd { left, .. } => left.qname(),
        }
    }
}

impl Evidence for PairedEndEvidence {}

impl PartialEq for PairedEndEvidence {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (PairedEndEvidence::SingleEnd(a), PairedEndEvidence::SingleEnd(b)) => {
                a.qname() == b.qname()
            }
            (
                PairedEndEvidence::PairedEnd { left: a, .. },
                PairedEndEvidence::PairedEnd { left: b, .. },
            ) => a.qname() == b.qname(),
            _ => false,
        }
    }
}

impl Hash for PairedEndEvidence {
    fn hash<H: Hasher>(&self, state: &mut H) {
        match self {
            PairedEndEvidence::SingleEnd(a) => a.qname().hash(state),
            PairedEndEvidence::PairedEnd { left: a, .. } => a.qname().hash(state),
        }
    }
}
