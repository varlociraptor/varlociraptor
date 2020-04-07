// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::f64;
use std::ops::Deref;
use std::collections::HashSet;
use std::hash::{Hash, Hasher};

use anyhow::Result;
use bio::stats::{LogProb, Prob};
use rgsl::randist::poisson::poisson_pdf;
use rust_htslib::bam;
use serde::ser::{SerializeStruct, Serializer};
use serde::Serialize;
// use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::bayesian::BayesFactor;
use itertools::Itertools;
use vec_map::VecMap;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::model::sample;
use crate::variants::Variant;

/// Calculate expected value of sequencing depth, considering mapping quality.
pub fn expected_depth(obs: &[Observation]) -> u32 {
    LogProb::ln_sum_exp(&obs.iter().map(|o| o.prob_mapping).collect_vec())
        .exp()
        .round() as u32
}

/// An observation for or against a variant.
#[derive(Clone, Debug, Builder, Default)]
pub struct Observation {
    /// Posterior probability that the read/read-pair has been mapped correctly (1 - MAPQ).
    #[builder(private)]
    pub prob_mapping: LogProb,
    /// Posterior probability that the read/read-pair has been mapped incorrectly (MAPQ).
    #[builder(private)]
    pub prob_mismapping: LogProb,
    /// Probability that the read/read-pair comes from the alternative allele.
    pub prob_alt: LogProb,
    /// Probability that the read/read-pair comes from the reference allele.
    pub prob_ref: LogProb,
    /// Probability that the read/read-pair comes from an unknown allele at an unknown true
    /// locus (in case it is mismapped). This should usually be set as the product of the maxima
    /// of prob_ref and prob_alt per read.
    pub prob_missed_allele: LogProb,
    /// Probability to sample the alt allele
    pub prob_sample_alt: LogProb,
    /// Probability to overlap with both strands
    #[builder(private)]
    pub prob_double_overlap: LogProb,
    /// Probability to overlap with one strand only (1-prob_double_overlap)
    #[builder(private)]
    pub prob_single_overlap: LogProb,
    /// Probability to observe any overlapping strand combination (when not associated with alt allele)
    pub prob_any_strand: LogProb,
    /// Observation relies on forward strand evidence
    pub forward_strand: bool,
    /// Observation relies on reverse strand evidence
    pub reverse_strand: bool,
}

impl ObservationBuilder {
    pub fn prob_mapping_mismapping(&mut self, prob_mapping: LogProb) -> &mut Self {
        self.prob_mapping(prob_mapping)
            .prob_mismapping(prob_mapping.ln_one_minus_exp())
    }

    pub fn prob_overlap(&mut self, prob_double_overlap: LogProb) -> &mut Self {
        self.prob_double_overlap(prob_double_overlap)
            .prob_single_overlap(prob_double_overlap.ln_one_minus_exp())
    }
}

impl Observation {
    pub fn bayes_factor_alt(&self) -> BayesFactor {
        BayesFactor::new(self.prob_alt, self.prob_ref)
    }

    pub fn bayes_factor_ref(&self) -> BayesFactor {
        BayesFactor::new(self.prob_ref, self.prob_alt)
    }

    #[inline]
    pub fn prob_alt_forward(&self) -> LogProb {
        if self.forward_strand {
            self.prob_alt
        } else {
            LogProb::ln_zero()
        }
    }

    #[inline]
    pub fn prob_alt_reverse(&self) -> LogProb {
        if self.reverse_strand {
            self.prob_alt
        } else {
            LogProb::ln_zero()
        }
    }

    #[inline]
    pub fn prob_ref_forward(&self) -> LogProb {
        if self.forward_strand {
            self.prob_ref
        } else {
            LogProb::ln_zero()
        }
    }

    #[inline]
    pub fn prob_ref_reverse(&self) -> LogProb {
        if self.reverse_strand {
            self.prob_ref
        } else {
            LogProb::ln_zero()
        }
    }
}

pub fn poisson_pmf(count: u32, mu: f64) -> LogProb {
    if mu == 0.0 {
        if count == 0 {
            LogProb::ln_one()
        } else {
            LogProb::ln_zero()
        }
    } else {
        LogProb(poisson_pdf(count, mu).ln())
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
pub trait Observable<'a, E>: Variant<'a, Evidence=E>
where
    E: Evidence<'a> + Eq + Hash,
{
    fn extract_observations(
        &self,
        buffer: &'a mut sample::RecordBuffer,
        alignment_properties: &mut AlignmentProperties,
        max_depth: usize,
    ) -> Result<Vec<Observation>>;

    /// Convert MAPQ (from read mapper) to LogProb for the event that the read maps
    /// correctly.
    fn prob_mapping(&self, evidence: &E) -> LogProb;

    /// Calculate strand information.
    fn strand(&self, evidence: &E) -> Strand;

    /// Calculate an observation from the given evidence.
    fn evidence_to_observation(&self, evidence: &E, alignment_properties: &AlignmentProperties) -> Result<Option<Observation>> {
        Ok(match self.prob_alleles(evidence)? {
            Some(prob_alleles) => {
                let strand = self.strand(evidence);

                Some(ObservationBuilder::default()
                    .prob_mapping_mismapping(self.prob_mapping(evidence))
                    .prob_alt(prob_alleles.alt_allele())
                    .prob_ref(prob_alleles.ref_allele())
                    .prob_sample_alt(
                        self.prob_sample_alt(evidence, alignment_properties),
                    )
                    .prob_missed_allele(prob_alleles.missed_allele())
                    .prob_overlap(LogProb::ln_zero()) // no double overlap possible
                    .prob_any_strand(LogProb::from(Prob(0.5)))
                    .forward_strand(strand.forward())
                    .reverse_strand(strand.reverse())
                    .build()
                    .unwrap())
            }
            None => None,
        })
    }
}

#[derive(Builder, Default, CopyGetters)]
#[getset(get_copy = "pub")]
pub struct Strand {
    forward: bool,
    reverse: bool,
}

pub trait Evidence<'a> {}

#[derive(new, Clone, Copy, Eq)]
pub struct SingleEndEvidence<'a> {
    inner: &'a bam::Record,
}

impl<'a> Deref for SingleEndEvidence<'a> {
    type Target = bam::Record;

    fn deref(&self) -> &bam::Record {
        self.inner
    }
}

impl<'a> Evidence<'a> for SingleEndEvidence<'a> {}

impl<'a> PartialEq for SingleEndEvidence<'a> {
    fn eq(&self, other: &Self) -> bool {
        self.qname() == other.qname()
    }
}

impl<'a> Hash for SingleEndEvidence<'a> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.qname().hash(state);
    }
}

#[derive(Clone, Copy, Eq)]
pub enum PairedEndEvidence<'a> {
    SingleEnd(&'a bam::Record),
    PairedEnd {
        left: &'a bam::Record,
        right: &'a bam::Record,
    },
}

impl<'a> Evidence<'a> for PairedEndEvidence<'a> {}

impl<'a> PartialEq for PairedEndEvidence<'a> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (PairedEndEvidence::SingleEnd(a), PairedEndEvidence::SingleEnd(b)) => a.qname() == b.qname(),
            (PairedEndEvidence::PairedEnd { left: a, .. }, PairedEndEvidence::PairedEnd { left: b, .. }) => a.qname() == b.qname(),
            _ => false,
        }
    }
}


impl<'a> Hash for PairedEndEvidence<'a> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        match self {
            PairedEndEvidence::SingleEnd(a) => a.qname().hash(state),
            PairedEndEvidence::PairedEnd { left: a, .. } => a.qname().hash(state),
        }
    }
}