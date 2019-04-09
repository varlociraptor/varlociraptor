// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::f64;
use std::str;

use derive_builder::Builder;
use rgsl::randist::poisson::poisson_pdf;
use serde::ser::{SerializeStruct, Serializer};
use serde::Serialize;

use bio::stats::LogProb;
// use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::bayesian::BayesFactor;
use rust_htslib::bam;
use rust_htslib::bam::record::CigarString;
use itertools::Itertools;


/// Calculate expected value of sequencing depth, considering mapping quality.
pub fn expected_depth(obs: &[Observation]) -> u32 {
    LogProb::ln_sum_exp(&obs.iter().map(|o| o.prob_mapping).collect_vec()).exp().round() as u32
}

/// An observation for or against a variant.
#[derive(Clone, Debug, Builder)]
pub struct Observation {
    /// Posterior probability that the read/read-pair has been mapped correctly (1 - MAPQ).
    pub prob_mapping: LogProb,
    /// Posterior probability that the read/read-pair has been mapped incorrectly (MAPQ).
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
    /// Observation relies on forward strand evidence
    pub forward_strand: bool,
    /// Observation relies on reverse strand evidence
    pub reverse_strand: bool,
    /// Type of evidence.
    pub evidence: Evidence,
}

impl Observation {
    pub fn new(
        prob_mapping: LogProb,
        prob_mismapping: LogProb,
        prob_alt: LogProb,
        prob_ref: LogProb,
        prob_missed_allele: LogProb,
        prob_sample_alt: LogProb,
        forward_strand: bool,
        reverse_strand: bool,
        evidence: Evidence,
    ) -> Self {
        Observation {
            prob_mapping: prob_mapping,
            prob_mismapping: prob_mismapping,
            prob_alt: prob_alt,
            prob_ref: prob_ref,
            prob_missed_allele: prob_missed_allele,
            prob_sample_alt: prob_sample_alt,
            forward_strand: forward_strand,
            reverse_strand: reverse_strand,
            evidence: evidence,
        }
    }

    pub fn is_alignment_evidence(&self) -> bool {
        if let Evidence::Alignment(_) = self.evidence {
            true
        } else {
            false
        }
    }

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
        s.serialize_field("evidence", &self.evidence)?;
        s.end()
    }
}

/// Types of evidence that lead to an observation.
/// The contained information is intended for debugging and will be printed together with
/// observations.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum Evidence {
    /// Insert size of fragment
    InsertSize(String),
    /// Alignment of a single read
    Alignment(String),
}

impl Evidence {
    /// Create a dummy alignment.
    pub fn dummy_alignment() -> Self {
        Evidence::Alignment("Dummy-Alignment".to_owned())
    }

    /// Create dummy insert size evidence.
    pub fn dummy_insert_size(insert_size: u32) -> Self {
        Evidence::InsertSize(format!("insert-size={}", insert_size))
    }

    /// Create insert size evidence.
    pub fn insert_size(
        insert_size: u32,
        left: &CigarString,
        right: &CigarString,
        left_record: &bam::Record,
        right_record: &bam::Record,
        p_left_ref: LogProb,
        p_left_alt: LogProb,
        p_right_ref: LogProb,
        p_right_alt: LogProb,
        p_isize_ref: LogProb,
        p_isize_alt: LogProb,
    ) -> Self {
        Evidence::InsertSize(format!(
            "left: cigar={} ({:e} vs {:e}), right: cigar={} ({:e} vs {:e}), insert-size={} ({:e} vs {:e}), qname={}, left: AS={:?}, XS={:?}, right: AS={:?}, XS={:?}",
            left, p_left_ref.exp(), p_left_alt.exp(),
            right, p_right_ref.exp(), p_right_alt.exp(),
            insert_size, p_isize_ref.exp(), p_isize_alt.exp(),
            str::from_utf8(left_record.qname()).unwrap(),
            left_record.aux(b"AS").map(|a| a.integer()),
            left_record.aux(b"XS").map(|a| a.integer()),
            right_record.aux(b"AS").map(|a| a.integer()),
            right_record.aux(b"XS").map(|a| a.integer())
        ))
    }

    /// Create alignment evidence.
    pub fn alignment(cigar: &CigarString, record: &bam::Record) -> Self {
        Evidence::Alignment(format!(
            "cigar={}, qname={}, AS={:?}, XS={:?}",
            cigar,
            str::from_utf8(record.qname()).unwrap(),
            record.aux(b"AS").map(|a| a.integer()),
            record.aux(b"XS").map(|a| a.integer())
        ))
    }
}
