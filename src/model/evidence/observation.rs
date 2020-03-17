// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::f64;

use anyhow::Result;
use rgsl::randist::poisson::poisson_pdf;
use serde::ser::{SerializeStruct, Serializer};
use serde::Serialize;

use bio::stats::LogProb;
// use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::bayesian::BayesFactor;
use itertools::Itertools;

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
