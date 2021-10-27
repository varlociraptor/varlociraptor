use std::cmp;
use std::collections::HashMap;
use std::hash::Hash;
use std::rc::Rc;

use anyhow::Result;
use bio::stats::probs::LogProb;
use bio::stats::Prob;
use bio_types::genome;
use ordered_float::NotNan;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::reference;
use crate::utils::homopolymers::homopolymer_indel_len;
use crate::variants::evidence::observation::{Observation, ReadPosition};
use crate::variants::model::bias::Bias;
use crate::variants::model::{Variant, VariantType};

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub(crate) enum HomopolymerError {
    None,
    Some,
}

impl HomopolymerError {
    pub(crate) fn values() -> Vec<Self> {
        vec![HomopolymerError::None, HomopolymerError::Some]
    }
}

impl Default for HomopolymerError {
    fn default() -> Self {
        HomopolymerError::None
    }
}

impl Bias for HomopolymerError {
    fn prob(&self, observation: &Observation<ReadPosition>) -> LogProb {
        match self {
            HomopolymerError::None => observation
                .prob_wildtype_homopolymer_error
                .unwrap_or(LogProb::ln_one()),
            HomopolymerError::Some => observation
                .prob_artifact_homopolymer_error
                .unwrap_or(LogProb::ln_zero()),
        }
    }

    fn prob_any(&self, _observation: &Observation<ReadPosition>) -> LogProb {
        LogProb::ln_one() // TODO check this
    }

    fn is_artifact(&self) -> bool {
        *self != HomopolymerError::None
    }

    fn is_informative(&self, pileups: &[Vec<Observation<ReadPosition>>]) -> bool {
        if !self.is_artifact() {
            return true;
        }
        // METHOD: this bias is only relevant if there is at least one recorded indel operation (indel operations are only recorded for some variants).
        pileups
            .iter()
            .any(|pileup| pileup.iter().any(|obs| self.is_bias_evidence(obs)))
    }

    fn is_possible(&self, pileups: &[Vec<Observation<ReadPosition>>]) -> bool {
        self.is_informative(pileups)
    }

    fn is_bias_evidence(&self, observation: &Observation<ReadPosition>) -> bool {
        observation
            .prob_artifact_homopolymer_error
            .map(|prob| prob != LogProb::ln_zero())
            .unwrap_or(false)
    }

    fn min_strong_evidence_ratio(&self) -> f64 {
        0.01 // TODO possibly tune, this is very permissive
    }
}
