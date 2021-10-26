use std::cmp;
use std::collections::HashMap;
use std::rc::Rc;
use std::hash::Hash;

use bio::stats::Prob;
use bio::stats::probs::LogProb;
use bio_types::genome;
use ordered_float::NotNan;
use anyhow::Result;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::reference;
use crate::utils::homopolymers::homopolymer_indel_len;
use crate::variants::evidence::observation::{Observation, ReadPosition};
use crate::variants::model::{Variant, VariantType};
use crate::variants::model::bias::Bias;

pub(crate) type HomopolymerErrorModel = HashMap<i8, LogProb>;

#[derive(Clone, Debug)]
pub(crate) enum HomopolymerError {
    None { error_model: Option<HomopolymerErrorModel> },
    Some { error_model: Option<HomopolymerErrorModel> },
}

impl PartialEq for HomopolymerError {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::None { .. }, Self::None { .. }) => true,
            (Self::Some { .. }, Self::Some { .. }) => true,
            _ => false,
        }
    }
}

impl Eq for HomopolymerError {}

impl Hash for HomopolymerError {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        core::mem::discriminant(self).hash(state);
    }
}

impl HomopolymerError {
    pub(crate) fn values() -> Vec<Self> {
        vec![
            HomopolymerError::None { error_model: None },
            HomopolymerError::Some { error_model: None },
        ]
    }
}

impl Default for HomopolymerError {
    fn default() -> Self {
        HomopolymerError::None { error_model: None }
    }
}

impl Bias for HomopolymerError {
    fn prob(&self, observation: &Observation<ReadPosition>) -> LogProb {
        match self {
            HomopolymerError::None { error_model: Some(error_model) } => {
                error_model.get(&observation.homopolymer_indel_len.unwrap_or(0)).cloned().unwrap_or(LogProb::ln_zero())
            }
            HomopolymerError::Some { error_model: Some(error_model) } => {
                error_model.get(&observation.homopolymer_indel_len.unwrap_or(0)).cloned().unwrap_or(LogProb::ln_zero())
            }
            _ => {
                unreachable!("bug: no error model learned for homopolymer errors.");
            }
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
            .any(|pileup| pileup.iter().any(|obs| obs.homopolymer_indel_len))
    }

    fn is_possible(&self, pileups: &[Vec<Observation<ReadPosition>>]) -> bool {
        if !self.is_artifact() {
            return true;
        }

        pileups.iter().any(|pileup| {
            pileup.iter().any(|observation| match self {
                HomopolymerError::Some { .. } => observation.homopolymer_indel_len,
                HomopolymerError::None => self.prob(observation) != LogProb::ln_zero(),
            })
        })
    }

    fn is_bias_evidence(&self, observation: &Observation<ReadPosition>) -> bool {
        observation.homopolymer_indel_len
    }

    fn learn_parameters(&mut self, pileups: &[Vec<Observation<ReadPosition>>], alignment_properties: &AlignmentProperties, variant: &Variant, locus: &genome::Locus, reference_buffer: &reference::Buffer) -> Result<()> {
        match self {
            HomopolymerError::None { ref mut error_model } => {
                error_model.replace(alignment_properties.homopolymer_error_model.iter().map(|(error, prob)| (*error, LogProb::from(Prob(*prob)))).collect());
            }
            HomopolymerError::Some { ref mut error_model } => {
                let indel_len = homopolymer_indel_len(variant, locus, reference_buffer)?;
                if let Some(indel_len) = indel_len {
                    error_model.replace(alignment_properties.homopolymer_error_model.iter().map(|(error, prob)| {
                        (*error - indel_len, LogProb::from(Prob(*prob)))
                    }).collect());
                } else {
                    *error_model = None;
                }
                
            }
        }
        Ok(())
    }
}
