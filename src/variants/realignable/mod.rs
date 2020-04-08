use std::cmp;
use std::fmt::Debug;

use anyhow::Result;
use bio::stats::{LogProb, PHREDProb, Prob};

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::model::Variant;

pub mod edit_distance;
pub mod pairhmm;

pub trait Realignable {

}
