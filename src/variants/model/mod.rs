// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::{Deref, Range};
use std::str;

use ordered_float::NotNan;
use strum_macros::{EnumIter, EnumString, IntoStaticStr};

use crate::grammar;
use crate::variants::model::bias::Biases;

use super::sample::Pileup;

pub(crate) mod bias;
pub(crate) mod likelihood;
pub(crate) mod modes;
pub(crate) mod prior;

#[derive(Debug, Clone)]
pub(crate) struct Contamination {
    pub(crate) by: usize,
    pub(crate) fraction: f64,
}

#[derive(Eq, PartialEq, Clone, Debug, Hash)]
pub(crate) struct Event {
    pub(crate) name: String,
    pub(crate) kind: EventKind,
}

impl Event {
    pub(crate) fn learn_parameters(&mut self, pileups: &[Pileup]) -> () {
        if let Event {
            kind: EventKind::Biases(biases),
            ..
        } = self
        {
            for bias in biases {
                bias.learn_parameters(pileups);
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) enum EventKind {
    VAFTree(grammar::VAFTree),
    Biases(Vec<Biases>),
}

pub(crate) type AlleleFreq = NotNan<f64>;

#[allow(non_snake_case)]
pub(crate) fn AlleleFreq(af: f64) -> AlleleFreq {
    NotNan::new(af).unwrap()
}

pub(crate) trait AlleleFreqs: Debug {
    fn is_absent(&self) -> bool;
}
impl AlleleFreqs for DiscreteAlleleFreqs {
    fn is_absent(&self) -> bool {
        self.inner.len() == 1 && self.inner[0] == AlleleFreq(0.0)
    }
}
impl AlleleFreqs for ContinuousAlleleFreqs {
    fn is_absent(&self) -> bool {
        self.is_singleton() && self.start == AlleleFreq(0.0)
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct DiscreteAlleleFreqs {
    inner: Vec<AlleleFreq>,
}

impl Deref for DiscreteAlleleFreqs {
    type Target = Vec<AlleleFreq>;

    fn deref(&self) -> &Vec<AlleleFreq> {
        &self.inner
    }
}

/// An allele frequency range
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ContinuousAlleleFreqs {
    inner: Range<AlleleFreq>,
    pub(crate) left_exclusive: bool,
    pub(crate) right_exclusive: bool,
    /// offset to add when calculating the smallest observable value for a left-exclusive 0.0 bound
    zero_offset: NotNan<f64>,
}

impl ContinuousAlleleFreqs {
    pub(crate) fn absent() -> Self {
        Self::singleton(0.0)
    }

    pub(crate) fn singleton(value: f64) -> Self {
        ContinuousAlleleFreqs {
            inner: AlleleFreq(value)..AlleleFreq(value),
            left_exclusive: false,
            right_exclusive: false,
            zero_offset: NotNan::from(1.0),
        }
    }

    pub(crate) fn is_singleton(&self) -> bool {
        self.start == self.end
    }
}

impl Default for ContinuousAlleleFreqs {
    fn default() -> Self {
        Self::absent()
    }
}

impl Ord for ContinuousAlleleFreqs {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.inner.start.cmp(&other.start) {
            Ordering::Equal => self.inner.end.cmp(&other.end),
            ord => ord,
        }
    }
}

impl PartialOrd for ContinuousAlleleFreqs {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Deref for ContinuousAlleleFreqs {
    type Target = Range<AlleleFreq>;

    fn deref(&self) -> &Range<AlleleFreq> {
        &self.inner
    }
}

#[derive(
    Display,
    Debug,
    Clone,
    Serialize,
    Deserialize,
    EnumString,
    EnumIter,
    IntoStaticStr,
    EnumVariantNames,
)]
pub enum VariantType {
    #[strum(serialize = "INS")]
    Insertion(Option<Range<u64>>),
    #[strum(serialize = "DEL")]
    Deletion(Option<Range<u64>>),
    #[strum(serialize = "SNV")]
    Snv,
    #[strum(serialize = "MNV")]
    Mnv,
    #[strum(serialize = "BND")]
    Breakend,
    #[strum(serialize = "INV")]
    Inversion,
    #[strum(serialize = "DUP")]
    Duplication,
    #[strum(serialize = "REP")]
    Replacement,
    #[strum(serialize = "REF")]
    None, // site with no suggested alternative allele
}

impl From<&str> for VariantType {
    fn from(string: &str) -> VariantType {
        match string {
            "INS" => VariantType::Insertion(None),
            "DEL" => VariantType::Deletion(None),
            "SNV" => VariantType::Snv,
            "REF" => VariantType::None,
            "INV" => VariantType::Inversion,
            "DUP" => VariantType::Duplication,
            "REP" => VariantType::Replacement,
            "BND" => VariantType::Breakend,
            _ => panic!("bug: given string does not describe a valid variant type"),
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) enum Variant {
    Deletion(u64),
    Insertion(Vec<u8>),
    Snv(u8),
    Mnv(Vec<u8>),
    Breakend {
        ref_allele: Vec<u8>,
        spec: Vec<u8>,
        event: Vec<u8>,
    },
    Inversion(u64),
    Duplication(u64),
    Replacement {
        ref_allele: Vec<u8>,
        alt_allele: Vec<u8>,
    },
    None,
}

impl Variant {
    pub(crate) fn is_breakend(&self) -> bool {
        matches!(self, Variant::Breakend { .. })
    }

    pub(crate) fn is_type(&self, vartype: &VariantType) -> bool {
        match (self, vartype) {
            (&Variant::Deletion(l), &VariantType::Deletion(Some(ref range))) => {
                l >= range.start && l < range.end
            }
            (&Variant::Insertion(_), &VariantType::Insertion(Some(ref range))) => {
                self.len() >= range.start && self.len() < range.end
            }
            (&Variant::Deletion(_), &VariantType::Deletion(None)) => true,
            (&Variant::Insertion(_), &VariantType::Insertion(None)) => true,
            (&Variant::Snv(_), &VariantType::Snv) => true,
            (&Variant::Mnv(_), &VariantType::Mnv) => true,
            (&Variant::None, &VariantType::None) => true,
            (&Variant::Breakend { .. }, &VariantType::Breakend) => true,
            (&Variant::Inversion { .. }, &VariantType::Inversion) => true,
            (&Variant::Duplication { .. }, &VariantType::Duplication) => true,
            (&Variant::Replacement { .. }, &VariantType::Replacement) => true,
            _ => false,
        }
    }

    pub(crate) fn to_type(&self) -> VariantType {
        match self {
            Variant::Deletion(_) => VariantType::Deletion(None),
            Variant::Insertion(_) => VariantType::Insertion(None),
            Variant::Snv(_) => VariantType::Snv,
            Variant::Mnv(_) => VariantType::Mnv,
            Variant::Breakend { .. } => VariantType::Breakend,
            Variant::Inversion(_) => VariantType::Inversion,
            Variant::Duplication(_) => VariantType::Duplication,
            Variant::Replacement { .. } => VariantType::Replacement,
            Variant::None => VariantType::None,
        }
    }

    pub(crate) fn len(&self) -> u64 {
        match *self {
            Variant::Deletion(l) => l,
            Variant::Insertion(ref s) => s.len() as u64,
            Variant::Snv(_) => 1,
            Variant::Mnv(ref alt) => alt.len() as u64,
            Variant::Breakend { .. } => 1,
            Variant::Inversion(l) => l,
            Variant::Duplication(l) => l,
            Variant::Replacement { ref alt_allele, .. } => alt_allele.len() as u64,
            Variant::None => 1,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::variants::evidence::observation::{
        Observation, ObservationBuilder, ReadPosition, Strand,
    };
    use bio_types::sequence::SequenceReadPairOrientation;

    use bio::stats::LogProb;

    pub(crate) fn observation(
        prob_mapping: LogProb,
        prob_alt: LogProb,
        prob_ref: LogProb,
    ) -> Observation<ReadPosition> {
        ObservationBuilder::default()
            .prob_mapping_mismapping(prob_mapping)
            .prob_alt(prob_alt)
            .prob_ref(prob_ref)
            .prob_missed_allele(prob_ref.ln_add_exp(prob_alt) - LogProb(2.0_f64.ln()))
            .prob_sample_alt(LogProb::ln_one())
            .prob_overlap(LogProb::ln_one())
            .read_orientation(SequenceReadPairOrientation::None)
            .read_position(ReadPosition::Some)
            .strand(Strand::Both)
            .softclipped(false)
            .has_alt_indel_operations(false)
            .paired(true)
            .prob_hit_base(LogProb::from(0.01f64.ln()))
            .build()
            .unwrap()
    }
}
