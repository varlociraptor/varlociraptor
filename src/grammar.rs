use std::collections::HashMap;
use std::ops;

use crate::model::ContinuousAlleleFreqs;

pub struct Grammar {
    // list of events
    events: Vec<Event>,
    // list of samples
    samples: Vec<Sample>,
}

pub struct Event {
    /// name of event
    name: String,
    /// Map from sample/group name to allele freq range
    vafs: HashMap<String, Range>,
}

pub enum Range {
    Singleton(f64),
    Exclusive(ops::Range<f64>),
    LeftExclusive(ops::Range<f64>),
    RightExclusive(ops::Range<f64>),
    Inclusive(ops::Range<f64>),
}

impl Into<ContinuousAlleleFreqs> for Range {
    fn into(self) -> ContinuousAlleleFreqs {
        match self {
            Range::Singleton(vaf) => ContinuousAlleleFreqs::singleton(vaf),
            Range::Exclusive(vafs) => ContinuousAlleleFreqs::exclusive(vafs),
            Range::Inclusive(vafs) => ContinuousAlleleFreqs::inclusive(vafs),
            Range::LeftExclusive(vafs) => ContinuousAlleleFreqs::left_exclusive(vafs),
            Range::RightExclusive(vafs) => ContinuousAlleleFreqs::right_exclusive(vafs),
        }
    }
}

pub struct Sample {
    /// sample name
    name: String,
    /// optional group name
    group: Option<String>,
    /// optional contamination
    contamination: Option<Contamination>,
    /// grid point resolution for integration over continuous allele frequency ranges
    resolution: usize,
}

pub struct Contamination {
    /// name of contaminating sample
    by: String,
    /// fraction of contamination
    fraction: f64,
}
