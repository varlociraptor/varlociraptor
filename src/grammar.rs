use std::ops;

pub struct Grammar {
    events: Vec<Event>,
    samples: Vec<Sample>,
}

pub struct Event {
    name: String,
    vafs: Vec<Range>,
}

pub enum Range {
    Singleton(f64),
    Exclusive(ops::Range<f64>),
    LeftExclusive(ops::Range<f64>),
    RightExclusive(ops::Range<f64>),
    Inclusive(ops::Range<f64>),
}

impl Into<ContinuousAlleleFreqs> for Range {
    fn into(&self) -> ContinuousAlleleFreqs {
        match self {
            Range::Singleton(vaf) => ContinuousAlleleFreqs::singleton(vaf),
            Range::Exclusive(vafs) => ContinuousAlleleFreqs::exclusive(vafs),
            Range::Inclusive(vafs) => ContinuousAlleleFreqs::inclusive(vafs),
            Range::LeftExclusive(vafs) => ContinuousAlleleFreqs::left_exclusive(vafs),
            Range::RightExclusive(vafs) => ContinuousAlleleFreqs::right_exclusive(vafs),
        }
    }
}

impl Into<Vec<ContinuousAlleleFreqs>> for Vec<Range> {
    fn into(&self) -> Vec<ContinuousAlleleFreqs> {
        self.iter().map(|range| range.into()).collect()
    }
}

pub struct Sample {
    name: String,
    contamination: Option<Contamination>,
}

pub struct Contamination {
    by: String,
    fraction: f64,
}
