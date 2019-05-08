use std::collections::HashMap;
use std::ops::Deref;

pub mod vafrange;

use crate::grammar::vafrange::VAFRange;

#[derive(Deserialize, Getters)]
#[get = "pub"]
pub struct Scenario {
    // map of events
    events: HashMap<String, Event>,
    // map of samples
    samples: HashMap<String, Sample>,
}

#[derive(Deserialize)]
pub struct Event (
    /// Map from sample/group name to allele freq range
    HashMap<String, VAFRange>,
);

impl Deref for Event {
    type Target = HashMap<String, VAFRange>;

    fn deref(&self) -> &HashMap<String, VAFRange> {
        &self.0
    }
}

#[derive(Deserialize, Getters)]
#[get = "pub"]
pub struct Sample {
    /// optional group name
    group: Option<String>,
    /// optional contamination
    contamination: Option<Contamination>,
    /// grid point resolution for integration over continuous allele frequency ranges
    resolution: usize,
}

#[derive(Deserialize, Getters)]
#[get = "pub"]
pub struct Contamination {
    /// name of contaminating sample
    by: String,
    /// fraction of contamination
    fraction: f64,
}
