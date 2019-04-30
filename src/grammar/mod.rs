use std::collections::HashMap;

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

#[derive(Deserialize, Getters)]
#[get = "pub"]
pub struct Event {
    /// Map from sample/group name to allele freq range
    vafs: HashMap<String, VAFRange>,
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
