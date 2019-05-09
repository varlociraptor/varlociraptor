use std::collections::HashMap;

pub mod formula;

use crate::grammar::formula::Formula;

#[derive(Deserialize, Getters)]
#[get = "pub"]
pub struct Scenario {
    // map of events
    events: HashMap<String, Formula>,
    // map of samples
    samples: HashMap<String, Sample>,
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
