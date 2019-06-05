use std::cell::RefCell;
use std::collections::HashMap;
use std::convert::TryFrom;

use serde_yaml;

pub mod formula;
pub mod vaftree;

use crate::errors;
pub use crate::grammar::formula::{Formula, VAFRange, VAFSpectrum, VAFUniverse};
pub use crate::grammar::vaftree::VAFTree;
use crate::model;

#[derive(Deserialize, Getters)]
#[get = "pub"]
pub struct Scenario {
    // map of events
    events: HashMap<String, Formula>,
    // map of samples
    samples: HashMap<String, Sample>,
    #[serde(skip)]
    sample_idx: RefCell<Option<HashMap<String, usize>>>,
}

impl Scenario {
    pub fn idx(&self, sample: &str) -> Option<usize> {
        if self.sample_idx.borrow().is_none() {
            self.sample_idx.borrow_mut().get_or_insert(
                self.samples()
                    .keys()
                    .enumerate()
                    .map(|(i, s)| (s.to_owned(), i))
                    .collect(),
            );
        }
        self.sample_idx
            .borrow()
            .as_ref()
            .unwrap()
            .get(sample)
            .map(|idx| *idx)
    }

    pub fn sort_samples_by_idx(&self, samples: &mut Vec<model::sample::Sample>) {
        samples.sort_by_key(|sample| self.idx(sample.name()));
    }

    pub fn vaftrees(&self) -> Result<HashMap<String, VAFTree>, errors::FormulaError> {
        self.events()
            .iter()
            .map(|(name, formula)| {
                let normalized = formula.normalize(self)?;
                let vaftree = VAFTree::new(&normalized, self)?;
                Ok((name.to_owned(), vaftree))
            })
            .collect()
    }
}

impl TryFrom<&'static str> for Scenario {
    type Error = serde_yaml::Error;

    fn try_from(yaml: &str) -> Result<Self, Self::Error> {
        serde_yaml::from_str(yaml)
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
    /// possible VAFs of given sample
    universe: VAFUniverse,
}

#[derive(Deserialize, Getters)]
#[get = "pub"]
pub struct Contamination {
    /// name of contaminating sample
    by: String,
    /// fraction of contamination
    fraction: f64,
}
