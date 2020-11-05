use std::collections::{BTreeMap, HashMap};
use std::convert::TryFrom;
use std::ops::{Deref, DerefMut};
use std::sync::Mutex;

use anyhow::Result;
use vec_map::VecMap;

pub(crate) mod formula;
pub(crate) mod vaftree;

use crate::errors;
pub(crate) use crate::grammar::formula::{Formula, VAFSpectrum, VAFUniverse};
pub(crate) use crate::grammar::vaftree::VAFTree;

/// Container for arbitrary sample information.
/// Use `varlociraptor::grammar::Scenario::sample_info()` to create it.
#[derive(Clone, Debug)]
pub(crate) struct SampleInfo<T> {
    inner: Vec<T>,
}

impl<T> SampleInfo<T> {
    /// Map to other value type.
    pub(crate) fn map<U, F: Fn(&T) -> U>(&self, f: F) -> SampleInfo<U> {
        SampleInfo {
            inner: self.inner.iter().map(f).collect(),
        }
    }

    pub(crate) fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.inner.iter_mut()
    }
}

impl<T> SampleInfo<Option<T>> {
    pub(crate) fn first_not_none(&self) -> Result<&T> {
        self.iter_not_none()
            .next()
            .ok_or(errors::Error::EmptyObservations.into())
    }

    pub(crate) fn first_not_none_mut(&mut self) -> Result<&mut T> {
        self.iter_not_none_mut()
            .next()
            .ok_or(errors::Error::EmptyObservations.into())
    }

    pub(crate) fn iter_not_none(&self) -> impl Iterator<Item = &T> {
        self.inner.iter().filter_map(|item| item.as_ref())
    }

    pub(crate) fn iter_not_none_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.inner.iter_mut().filter_map(|item| item.as_mut())
    }
}

impl<T> Default for SampleInfo<T> {
    fn default() -> Self {
        SampleInfo {
            inner: Vec::default(),
        }
    }
}

impl<T> Deref for SampleInfo<T> {
    type Target = Vec<T>;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<T> DerefMut for SampleInfo<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.inner
    }
}

/// Builder for `SampleInfo`.
#[derive(new, Debug)]
pub(crate) struct SampleInfoBuilder<T> {
    #[new(default)]
    inner: VecMap<T>,
    sample_idx: HashMap<String, usize>,
}

impl<T> SampleInfoBuilder<T> {
    pub(crate) fn push(mut self, sample_name: &str, value: T) -> Self {
        let idx = *self
            .sample_idx
            .get(sample_name)
            .expect("unknown sample name, it does not occur in the scenario");
        self.inner.insert(idx, value);

        self
    }

    pub(crate) fn build(self) -> SampleInfo<T> {
        SampleInfo {
            inner: self.inner.into_iter().map(|(_, v)| v).collect(),
        }
    }
}

#[derive(Deserialize, Getters)]
#[get = "pub"]
pub(crate) struct Scenario {
    // map of events
    events: BTreeMap<String, Formula>,
    // map of samples
    samples: BTreeMap<String, Sample>,
    #[serde(skip)]
    sample_idx: Mutex<Option<HashMap<String, usize>>>,
}

impl Scenario {
    pub(crate) fn sample_info<T>(&self) -> SampleInfoBuilder<T> {
        let mut sample_idx = self.sample_idx.lock().unwrap();
        if sample_idx.is_none() {
            sample_idx.get_or_insert(
                self.samples()
                    .keys()
                    .enumerate()
                    .map(|(i, s)| (s.to_owned(), i))
                    .collect(),
            );
        }
        SampleInfoBuilder::new(sample_idx.as_ref().unwrap().clone())
    }

    pub(crate) fn idx(&self, sample: &str) -> Option<usize> {
        let mut sample_idx = self.sample_idx.lock().unwrap();
        if sample_idx.is_none() {
            sample_idx.get_or_insert(
                self.samples()
                    .keys()
                    .enumerate()
                    .map(|(i, s)| (s.to_owned(), i))
                    .collect(),
            );
        }
        sample_idx.as_ref().unwrap().get(sample).copied()
    }

    pub(crate) fn vaftrees(&self, contig: &str) -> Result<HashMap<String, VAFTree>> {
        self.events()
            .iter()
            .map(|(name, formula)| {
                let normalized = formula.normalize(self, contig)?;
                let vaftree = VAFTree::new(&normalized, self, contig)?;
                Ok((name.to_owned(), vaftree))
            })
            .collect()
    }
}

impl<'a> TryFrom<&'a str> for Scenario {
    type Error = serde_yaml::Error;

    fn try_from(yaml: &str) -> Result<Self, Self::Error> {
        serde_yaml::from_str(yaml)
    }
}

#[derive(Deserialize, Getters)]
#[get = "pub"]
pub(crate) struct Sample {
    /// optional contamination
    contamination: Option<Contamination>,
    /// grid point resolution for integration over continuous allele frequency ranges
    resolution: usize,
    /// possible VAFs of given sample
    universe: UniverseDefinition,
}

impl Sample {
    pub(crate) fn contig_universe(&self, contig: &str) -> Result<&VAFUniverse> {
        Ok(match self.universe {
            UniverseDefinition::Simple(ref universe) => universe,
            UniverseDefinition::Map(ref map) => match map.get(contig) {
                Some(universe) => universe,
                None => map
                    .get("all")
                    .ok_or_else(|| errors::Error::UniverseContigNotFound {
                        contig: contig.to_owned(),
                    })?,
            },
        })
    }
}

#[derive(Deserialize, Getters)]
#[get = "pub"]
pub(crate) struct Contamination {
    /// name of contaminating sample
    by: String,
    /// fraction of contamination
    fraction: f64,
}

#[derive(Deserialize)]
#[serde(untagged)]
pub(crate) enum UniverseDefinition {
    Map(BTreeMap<String, VAFUniverse>),
    Simple(VAFUniverse),
}
