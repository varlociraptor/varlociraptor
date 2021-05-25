use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::convert::TryFrom;
use std::fs::File;
use std::io::Read;
use std::ops::{Deref, DerefMut};
use std::path::Path;
use std::string::ToString;
use std::sync::Mutex;

use anyhow::{Context, Result};
use serde_yaml;
use vec_map::VecMap;

pub(crate) mod formula;
pub(crate) mod vaftree;

use crate::errors;
pub(crate) use crate::grammar::formula::{Formula, VAFRange, VAFSpectrum, VAFUniverse};
pub(crate) use crate::grammar::vaftree::VAFTree;
use crate::variants::model::{AlleleFreq, VariantType};

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
            .ok_or_else(|| errors::Error::EmptyObservations.into())
    }

    pub(crate) fn first_not_none_mut(&mut self) -> Result<&mut T> {
        self.iter_not_none_mut()
            .next()
            .ok_or_else(|| errors::Error::EmptyObservations.into())
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

impl<T> From<Vec<T>> for SampleInfo<T> {
    fn from(vec: Vec<T>) -> Self {
        SampleInfo { inner: vec }
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

#[derive(Derefable, PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Hash, Deserialize)]
pub(crate) struct ExpressionIdentifier(#[deref] String);

impl ToString for ExpressionIdentifier {
    fn to_string(&self) -> String {
        self.0.clone()
    }
}

#[derive(Deserialize, Getters)]
#[get = "pub(crate)"]
#[serde(deny_unknown_fields)]
pub(crate) struct Scenario {
    // map of reusable expressions
    #[serde(default)]
    expressions: HashMap<ExpressionIdentifier, Formula>,
    // map of events
    events: BTreeMap<String, Formula>,
    // map of samples
    samples: BTreeMap<String, Sample>,
    #[serde(skip)]
    sample_idx: Mutex<Option<HashMap<String, usize>>>,
    #[serde(default)]
    species: Option<Species>,
}

impl Scenario {
    pub(crate) fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut scenario_content = String::new();
        File::open(path)?.read_to_string(&mut scenario_content)?;

        let mut scenario: Self = serde_yaml::from_str(&scenario_content)?;

        let mut event_expressions = HashMap::new();

        // register all events as expressions
        for (name, formula) in scenario.events() {
            let identifier = ExpressionIdentifier(name.clone());
            if !scenario.expressions().contains_key(&identifier) {
                event_expressions.insert(identifier, formula.clone());
            }
        }
        let absent_identifier = ExpressionIdentifier("absent".to_owned());
        if !scenario.expressions.contains_key(&absent_identifier) {
            event_expressions.insert(absent_identifier, Formula::absent(&scenario));
        }
        scenario.expressions.extend(event_expressions);

        Ok(scenario)
    }

    pub(crate) fn variant_type_fractions(&self) -> VariantTypeFraction {
        self.species()
            .as_ref()
            .map_or(VariantTypeFraction::default(), |species| {
                species.variant_type_fractions().clone()
            })
    }

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
        info!("Preprocessing events for contig {}", contig);
        self.events()
            .iter()
            .map(|(name, formula)| {
                let normalized = formula
                    .normalize(self, contig)
                    .with_context(|| format!("invalid event definition for {}", name))?;
                info!("    {}: {}", name, normalized);
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

#[derive(Deserialize)]
#[serde(untagged)]
pub(crate) enum PloidyDefinition {
    Simple(u32),
    Map(HashMap<String, u32>),
}

impl PloidyDefinition {
    pub(crate) fn contig_ploidy(&self, contig: &str) -> Result<u32> {
        Ok(match self {
            PloidyDefinition::Simple(ploidy) => *ploidy,
            PloidyDefinition::Map(map) => match map.get(contig) {
                Some(ploidy) => *ploidy,
                None => map.get("all").map(|ploidy| *ploidy).ok_or_else(|| {
                    errors::Error::PloidyContigNotFound {
                        contig: contig.to_owned(),
                    }
                })?,
            },
        })
    }
}

#[derive(Deserialize)]
#[serde(untagged)]
pub(crate) enum SexPloidyDefinition {
    Generic(PloidyDefinition),
    Specific(HashMap<Sex, PloidyDefinition>),
}

impl SexPloidyDefinition {
    pub(crate) fn contig_ploidy(&self, sex: Option<Sex>, contig: &str) -> Result<u32> {
        match (self, sex) {
            (SexPloidyDefinition::Generic(p), _) => p.contig_ploidy(contig),
            (SexPloidyDefinition::Specific(p), Some(s)) => p.get(&s).map_or_else(
                || {
                    Err(errors::Error::InvalidPriorConfiguration {
                        msg: format!("ploidy definition for {} not found", s),
                    }
                    .into())
                },
                |p| p.contig_ploidy(contig),
            ),
            (SexPloidyDefinition::Specific(_), None) => {
                Err(errors::Error::InvalidPriorConfiguration {
                    msg: "sex specific ploidy definition found but no sex specified in sample"
                        .to_owned(),
                }
                .into())
            }
        }
    }
}

#[derive(Deserialize, Getters)]
#[get = "pub(crate)"]
#[serde(deny_unknown_fields)]
pub(crate) struct Species {
    #[serde(default)]
    heterozygosity: Option<f64>,
    #[serde(default, rename = "germline-mutation-rate")]
    germline_mutation_rate: Option<f64>,
    #[serde(default, rename = "somatic-effective-mutation-rate")]
    somatic_effective_mutation_rate: Option<f64>,
    #[serde(default, rename = "variant-fractions")]
    variant_type_fractions: VariantTypeFraction,
    #[serde(default)]
    ploidy: Option<SexPloidyDefinition>,
    #[serde(default)]
    #[serde(rename = "genome-size")]
    genome_size: Option<f64>,
}

impl Species {
    pub(crate) fn contig_ploidy(&self, contig: &str, sex: Option<Sex>) -> Result<Option<u32>> {
        if let Some(ploidy) = &self.ploidy {
            Ok(Some(ploidy.contig_ploidy(sex, contig)?))
        } else {
            Ok(None)
        }
    }
}

fn default_indel_fraction() -> f64 {
    0.0125
}

fn default_mnv_fraction() -> f64 {
    0.001
}

fn default_sv_fraction() -> f64 {
    0.01
}

/// Mutation rate reduciton factors (relative to SNVs) for variant types.
/// Defaults:
/// * mnvs: 0.001 (see https://www.nature.com/articles/s41467-019-12438-5)
/// * indels: 0.0125 (see https://gatk.broadinstitute.org/hc/en-us/articles/360036826431-HaplotypeCaller, reduction in heterozygosity)
/// * svs: 0.001 (predicted several hundred times less frequent that SNVs: https://doi.org/10.1038/s41588-018-0107-y)
#[derive(Deserialize, Getters, Clone, Debug)]
#[get = "pub(crate)"]
#[serde(deny_unknown_fields)]
pub(crate) struct VariantTypeFraction {
    #[serde(default = "indel")]
    indel: f64,
    #[serde(default = "mnv")]
    mnv: f64,
    #[serde(default = "sv")]
    sv: f64,
}

impl VariantTypeFraction {
    pub(crate) fn get(&self, variant_type: &VariantType) -> f64 {
        match variant_type {
            VariantType::Insertion(_) | VariantType::Deletion(_) | VariantType::Replacement => {
                self.indel
            }
            VariantType::MNV => self.mnv,
            VariantType::Inversion | VariantType::Breakend | VariantType::Duplication => self.sv,
            _ => 1.0,
        }
    }
}

impl Default for VariantTypeFraction {
    fn default() -> Self {
        VariantTypeFraction {
            indel: default_indel_fraction(),
            mnv: default_mnv_fraction(),
            sv: default_sv_fraction(),
        }
    }
}

fn default_resolution() -> usize {
    100
}

#[derive(Deserialize, Getters)]
#[serde(deny_unknown_fields)]
pub(crate) struct Sample {
    /// optional contamination
    #[get = "pub(crate)"]
    contamination: Option<Contamination>,
    /// grid point resolution for integration over continuous allele frequency ranges
    #[serde(default = "default_resolution")]
    #[get = "pub(crate)"]
    resolution: usize,
    /// possible VAFs of given sample
    #[serde(default)]
    #[get = "pub(crate)"]
    universe: Option<UniverseDefinition>,
    #[serde(default, rename = "somatic-effective-mutation-rate")]
    somatic_effective_mutation_rate: Option<f64>,
    #[serde(default, rename = "germline-mutation-rate")]
    germline_mutation_rate: Option<f64>,
    #[serde(default)]
    #[get = "pub(crate)"]
    ploidy: Option<PloidyDefinition>,
    #[get = "pub(crate)"]
    inheritance: Option<Inheritance>,
    #[serde(default)]
    sex: Option<Sex>,
}

impl Sample {
    pub(crate) fn has_uniform_prior(&self) -> bool {
        self.universe.is_some()
    }

    pub(crate) fn contig_universe(
        &self,
        contig: &str,
        species: &Option<Species>,
    ) -> Result<VAFUniverse> {
        if let Some(universe) = &self.universe {
            Ok(match universe {
                UniverseDefinition::Simple(ref universe) => universe.clone(),
                UniverseDefinition::Map(ref map) => match map.get(contig) {
                    Some(universe) => universe.clone(),
                    None => map
                        .get("all")
                        .ok_or_else(|| errors::Error::UniverseContigNotFound {
                            contig: contig.to_owned(),
                        })?
                        .clone(),
                },
            })
        } else {
            let ploidy_derived_spectrum = |ploidy| -> BTreeSet<AlleleFreq> {
                (0..=ploidy)
                    .map(|n_alt| {
                        if ploidy > 0 {
                            AlleleFreq(n_alt as f64 / ploidy as f64)
                        } else {
                            assert_eq!(n_alt, 0);
                            AlleleFreq(0.0)
                        }
                    })
                    .collect()
            };
            Ok(
                match (
                    self.contig_ploidy(contig, species)?,
                    self.somatic_effective_mutation_rate.is_some(),
                ) {
                    (Some(ploidy), false) => {
                        let mut universe = VAFUniverse::default();
                        universe.insert(VAFSpectrum::Set(ploidy_derived_spectrum(ploidy)));
                        universe
                    }
                    (Some(ploidy), true) => {
                        let ploidy_spectrum = ploidy_derived_spectrum(ploidy);

                        let mut universe = VAFUniverse::default();

                        let mut last = ploidy_spectrum.iter().next().unwrap();
                        for vaf in ploidy_spectrum.iter().skip(1) {
                            universe.insert(VAFSpectrum::Range(VAFRange::builder()
                                .inner(*last..*vaf)
                                .left_exclusive(true)
                                .right_exclusive(true)
                                .build()
                            ));
                            last = vaf;
                        }
                        universe.insert(VAFSpectrum::Set(ploidy_spectrum));
                        universe
                    }
                    (None, true) => {
                        let mut universe = VAFUniverse::default();
                        universe.insert(VAFSpectrum::Range(VAFRange::builder()
                            .inner(AlleleFreq(0.0)..AlleleFreq(1.0))
                            .left_exclusive(false)
                            .right_exclusive(false)
                            .build()
                        ));
                        universe
                    }
                    (None, false) => return Err(errors::Error::InvalidPriorConfiguration{
                        msg: "sample needs to define either universe, ploidy or somatic_mutation_rate".to_owned(),
                    }
                    .into()),
                },
            )
        }
    }

    pub(crate) fn contig_ploidy(
        &self,
        contig: &str,
        species: &Option<Species>,
    ) -> Result<Option<u32>> {
        if let Some(ploidy) = &self.ploidy {
            Ok(Some(ploidy.contig_ploidy(contig)?))
        } else {
            species
                .as_ref()
                .map_or(Ok(None), |species| species.contig_ploidy(contig, self.sex))
        }
    }

    pub(crate) fn germline_mutation_rate(&self, species: &Option<Species>) -> Option<f64> {
        if let Some(rate) = self.germline_mutation_rate {
            Some(rate)
        } else {
            species
                .as_ref()
                .map_or(None, |species| species.germline_mutation_rate)
        }
    }

    pub(crate) fn somatic_effective_mutation_rate(&self, species: &Option<Species>) -> Option<f64> {
        if let Some(rate) = self.somatic_effective_mutation_rate {
            Some(rate)
        } else {
            species
                .as_ref()
                .map_or(None, |species| species.somatic_effective_mutation_rate)
        }
    }
}

#[derive(Deserialize, Getters)]
#[get = "pub(crate)"]
#[serde(deny_unknown_fields)]
pub(crate) struct Contamination {
    /// name of contaminating sample
    by: String,
    /// fraction of contamination
    fraction: f64,
}

#[derive(
    Display,
    Debug,
    Clone,
    Serialize,
    Deserialize,
    IntoStaticStr,
    EnumVariantNames,
    PartialEq,
    Eq,
    Hash,
)]
#[strum(serialize_all = "kebab_case")]
pub(crate) enum Inheritance {
    #[serde(rename = "mendelian")]
    Mendelian { from: (String, String) },
    #[serde(rename = "clonal")]
    Clonal { from: String, somatic: bool },
    #[serde(rename = "subclonal")]
    Subclonal {
        from: String,
        origin: SubcloneOrigin,
    },
}

#[derive(
    Display,
    Debug,
    Clone,
    Copy,
    Serialize,
    Deserialize,
    IntoStaticStr,
    EnumVariantNames,
    PartialEq,
    Eq,
    Hash,
)]
#[strum(serialize_all = "kebab_case")]
pub(crate) enum SubcloneOrigin {
    #[serde(rename = "single-cell")]
    SingleCell,
    #[serde(rename = "multi-cell")]
    MultiCell,
}

#[derive(
    Display,
    Debug,
    Clone,
    Copy,
    Serialize,
    Deserialize,
    IntoStaticStr,
    EnumVariantNames,
    PartialEq,
    Eq,
    Hash,
)]
pub(crate) enum Sex {
    #[serde(rename = "male")]
    Male,
    #[serde(rename = "female")]
    Female,
}

#[derive(Deserialize, Debug)]
#[serde(untagged)]
pub(crate) enum UniverseDefinition {
    Map(BTreeMap<String, VAFUniverse>),
    Simple(VAFUniverse),
}
