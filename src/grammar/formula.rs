use std::ops;
use std::collections::{BTreeSet, HashMap, VecDeque};
use std::fmt;
use std::cmp::{Ordering, Ord, PartialOrd};
use std::hash::Hash;
use std::fmt::Display;

use itertools::Itertools;
use pest::iterators::{Pair,Pairs};
use pest::Parser;
use serde::de;
use serde::Deserialize;

use crate::model::AlleleFreq;
use crate::errors;


#[derive(Parser)]
#[grammar = "grammar/formula.pest"]
pub struct FormulaParser;

#[derive(PartialEq, PartialOrd, Ord, Eq, Clone, Debug)]
pub enum Formula<T> {
    Conjunction {
        operands: Vec<Box<Formula<T>>>
    },
    Disjunction {
        operands: Vec<Box<Formula<T>>>
    },
    Negation {
        operand: Box<Formula<T>>
    },
    Atom {
        sample: T,
        vafs: VAFSpectrum
    }
}

impl<T: Hash + Eq + Display + Clone> Formula<T> {

    /// Negate formula.
    pub fn negate(self, vaf_universe: &HashMap<T, VAFUniverse>) -> Formula<T> {
        match self {
            Formula::Conjunction { operands } => Formula::Disjunction { operands: operands.iter().map(|o| Box::new(o.negate(vaf_universe))).collect() },
            Formula::Disjunction { operands } => Formula::Conjunction { operands: operands.iter().map(|o| Box::new(o.negate(vaf_universe))).collect() },
            Formula::Negation { operand } => *operand,
            Formula::Atom { sample, vafs } => {
                let universe = vaf_universe.get(&sample).expect(&format!("bug: no VAF-Universe given for sample {}", sample));
                let mut disjunction = Vec::new();
                match vafs {
                    spectrum @ VAFSpectrum::Set(vafs) => {
                        let uvaf_stack: VecDeque<_> = universe.iter().cloned().collect();
                        while let Some(uvafs) = uvaf_stack.pop_front() {
                            match uvafs {
                                VAFSpectrum::Set(uvafs) => {
                                    let mut difference: BTreeSet<_> = uvafs.difference(&vafs).cloned().collect();
                                    if !difference.is_empty() {
                                        disjunction.push(VAFSpectrum::Set(difference));
                                    }
                                }
                                uspectrum @ VAFSpectrum::Range(urange) => {
                                    for vaf in vafs {
                                        if urange.contains(vaf) {
                                            let (left_urange, right_urange) = urange.split_at(vaf);
                                            uvaf_stack.push_back(VAFSpectrum::Range(right_urange));
                                            disjunction.push(VAFSpectrum::Range(left_urange));
                                        } else {
                                            disjunction.push(uspectrum.clone());
                                        }
                                    }
                                }
                            }
                        }
                    },
                    VAFSpectrum::Range(range) => {
                        for uvafs in universe.iter() {
                            match uvafs {
                                uspectrum @ VAFSpectrum::Set(uvafs) => {
                                    for uvaf in uvafs {
                                        if !range.contains(*uvaf) {
                                            disjunction.push(uspectrum.clone());
                                        }
                                    }
                                }
                                uspectrum @ VAFSpectrum::Range(urange) => {
                                    match range.overlap(urange) {
                                        VAFRangeOverlap::Contained => {
                                            let left = urange.split_at(range.start).0;
                                            let right = urange.split_at(range.end).1;
                                            disjunction.push(VAFSpectrum::Range(left));
                                            disjunction.push(VAFSpectrum::Range(right));
                                        }
                                        VAFRangeOverlap::End => {
                                            disjunction.push(VAFSpectrum::Range(urange.split_at(range.end).1));
                                        }
                                        VAFRangeOverlap::Start => {
                                            disjunction.push(VAFSpectrum::Range(urange.split_at(range.start).0));
                                        }
                                        VAFRangeOverlap::None => {
                                            disjunction.push(uspectrum.clone())
                                        }
                                        VAFRangeOverlap::Contains => {
                                            ()
                                        }
                                    }
                                }
                            }
                        }

                    }
                }
                Formula::Disjunction { operands: disjunction.into_iter().map(|vafs| Box::new(Formula::Atom { sample: sample.clone(), vafs })).collect() }
            }
        }
    }

    pub fn atomize_negations(self, vaf_universe: &HashMap<T, VAFUniverse>) -> Formula<T> {
        match self {
            Formula::Negation { operand } => operand.negate(vaf_universe).atomize_negations(vaf_universe),
            Formula::Atom { .. } => self,
            Formula::Conjunction { operands } => Formula::Conjunction { operands: operands.into_iter().map(|o| Box::new(o.atomize_negations(vaf_universe))).collect() },
            Formula::Disjunction { operands } => Formula::Disjunction { operands: operands.into_iter().map(|o| Box::new(o.atomize_negations(vaf_universe))).collect() },
        }
    }
}

impl Formula<String> {
    pub fn to_sample_idx(&self, sample_idx: &HashMap<String, usize>) -> Result<Formula<usize>, errors::FormulaError> {
        let recurse_to_operands = |operands: &[Box<Formula<String>>]| {
            let operands: Result<Vec<_>, _> = operands.iter().map(|f| Ok(Box::new(f.to_sample_idx(sample_idx)?))).collect();
            operands
        };
        Ok(match self {
            &Formula::Conjunction { ref operands } => Formula::Conjunction { operands: recurse_to_operands(operands)? },
            &Formula::Disjunction { ref operands } => Formula::Disjunction { operands: recurse_to_operands(operands)? },
            &Formula::Negation { ref operand } => Formula::Negation { operand: Box::new(operand.to_sample_idx(sample_idx)?) },
            &Formula::Atom { ref sample, ref vafs } => {
                if let Some(idx) = sample_idx.get(sample) {
                    Formula::Atom { sample: *idx, vafs: vafs.clone() }
                } else {
                    return Err(errors::FormulaError::InvalidSampleName { name: sample.to_owned() });
                }
            },
        })
    }
}

impl Formula<usize> {
    pub fn absent(n_samples: usize) -> Formula<usize> {
        Formula::Conjunction { operands: (0..n_samples).into_iter().map(|sample| Box::new(Formula::Atom { sample, vafs: VAFSpectrum::singleton(AlleleFreq(0.0)) })).collect_vec() }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord) ]
pub enum VAFSpectrum {
    Set(BTreeSet<AlleleFreq>),
    Range(VAFRange)
}

impl VAFSpectrum {
    pub fn singleton(vaf: AlleleFreq) -> Self {
        let mut set = BTreeSet::new();
        set.insert(vaf);
        VAFSpectrum::Set(set)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VAFRange {
    inner: ops::Range<AlleleFreq>,
    left_exclusive: bool,
    right_exclusive: bool,
}

pub enum VAFRangeOverlap {
    Contained,
    Contains,
    End,
    Start,
    None
}

impl VAFRange {
    pub fn contains(&self, vaf: AlleleFreq) -> bool {
        match (self.left_exclusive, self.right_exclusive) {
            (true, true) => self.start < vaf && self.end > vaf,
            (true, false) => self.start < vaf && self.end >= vaf,
            (false, true) => self.start <= vaf && self.end > vaf,
            (true, true) => self.start <= vaf && self.end >= vaf,
        }
    }

    pub fn split_at(&self, vaf: AlleleFreq) -> (VAFRange, VAFRange) {
        assert!(self.contains(vaf), "bug: split_at is only defined if given VAF is contained in range");
        let left = VAFRange { inner: self.start..vaf, left_exclusive: self.left_exclusive, right_exclusive: true };
        let right = VAFRange { inner: vaf..self.end, left_exclusive: true, right_exclusive: self.right_exclusive };
        (left, right)
    }

    pub fn overlap(&self, vafs: &VAFRange) -> VAFRangeOverlap {
        let range = *self;
        let other_range = **vafs;
        let start_is_right_of_start = match (self.left_exclusive, self.right_exclusive) {
            (true, true)  => range.start >= other_range.start,
            (true, false) => range.start >= other_range.start,
            (false, true) => range.start > other_range.start,
            (false, false) => range.start >= other_range.start,
        };
        let end_is_left_of_end = match (self.left_exclusive, self.right_exclusive) {
            (true, true)  => range.end <= other_range.end,
            (true, false) => range.end <= other_range.end,
            (false, true) => range.end < other_range.end,
            (false, false) => range.end <= other_range.end,
        };
        if range.end < other_range.start || range.start >= other_range.end {
            VAFRangeOverlap::None
        } else {
            match (start_is_right_of_start, end_is_left_of_end) {
                (true, true) => VAFRangeOverlap::Contained,
                (true, false) => VAFRangeOverlap::Start,
                (false, true) => VAFRangeOverlap::End,
                (false, false) => VAFRangeOverlap::Contains,
            }
        }
    }

    pub fn observable_min(&self, n_obs: usize) -> AlleleFreq {
        if n_obs < 10 {
            self.start
        } else {
            let obs_count = Self::expected_observation_count(self.start, n_obs);
            let adjust_allelefreq = |obs_count: f64| AlleleFreq(obs_count.ceil() / n_obs as f64);

            if self.left_exclusive && obs_count % 1.0 == 0.0 {
                // We are left exclusive and need to find a supremum from the right.

                let adjusted_end = self.observable_max(n_obs);

                for offset in &[1.0, 0.0] {
                    let adjusted_obs_count = obs_count + offset;
                    let adjusted_start = adjust_allelefreq(adjusted_obs_count);
                    if *adjusted_start <= 1.0 && adjusted_start <= adjusted_end {
                        return adjusted_start;
                    }
                }
            }

            adjust_allelefreq(obs_count)
        }
    }

    pub fn observable_max(&self, n_obs: usize) -> AlleleFreq {
        assert!(
            *self.end != 0.0,
            "bug: observable_max may not be called if end=0.0."
        );
        if n_obs < 10 {
            self.end
        } else {
            let mut obs_count = Self::expected_observation_count(self.end, n_obs);
            if self.right_exclusive && obs_count % 1.0 == 0.0 {
                obs_count -= 1.0;
            }
            AlleleFreq(obs_count.floor() / n_obs as f64)
        }
    }

    fn expected_observation_count(freq: AlleleFreq, n_obs: usize) -> f64 {
        n_obs as f64 * *freq
    }
}

impl ops::Deref for VAFRange {
    type Target = ops::Range<AlleleFreq>;

    fn deref(&self) -> &ops::Range<AlleleFreq> {
        &self.inner
    }
}


impl Ord for VAFRange {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.start.cmp(&other.start) {
            Ordering::Equal => self.end.cmp(&other.end),
            ord @ _ => ord,
        }
    }
}

impl PartialOrd for VAFRange {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone)]
pub struct VAFUniverse(BTreeSet<VAFSpectrum>);

impl ops::Deref for VAFUniverse {
    type Target = BTreeSet<VAFSpectrum>;

    fn deref(&self) -> &BTreeSet<VAFSpectrum> {
        &self.0
    }
}

impl<'de> Deserialize<'de> for VAFUniverse {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de>,
    {
        deserializer.deserialize_string(VAFUniverseVisitor)
    }
}


struct VAFUniverseVisitor;

impl<'de> de::Visitor<'de> for VAFUniverseVisitor {
    type Value = VAFUniverse;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str(
            "a valid disjunction of possible VAFs (see https://varlociraptor.github.io/docs/calling)",
        )
    }

    fn visit_str<E>(self, v: &str) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        let res = FormulaParser::parse(Rule::universe, v);
        if let Ok(mut pairs) = res {
            let pair = pairs.next().unwrap();
            match pair.as_rule() {
                Rule::universe => {
                    let inner = pair.into_inner();
                    let mut operands = BTreeSet::new();
                    loop {
                        let pair = inner.next().unwrap();
                        match pair.as_rule() {
                            Rule::vaf => {
                                let inner = pair.into_inner();
                                operands.insert(parse_vaf(inner));
                            }
                            Rule::vafrange => {
                                let inner = pair.into_inner();
                                operands.insert(parse_vafrange(inner));
                            }
                            _ => unreachable!()
                        }
                        if inner.next().is_none() {
                            break;
                        }
                    }
                    Ok(VAFUniverse(operands))
                }
            }

        } else {
            Err(de::Error::invalid_value(
                serde::de::Unexpected::Other("invalid VAF formula"),
                &self,
            ))
        }
    }
}


impl<'de> Deserialize<'de> for Formula<String> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de>,
    {
        deserializer.deserialize_string(FormulaVisitor)
    }
}

struct FormulaVisitor;

impl<'de> de::Visitor<'de> for FormulaVisitor {
    type Value = Formula<String>;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str(
            "a valid VAF formula (see https://varlociraptor.github.io/docs/calling)",
        )
    }

    fn visit_str<E>(self, v: &str) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        let res = FormulaParser::parse(Rule::formula, v);
        if let Ok(mut pairs) = res {
            let pair = pairs.next().unwrap();
            parse_formula(pair)
        } else {
            Err(de::Error::invalid_value(
                serde::de::Unexpected::Other("invalid VAF formula"),
                &self,
            ))
        }
    }
}

fn parse_vaf(inner: Pairs<Rule>) -> VAFSpectrum {
    let vaf = inner.next().unwrap().as_str().parse().unwrap();
    VAFSpectrum::singleton(AlleleFreq(vaf))
}

fn parse_vafrange(inner: Pairs<Rule>) -> VAFSpectrum {
    let left = inner.next().unwrap().as_str();
    let lower = inner.next().unwrap().as_str().parse().unwrap();
    let upper = inner.next().unwrap().as_str().parse().unwrap();
    let right = inner.next().unwrap().as_str();

    let range = VAFRange {
        inner: lower..upper,
        left_exclusive: left == "]",
        right_exclusive: right == "[",
    };

    VAFSpectrum::Range(range)
}

fn parse_formula<E>(pair: Pair<Rule>) -> Result<Formula<String>, E>
where
    E: de::Error
{
    Ok(match pair.as_rule() {
        Rule::vaf => {
            let inner = pair.into_inner();
            let sample = inner.next().unwrap().as_str().to_owned();
            // the colon
            inner.next().unwrap();
            Formula::Atom { sample, vafs: parse_vaf(inner) }
        }
        Rule::vafrange => {
            let inner = pair.into_inner();
            let sample = inner.next().unwrap().as_str().to_owned();
            // the colon
            inner.next().unwrap();
            Formula::Atom { sample, vafs: parse_vafrange(inner) }
        }
        Rule::conjunction => {
            let inner = pair.into_inner();
            let mut operands = Vec::new();
            loop {
                operands.push(Box::new(parse_formula(inner.next().unwrap())?));
                if inner.next().is_none() {
                    break;
                }
            }
            Formula::Conjunction { operands: operands }
        }
        Rule::disjunction => {
            let inner = pair.into_inner();
            let mut operands = Vec::new();
            loop {
                operands.push(Box::new(parse_formula(inner.next().unwrap())?));
                if inner.next().is_none() {
                    break;
                }
            }
            Formula::Disjunction { operands: operands }
        }
        Rule::negation => {
            let inner = pair.into_inner();
            Formula::Negation { operand: Box::new(parse_formula(inner.next().unwrap())?) }
        }
        Rule::formula => unreachable!(),
        Rule::subformula => unreachable!(),
        Rule::vafdef => unreachable!(),
    })
}
