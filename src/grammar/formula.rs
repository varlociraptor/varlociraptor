use std::ops;
use std::collections::{BTreeSet, HashMap};
use std::fmt;
use std::cmp::{Ordering, Ord, PartialOrd};

use itertools::Itertools;
use pest::iterators::Pair;
use pest::Parser;
use serde::de;
use serde::Deserialize;

use crate::model::AlleleFreq;

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

impl Formula<String> {
    pub fn to_sample_idx(&self, sample_idx: &HashMap<String, usize>) -> Formula<usize> {
        let recurse_to_operands = |operands: &[Box<Formula<String>>]| {
            operands.iter().map(|f| Box::new(f.to_sample_idx(sample_idx))).collect_vec()
        };
        match self {
            &Formula::Conjunction { ref operands } => Formula::Conjunction { operands: recurse_to_operands(operands) },
            &Formula::Disjunction { ref operands } => Formula::Disjunction { operands: recurse_to_operands(operands) },
            &Formula::Negation { ref operand } => Formula::Negation { operand: Box::new(operand.to_sample_idx(sample_idx)) },
            &Formula::Atom { ref sample, ref vafs } => Formula::Atom { sample: *sample_idx.get(sample).unwrap(), vafs: vafs.clone() },
        }
    }
}

impl Formula<usize> {
    pub fn absent(n_samples: usize) -> Formula<usize> {
        Formula::Conjunction { operands: (0..n_samples).into_iter().map(|sample| Box::new(Formula::Atom { sample, vafs: VAFSpectrum::Singleton(AlleleFreq(0.0)) })).collect_vec() }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord) ]
pub enum VAFSpectrum {
    Singleton(AlleleFreq),
    Set(BTreeSet<AlleleFreq>),
    Range(VAFRange)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum VAFRange {
    Exclusive(ops::Range<AlleleFreq>),
    LeftExclusive(ops::Range<AlleleFreq>),
    RightExclusive(ops::Range<AlleleFreq>),
    Inclusive(ops::Range<AlleleFreq>),
}

impl ops::Deref for VAFRange {
    type Target = ops::Range<AlleleFreq>;

    fn deref(&self) -> &ops::Range<AlleleFreq> {
        match self {
            VAFRange::Exclusive(ref range) => range,
            VAFRange::LeftExclusive(ref range) => range,
            VAFRange::RightExclusive(ref range) => range,
            VAFRange::Inclusive(ref range) => range,
        }
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

fn parse_formula<E>(pair: Pair<Rule>) -> Result<Formula<String>, E>
where
    E: de::Error
{
    Ok(match pair.as_rule() {
        Rule::vaf => {
            let inner = pair.into_inner();
            let sample = inner.next().unwrap().as_str().to_owned();
            let vaf = inner.next().unwrap().as_str().parse().unwrap();
            Formula::Atom { sample, vafs: VAFSpectrum::Singleton(AlleleFreq(vaf)) }
        }
        Rule::vafrange => {
            let inner = pair.into_inner();
            let sample = inner.next().unwrap().as_str().to_owned();
            let left = inner.next().unwrap().as_str();
            let lower = inner.next().unwrap().as_str().parse().unwrap();
            let upper = inner.next().unwrap().as_str().parse().unwrap();
            let right = inner.next().unwrap().as_str();

            let range = match (left, right) {
                ("[", "]") => VAFRange::Inclusive(lower..upper),
                ("]", "]") => VAFRange::LeftExclusive(lower..upper),
                ("[", "[") => VAFRange::RightExclusive(lower..upper),
                ("]", "[") => VAFRange::Exclusive(lower..upper),
                _ => unreachable!(),
            };
            Formula::Atom { sample, vafs: VAFSpectrum::Range(range) }
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
