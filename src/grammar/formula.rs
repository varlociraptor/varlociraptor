use std::ops;
use std::collections::{BTreeSet, HashMap, VecDeque};
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

impl<T: Hash + Eq> Formula<T> {

    /// Negate formula.
    pub fn negate(self, vaf_universe: &HashMap<T, BTreeSet<VAFSpectrum>>) -> Formula<T> {
        match self {
            Formula::Conjunction { operands } => Formula::Disjunction { operands },
            Formula::Disjunction { operands } => Formula::Conjunction { operands },
            Formula::Negation { operand } => *operand,
            Formula::Atom { sample, vafs } => {
                let universe = vaf_universe.get(&sample).expect(format!("No VAF-Universe given for sample {}", sample));
                let mut disjunction = Vec::new();
                match vafs {
                    spectrum @ VAFSpectrum::Set(vafs) => {
                        let uvaf_stack: VecDeque<_> = universe.iter().cloned().collect();
                        loop {
                            let uvafs = universe.pop();
                            match uvafs {
                                VAFSpectrum::Set(uvafs) => {
                                    let mut difference = uvafs.clone();
                                    difference.remove(vafs);
                                    if !difference.is_empty() {
                                        disjunction.push(VAFSpectrum::Set(difference));
                                    }
                                }
                                uspectrum @ VAFSpectrum::Range(urange) => {
                                    for vaf in vafs {
                                        if urange.contains(vaf) {
                                            let (left_urange, right_urange) = urange.split_at(vaf);
                                            uvaf_stack.push(VAFSpectrum::Range(right_urange));
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
                        for uvafs in universe {
                            match uvafs {
                                uspectrum @ VAFSpectrum::Set(uvafs) => {
                                    for uvaf in uvafs {
                                        if !range.contains(uvaf) {
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
            }
        }
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
            Formula::Atom { sample, vafs: VAFSpectrum::singleton(AlleleFreq(vaf)) }
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
