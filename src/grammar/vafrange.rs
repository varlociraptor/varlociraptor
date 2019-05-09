use std::fmt;
use std::ops;

use pest::Parser;
use serde::de;
use serde::Deserialize;

use crate::model::ContinuousAlleleFreqs;

#[derive(Parser)]
#[grammar = "grammar/vafrange.pest"]
pub struct VAFRangeParser;

#[derive(Clone, Debug)]
pub enum VAFRange {
    Singleton(f64),
    Exclusive(ops::Range<f64>),
    LeftExclusive(ops::Range<f64>),
    RightExclusive(ops::Range<f64>),
    Inclusive(ops::Range<f64>),
}

impl Into<ContinuousAlleleFreqs> for VAFRange {
    fn into(self) -> ContinuousAlleleFreqs {
        match self {
            VAFRange::Singleton(vaf) => ContinuousAlleleFreqs::singleton(vaf),
            VAFRange::Exclusive(vafs) => ContinuousAlleleFreqs::exclusive(vafs),
            VAFRange::Inclusive(vafs) => ContinuousAlleleFreqs::inclusive(vafs),
            VAFRange::LeftExclusive(vafs) => ContinuousAlleleFreqs::left_exclusive(vafs),
            VAFRange::RightExclusive(vafs) => ContinuousAlleleFreqs::right_exclusive(vafs),
        }
    }
}

impl<'de> Deserialize<'de> for VAFRange {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de>,
    {
        deserializer.deserialize_string(VAFRangeVisitor)
    }
}

struct VAFRangeVisitor;

impl<'de> de::Visitor<'de> for VAFRangeVisitor {
    type Value = VAFRange;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str(
            "a string of the form 'VAF', '[VAF, VAF]', ']VAF, VAF[', \
             '[VAF, VAF[', ']VAF, VAF]', with VAF between 0.0 and 1.0",
        )
    }

    fn visit_str<E>(self, v: &str) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        let res = VAFRangeParser::parse(Rule::vafdef, v);
        if let Ok(mut pairs) = res {
            let pair = pairs.next().unwrap();
            Ok(match pair.as_rule() {
                Rule::vaf => VAFRange::Singleton(pair.as_str().parse().unwrap()),
                Rule::vafrange => {
                    let mut inner = pair.into_inner();
                    let left = inner.next().unwrap().as_str();
                    let lower = inner.next().unwrap().as_str().parse().unwrap();
                    let upper = inner.next().unwrap().as_str().parse().unwrap();
                    let right = inner.next().unwrap().as_str();

                    match (left, right) {
                        ("[", "]") => VAFRange::Inclusive(lower..upper),
                        ("]", "]") => VAFRange::LeftExclusive(lower..upper),
                        ("[", "[") => VAFRange::RightExclusive(lower..upper),
                        ("]", "[") => VAFRange::Exclusive(lower..upper),
                        _ => unreachable!(),
                    }
                }
                Rule::bound => {
                    return Err(de::Error::invalid_value(
                        serde::de::Unexpected::Other("invalid VAF or VAF range"),
                        &self,
                    ))
                }
                Rule::vafdef => unreachable!("vafdef is not recursive"),
            })
        } else {
            Err(de::Error::invalid_value(
                serde::de::Unexpected::Other("invalid VAF or VAF range"),
                &self,
            ))
        }
    }
}
