// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Data types for strand information on annotations.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::fmt::{self, Display, Formatter};
use std::ops::Neg;
use std::str::FromStr;
use thiserror::Error;

/// Strand information.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

impl Strand {
    /// Returns a `Strand` enum representing the given char.
    ///
    /// The mapping is as follows:
    ///     * '+', 'f', or 'F' becomes `Strand::Forward`
    ///     * '-', 'r', or 'R' becomes `Strand::Reverse`
    ///     * '.', '?' becomes `Strand::Unknown`
    ///     * Any other inputs will return an `Err(StrandError::InvalidChar)`
    pub fn from_char(strand_char: &char) -> Result<Strand, StrandError> {
        match *strand_char {
            '+' | 'f' | 'F' => Ok(Strand::Forward),
            '-' | 'r' | 'R' => Ok(Strand::Reverse),
            '.' | '?' => Ok(Strand::Unknown),
            invalid => Err(StrandError::InvalidChar(invalid)),
        }
    }

    /// Symbol denoting the strand. By convention, in BED and GFF
    /// files, the forward strand is `+`, the reverse strand is `-`,
    /// and unknown or unspecified strands are `.`.
    pub fn strand_symbol(&self) -> &str {
        match *self {
            Strand::Forward => "+",
            Strand::Reverse => "-",
            Strand::Unknown => ".",
        }
    }

    pub fn is_unknown(&self) -> bool {
        matches!(*self, Strand::Unknown)
    }
}

#[allow(clippy::match_like_matches_macro)]
impl PartialEq for Strand {
    /// Returns true if both are `Forward` or both are `Reverse`, otherwise returns false.
    fn eq(&self, other: &Strand) -> bool {
        match (self, other) {
            (&Strand::Forward, &Strand::Forward) => true,
            (&Strand::Reverse, &Strand::Reverse) => true,
            _ => false,
        }
    }
}

impl Neg for Strand {
    type Output = Strand;
    fn neg(self) -> Strand {
        match self {
            Strand::Forward => Strand::Reverse,
            Strand::Reverse => Strand::Forward,
            Strand::Unknown => Strand::Unknown,
        }
    }
}

#[allow(clippy::match_like_matches_macro)]
impl Same for Strand {
    fn same(&self, s1: &Self) -> bool {
        match (*self, *s1) {
            (Strand::Forward, Strand::Forward) => true,
            (Strand::Reverse, Strand::Reverse) => true,
            (Strand::Unknown, Strand::Unknown) => true,
            _ => false,
        }
    }
}

impl Display for Strand {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        f.write_str(self.strand_symbol())
    }
}

impl FromStr for Strand {
    type Err = StrandError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" | "(+)" => Ok(Strand::Forward),
            "-" | "(-)" => Ok(Strand::Reverse),
            "." | "" => Ok(Strand::Unknown),
            _ => Err(StrandError::ParseError),
        }
    }
}

impl From<ReqStrand> for Strand {
    fn from(rstr: ReqStrand) -> Self {
        match rstr {
            ReqStrand::Forward => Strand::Forward,
            ReqStrand::Reverse => Strand::Reverse,
        }
    }
}

impl From<Option<ReqStrand>> for Strand {
    fn from(orstr: Option<ReqStrand>) -> Self {
        match orstr {
            Some(ReqStrand::Forward) => Strand::Forward,
            Some(ReqStrand::Reverse) => Strand::Reverse,
            None => Strand::Unknown,
        }
    }
}

impl From<NoStrand> for Strand {
    fn from(_: NoStrand) -> Self {
        Strand::Unknown
    }
}

/// Strand information for annotations that require a strand.
#[derive(Debug, Clone, Hash, PartialEq, Eq, Ord, PartialOrd, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ReqStrand {
    Forward,
    Reverse,
}

impl ReqStrand {
    /// Returns a `ReqStrand` enum representing the given char.
    ///
    /// The mapping is as follows:
    ///     * '+', 'f', or 'F' becomes `Strand::Forward`
    ///     * '-', 'r', or 'R' becomes `Strand::Reverse`
    ///     * Any other inputs will return an `Err(StrandError::InvalidChar)`
    pub fn from_char(strand_char: &char) -> Result<ReqStrand, StrandError> {
        match *strand_char {
            '+' | 'f' | 'F' => Ok(ReqStrand::Forward),
            '-' | 'r' | 'R' => Ok(ReqStrand::Reverse),
            invalid => Err(StrandError::InvalidChar(invalid)),
        }
    }

    /// Symbol denoting the strand. By convention, in BED and GFF
    /// files, the forward strand is `+` and the reverse strand is `-`.
    pub fn strand_symbol(&self) -> &str {
        match *self {
            ReqStrand::Forward => "+",
            ReqStrand::Reverse => "-",
        }
    }

    /// Convert the (optional) strand of some other annotation
    /// according to this strand. That is, reverse the strand of the
    /// other annotation for `ReqStrand::Reverse` and leave it
    /// unchanged for `ReqStrand::Forward`.
    ///
    /// # Arguments
    ///
    /// * `x` is the strand information from some other annotation.
    ///
    /// ```
    /// use bio_types::strand::{ReqStrand,Strand};
    /// assert_eq!(ReqStrand::Forward.on_strand(Strand::Reverse),
    ///            ReqStrand::Reverse.on_strand(Strand::Forward));
    /// ```
    pub fn on_strand<T>(&self, x: T) -> T
    where
        T: Neg<Output = T>,
    {
        match self {
            ReqStrand::Forward => x,
            ReqStrand::Reverse => -x,
        }
    }
}

impl Same for ReqStrand {
    fn same(&self, s1: &Self) -> bool {
        self == s1
    }
}

impl Display for ReqStrand {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        f.write_str(self.strand_symbol())
    }
}

impl FromStr for ReqStrand {
    type Err = StrandError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" | "(+)" => Ok(ReqStrand::Forward),
            "-" | "(-)" => Ok(ReqStrand::Reverse),
            _ => Err(StrandError::ParseError),
        }
    }
}

impl From<Strand> for Option<ReqStrand> {
    fn from(strand: Strand) -> Option<ReqStrand> {
        match strand {
            Strand::Forward => Some(ReqStrand::Forward),
            Strand::Reverse => Some(ReqStrand::Reverse),
            Strand::Unknown => None,
        }
    }
}

impl From<NoStrand> for Option<ReqStrand> {
    fn from(_: NoStrand) -> Option<ReqStrand> {
        None
    }
}

impl Neg for ReqStrand {
    type Output = ReqStrand;
    fn neg(self) -> ReqStrand {
        match self {
            ReqStrand::Forward => ReqStrand::Reverse,
            ReqStrand::Reverse => ReqStrand::Forward,
        }
    }
}

/// Strand information for annotations that definitively have no
/// strand information.
#[derive(Debug, Clone, Hash, PartialEq, Eq, Ord, PartialOrd, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum NoStrand {
    Unknown,
}

impl Neg for NoStrand {
    type Output = NoStrand;
    fn neg(self) -> NoStrand {
        match self {
            NoStrand::Unknown => NoStrand::Unknown,
        }
    }
}

impl Same for NoStrand {
    fn same(&self, _s1: &Self) -> bool {
        true
    }
}

impl FromStr for NoStrand {
    type Err = StrandError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Ok(NoStrand::Unknown),
            _ => Err(StrandError::ParseError),
        }
    }
}

impl Display for NoStrand {
    fn fmt(&self, _f: &mut Formatter) -> fmt::Result {
        Ok(())
    }
}

/// Equality-like operator for comparing strand information. Unknown
/// strands are not equal, but they are the "same" as other unknown
/// strands.
pub trait Same {
    /// Indicate when two strands are the "same" -- two
    /// unknown/unspecified strands are the "same" but are not equal.
    fn same(&self, other: &Self) -> bool;
}

impl<T> Same for Option<T>
where
    T: Same,
{
    fn same(&self, s1: &Self) -> bool {
        match (self, s1) {
            (&Option::None, &Option::None) => true,
            (&Option::Some(ref x), &Option::Some(ref x1)) => x.same(x1),
            (_, _) => false,
        }
    }
}

#[derive(Error, Debug)]
pub enum StrandError {
    #[error("invalid character for strand conversion: {0:?}: can not be converted to a Strand")]
    InvalidChar(char),
    #[error("error parsing strand")]
    ParseError,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strand() {
        assert_eq!(Strand::from_char(&'+').unwrap(), Strand::Forward);
        assert_eq!(Strand::from_char(&'-').unwrap(), Strand::Reverse);
        assert!(Strand::from_char(&'.').unwrap().is_unknown());
        assert!(Strand::from_char(&'o').is_err());
        assert_eq!(Strand::Forward.strand_symbol(), "+");
        assert_eq!(Strand::Reverse.strand_symbol(), "-");
        assert_eq!(Strand::Unknown.strand_symbol(), ".");
    }

    #[test]
    fn test_req_strand() {
        assert_eq!(ReqStrand::from_char(&'+').unwrap(), ReqStrand::Forward);
        assert_eq!(ReqStrand::from_char(&'-').unwrap(), ReqStrand::Reverse);
        assert!(ReqStrand::from_char(&'o').is_err());
        assert_eq!(ReqStrand::Forward.strand_symbol(), "+");
        assert_eq!(ReqStrand::Reverse.strand_symbol(), "-");
    }
}
