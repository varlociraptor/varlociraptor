use std::cmp;
use std::ops::Range;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

pub type Position = u64;
pub type Length = u64;

pub trait AbstractInterval {
    /// Identifier for a genomic contig, e.g., a chromosome
    fn contig(&self) -> &str;
    /// Interval on the contig
    fn range(&self) -> Range<Position>;
    /// Return true if interval contains given locus.
    fn contains<L>(&self, locus: L) -> bool
    where
        L: AbstractLocus,
    {
        self.contig() == locus.contig()
            && locus.pos() >= self.range().start
            && locus.pos() < self.range().end
    }
}
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(new, Debug, PartialEq, Eq, Clone, Hash)]
pub struct Interval {
    contig: String,
    range: Range<Position>,
}

impl PartialOrd for Interval {
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        Some(self.contig.cmp(&other.contig).then_with(|| {
            self.range
                .start
                .cmp(&other.range.start)
                .then_with(|| self.range.end.cmp(&other.range.end))
        }))
    }
}

impl Ord for Interval {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl Interval {
    /// Mutable reference to interval on the contig
    pub fn range_mut(&mut self) -> &mut Range<Position> {
        &mut self.range
    }
}

impl AbstractInterval for Interval {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn range(&self) -> Range<Position> {
        self.range.clone()
    }
}

pub trait AbstractLocus {
    /// Identifier for a genomic contig, e.g., a chromosome
    fn contig(&self) -> &str;
    /// Position on the contig
    fn pos(&self) -> Position;
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(new, Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Hash)]
pub struct Locus {
    contig: String,
    pos: Position,
}

impl Locus {
    /// Mutable reference to position.
    pub fn pos_mut(&mut self) -> &mut Position {
        &mut self.pos
    }
}

impl AbstractLocus for Locus {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn pos(&self) -> Position {
        self.pos
    }
}
