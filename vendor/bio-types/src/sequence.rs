// Copyright 2021 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use strum_macros::{AsRefStr, Display};
use SequenceReadPairOrientation::None;

/// A DNA base
pub type Base = u8;
/// An amino acid
pub type AminoAcid = u8;
/// A biological sequence
pub type Sequence = Vec<u8>;

pub trait SequenceRead {
    /// Read name.
    fn name(&self) -> &[u8];
    /// Base at position `i` in the read.
    fn base(&self, i: usize) -> u8;
    /// Base quality at position `i` in the read.
    fn base_qual(&self, i: usize) -> u8;
    /// Read length.
    fn len(&self) -> usize;
    /// Return `true` if read is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Representation of sequence read pair orientation
/// (e.g. F1R2 means that the forward read comes first on the reference contig,
/// followed by the reverse read, on the same contig).
///
/// This enum can be pretty-printed into a readable string repesentation:
///
/// ```rust
/// use bio_types::sequence::SequenceReadPairOrientation;
///
/// // format into string
/// println!("{}", SequenceReadPairOrientation::F1R2);
/// // obtain string via `AsRef<&'static str>`
/// assert_eq!(SequenceReadPairOrientation::R1F2.as_ref(), "R1F2");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, AsRefStr, Display)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum SequenceReadPairOrientation {
    F1R2,
    F2R1,
    R1F2,
    R2F1,
    F1F2,
    R1R2,
    F2F1,
    R2R1,
    None,
}

impl Default for SequenceReadPairOrientation {
    fn default() -> Self {
        None
    }
}
