use crate::genome;
use crate::sequence::{Base, Sequence};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A trait for providing variant information. This can e.g. be implemented by file readers.
pub trait AbstractVariant: genome::AbstractLocus {
    fn kind(&self) -> &Kind;
}

/// Possible genomic variants.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub enum Kind {
    SNV(Base),
    MNV(Sequence),
    Insertion(Sequence),
    Deletion(genome::Length),
    Duplication(genome::Length),
    Inversion(genome::Length),
    None,
}

impl Kind {
    /// Return variant length.
    pub fn len(&self) -> genome::Length {
        match *self {
            Kind::SNV(_) => 1,
            Kind::MNV(ref s) => s.len() as u64,
            Kind::Insertion(ref s) => s.len() as u64,
            Kind::Deletion(l) => l,
            Kind::Duplication(l) => l,
            Kind::Inversion(l) => l,
            Kind::None => 1,
        }
    }
    /// Check if length is zero
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
