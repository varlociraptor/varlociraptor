// Copyright 2017 Nicholas Ingolia
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Trait shared across sequence locations -- spliced, contiguous, or
//! single-position.

use std::ops::Neg;

use crate::annot::contig::Contig;
use crate::annot::pos::Pos;

use crate::strand::*;

/// A trait for a sequence location -- a defined region on a named
/// chromosome (or other reference sequence), which may also have
/// defined strand information. The trait is generic over the type of
/// identifier for the reference sequence (allowing owned strings,
/// sequence IDs, and other options) and the strand information
/// (allowing type-level distinction between stranded and unstranded
/// annotations).
pub trait Loc {
    type RefID;
    type Strand;

    /// Name of the reference sequence (chromosome name, etc.)
    fn refid(&self) -> &Self::RefID;
    /// Starting (lowest, left-most, 5'-most) position on the
    /// reference sequence (0-based).
    fn start(&self) -> isize;
    /// Length of the region
    fn length(&self) -> usize;
    /// `Strand` of the position
    fn strand(&self) -> Self::Strand
    where
        Self::Strand: Copy;

    /// Map a sequence position on a reference sequence _into_ a
    /// relative position within an annotated location on the
    /// reference sequence.
    ///
    /// The first position of the annotated location is mapped to a
    /// position at 0, the next position of the annotated location is
    /// mapped to a position at 1, and so forth. The annotated
    /// location must have a known strandedness, which is taken into
    /// account. For reverse-strand annotations, the 3'-most position
    /// on the reference sequence is mapped to 0, and the strandedness
    /// of the position is reversed. When the sequence position lies
    /// on a different named reference sequence than the annotated
    /// location, or doesn't fall within the annotated location, then
    /// `None` is returned.
    ///
    /// This function serves as an inverse of @pos_outof.
    fn pos_into<T>(&self, pos: &Pos<Self::RefID, T>) -> Option<Pos<(), T>>
    where
        Self::RefID: Eq,
        Self::Strand: Into<ReqStrand> + Copy,
        T: Neg<Output = T> + Copy;

    /// Map a relative position within an annotated location _out of_
    /// that location onto the enclosing reference sequence.
    ///
    /// Position 0 within the annotated location is mapped to the
    /// first position of the annotated location, position 1 is mapped
    /// to the subsequent position, and so forth. The annotated
    /// location must have a known strandedness, which is taken into
    /// account. For reverse-strand annotations, position 0 is mapped
    /// to the 3'-most position of the reference sequence. When the
    /// sequence position is either negative, or greater than the
    /// length of the annotated location, then `None` is returned. The
    /// reference name for the sequence position is discarded; the
    /// mapped position receives a clone of the annotation's reference
    /// sequence name.
    ///
    /// This function serves as an inverse of @pos_into.
    fn pos_outof<Q, T>(&self, pos: &Pos<Q, T>) -> Option<Pos<Self::RefID, T>>
    where
        Self::RefID: Clone,
        Self::Strand: Into<ReqStrand> + Copy,
        T: Neg<Output = T> + Copy;

    fn contig_intersection<T>(&self, other: &Contig<Self::RefID, T>) -> Option<Self>
    where
        Self: ::std::marker::Sized,
        Self::RefID: PartialEq + Clone,
        Self::Strand: Copy;

    /// Contiguous sequence location that fully covers the location.
    fn contig(&self) -> Contig<Self::RefID, Self::Strand>
    where
        Self::RefID: Clone,
        Self::Strand: Copy,
    {
        Contig::new(
            self.refid().clone(),
            self.start(),
            self.length(),
            self.strand(),
        )
    }

    /// The first `Pos` in a location, on the annotated strand.
    ///
    /// The first position in a zero-length annotation will be the
    /// starting postion. This is the same as the first position in a
    /// length-1 annotation, on either strand.
    fn first_pos(&self) -> Pos<Self::RefID, Self::Strand>
    where
        Self::RefID: Clone,
        Self::Strand: Into<ReqStrand> + Copy,
    {
        match self.strand().into() {
            ReqStrand::Forward => Pos::new(self.refid().clone(), self.start(), self.strand()),
            ReqStrand::Reverse => {
                if self.length() == 0 {
                    Pos::new(self.refid().clone(), self.start(), self.strand())
                } else {
                    Pos::new(
                        self.refid().clone(),
                        self.start() + (self.length() as isize) - 1,
                        self.strand(),
                    )
                }
            }
        }
    }

    /// The last `Pos` in a location, on the annotated strand.
    ///
    /// The last position in a zero-length annotation will be the
    /// starting postion. This is the same as the last position in a
    /// length-1 annotation, on either strand.
    fn last_pos(&self) -> Pos<Self::RefID, Self::Strand>
    where
        Self::RefID: Clone,
        Self::Strand: Into<ReqStrand> + Copy,
    {
        match self.strand().into() {
            ReqStrand::Forward => {
                if self.length() == 0 {
                    Pos::new(self.refid().clone(), self.start(), self.strand())
                } else {
                    Pos::new(
                        self.refid().clone(),
                        self.start() + (self.length() as isize) - 1,
                        self.strand(),
                    )
                }
            }
            ReqStrand::Reverse => Pos::new(self.refid().clone(), self.start(), self.strand()),
        }
    }
}
