//! Data types for positions and regions on named sequences
//! (e.g. chromosomes), useful for annotating features in a genome.
//! For example, these data types let you represent that _TMA22_ is on
//! chromosome X, positions 461,829-462,426, on the forward strand. They
//! also allow coordinate math on these annotations, e.g., that
//! position chrX:461,839 is +10 within _TMA22_ and vice versa.
//!
//! This module provides three concrete data types to represent a
//! single position ([`Pos`](pos/Pos.t.html)), a contiguous region
//! ([`Contig`](contig/Contig.t.html)), or a "spliced" region
//! ([`Spliced`](spliced/Spliced.t.html)) consisting of one or more
//! exons separated by introns.  All three data types implement a
//! location trait [`Loc`](loc/Loc.t.html).
//!
//! These data types are generic over the data type used to "name" the
//! annotated reference sequence (e.g., the chromosome name). It's
//! possible to use an owned `String`, an interned `Rc<String>`, or an
//! integer sequence identifier like the "target id" field in a BAM
//! file.
//!
//! These data types are also generic over the kind of strand
//! information in the annotation. This allows annotations with
//! _required_ strand annotation
//! ([`ReqStrand`](../strand/enum.ReqStrand.html)), _optional_ strand
//! annotation ([`Strand`](../strand/enum.Strand.html)), or _no_
//! strand annotation ([`NoStrand`](../strand/enum.NoStrand.html)).
//!
//! The example below shows how to create the _TMA22_ annotation and
//! find where chrX:461,839 falls within this gene.
//! ```
//! # use std::str::FromStr;
//! # use bio_types::annot::ParseAnnotError;
//! # fn try_main() -> Result<(), ParseAnnotError> {
//! use bio_types::annot::contig::Contig;
//! use bio_types::annot::loc::Loc;
//! use bio_types::annot::pos::Pos;
//! use bio_types::strand::{ReqStrand,NoStrand};
//! let tma22: Contig<String, ReqStrand> = Contig::from_str("chrX:461829-462426(+)")?;
//! let p0: Pos<String, NoStrand> = Pos::from_str("chrX:461839")?;
//! let p0_into = tma22.pos_into(&p0).unwrap_or_else(|| panic!("p0 not within TMA22"));
//! assert!(p0_into.pos() == 10);
//! # Ok(())
//! # }
//! # fn main() { try_main().unwrap(); }
//! ```

use crate::strand;
use thiserror::Error;

pub mod contig;
pub mod loc;
pub mod pos;
pub mod refids;
pub mod spliced;

// Errors that arise in parsing annotations.
#[derive(Error, Debug)]
pub enum ParseAnnotError {
    #[error("Annotation string does not match regex")]
    BadAnnot,
    #[error("Integer parsing error")]
    ParseInt(#[from] ::std::num::ParseIntError),
    #[error("Strand parsing error")]
    ParseStrand(#[from] strand::StrandError),
    #[error("Bad splicing structure")]
    Splicing(#[from] spliced::SplicingError),
    #[error("Ending position < starting position")]
    EndBeforeStart,
}

// Errors that arise in maniuplating annotations
#[derive(Error, Debug)]
pub enum AnnotError {
    #[error("No strand information")]
    NoStrand,
    #[error("Invalid splicing structure")]
    BadSplicing,
}
