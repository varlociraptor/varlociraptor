// Copyright 2017 Nicholas Ingolia
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Contiguous region on a named sequence, e.g., chromosome XI
//! 334,915-334,412.

use std::cmp::{max, min};
use std::convert::Into;
use std::fmt::{self, Display, Formatter};
use std::ops::Neg;
use std::str::FromStr;

use regex::Regex;

use crate::annot::loc::Loc;
use crate::annot::pos::Pos;
use crate::annot::*;
use crate::strand::*;

/// Contiguous sequence region on a particular, named sequence (e.g. a
/// chromosome)
///
/// Parameterized over the type of the reference sequence identifier
/// and over the strandedness of the position.
///
/// The display format for a `Contig` is _chr:start-end(+/-/.)_. The
/// boundaries are given as a half-open 0-based interval, like the
/// Rust `Range` and BED format.
///
/// ```
/// # use bio_types::annot::ParseAnnotError;
/// # fn try_main() -> Result<(), Box<ParseAnnotError>> {
/// use bio_types::annot::contig::Contig;
/// use bio_types::strand::ReqStrand;
/// let tma19 = Contig::new("chrXI".to_owned(), 334412, (334916 - 334412), ReqStrand::Reverse);
/// let tma19_str = tma19.to_string();
/// assert_eq!(tma19_str, "chrXI:334412-334916(-)");
/// let tma19_str_loc = tma19_str.parse()?;
/// assert_eq!(tma19, tma19_str_loc);
/// # Ok(())
/// # }
/// # fn main() { try_main().unwrap(); }
/// ```
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct Contig<R, S> {
    refid: R,
    start: isize,
    length: usize,
    strand: S,
}

impl<R, S> Contig<R, S> {
    /// Construct a new sequence contig location
    ///
    /// ```
    /// use std::rc::Rc;
    /// use bio_types::annot::contig::Contig;
    /// use bio_types::strand::ReqStrand;
    /// let chr = Rc::new("chrX".to_owned());
    /// let tma22 = Contig::new(chr, 461829, 462426 - 461829, ReqStrand::Forward);
    /// ```
    pub fn new(refid: R, start: isize, length: usize, strand: S) -> Self {
        Contig {
            refid,
            start,
            length,
            strand,
        }
    }

    /// Construct a new sequence contig location from a starting
    /// position and length.
    ///
    /// In general, the starting position must have a "strandedness",
    /// and reverse-strand starting positions will extend towards
    /// lower coordinates from the starting position.
    ///
    ///
    ///
    /// ```
    /// # use bio_types::annot::AnnotError;
    /// # fn try_main() -> Result<(), Box<AnnotError>> {
    /// use bio_types::annot::contig::Contig;
    /// use bio_types::annot::pos::Pos;
    /// use bio_types::strand::ReqStrand;
    ///
    /// let tma22_first = Pos::new("chrX".to_string(), 461829, ReqStrand::Forward);
    /// let tma22 = Contig::with_first_length(&tma22_first, 462426 - 461829)?;
    /// assert_eq!(tma22.to_string(), "chrX:461829-462426(+)");
    ///
    /// let tma19_first = Pos::new("chrXI".to_string(), 335015, ReqStrand::Reverse);
    /// let tma19 = Contig::with_first_length(&tma19_first, 335016 - 334412)?;
    /// assert_eq!(tma19.to_string(), "chrXI:334412-335016(-)");
    /// # Ok(())
    /// # }
    /// # fn main() { try_main().unwrap(); }
    /// ```
    pub fn with_first_length(pos: &Pos<R, S>, length: usize) -> Result<Self, AnnotError>
    where
        R: Clone,
        S: Into<Option<ReqStrand>> + Copy,
    {
        if length < 2 {
            Ok(Contig {
                refid: pos.refid().clone(),
                start: pos.start(),
                length,
                strand: pos.strand(),
            })
        } else {
            let start = match pos.strand().into() {
                None => Err(AnnotError::NoStrand),
                Some(ReqStrand::Forward) => Ok(pos.start()),
                Some(ReqStrand::Reverse) => Ok(1 + pos.start() - length as isize),
            }?;

            Ok(Contig {
                refid: pos.refid().clone(),
                start,
                length,
                strand: pos.strand(),
            })
        }
    }

    /// Convert into a stranded sequence location on the specified strand
    pub fn into_stranded(self, strand: ReqStrand) -> Contig<R, ReqStrand> {
        Contig {
            refid: self.refid,
            start: self.start,
            length: self.length,
            strand,
        }
    }
}

impl<R> Contig<R, ReqStrand> {
    /// Extend the annotation by `dist` in the upstream direction on the
    /// annotated strand.
    ///
    /// # Arguments
    ///
    /// * `dist` specifies the offset for sliding the position. The
    /// left, 5'-most end of the contig will expand for forward-strand
    /// annotations and the right, 3'-most end will expand for
    /// reverse-strand annotations.
    ///
    /// ```
    /// use bio_types::annot::contig::Contig;
    /// use bio_types::strand::ReqStrand;
    /// let mut tma22 = Contig::new("chrX".to_owned(), 461829, 462426 - 461829, ReqStrand::Forward);
    /// tma22.extend_upstream(100);
    /// assert_eq!(tma22.to_string(), "chrX:461729-462426(+)");
    /// let mut tma19 = Contig::new("chrXI".to_owned(), 334412, 334916 - 334412, ReqStrand::Reverse);
    /// tma19.extend_upstream(100);
    /// assert_eq!(tma19.to_string(), "chrXI:334412-335016(-)");
    /// ```
    pub fn extend_upstream(&mut self, dist: usize) {
        self.length += dist;
        if self.strand == ReqStrand::Forward {
            self.start -= dist as isize;
        }
    }

    /// Extend the annotation by `dist` in the downstream direction on the
    /// annotated strand.
    ///
    /// # Arguments
    ///
    /// * `dist` specifies the offset for sliding the position. The
    /// right, 3'-most end of the contig will expand for
    /// forward-strand annotations and the left, 5'-most end will
    /// expand for reverse-strand annotations.
    ///
    /// ```
    /// use bio_types::annot::contig::Contig;
    /// use bio_types::strand::ReqStrand;
    /// let mut tma22 = Contig::new("chrX".to_owned(), 461829, 462426 - 461829, ReqStrand::Forward);
    /// tma22.extend_downstream(100);
    /// assert_eq!(tma22.to_string(), "chrX:461829-462526(+)");
    /// let mut tma19 = Contig::new("chrXI".to_owned(), 334412, 334916 - 334412, ReqStrand::Reverse);
    /// tma19.extend_downstream(100);
    /// assert_eq!(tma19.to_string(), "chrXI:334312-334916(-)");
    /// ```
    pub fn extend_downstream(&mut self, dist: usize) {
        self.length += dist;
        if self.strand == ReqStrand::Reverse {
            self.start -= dist as isize;
        }
    }
}

impl<R, S> Loc for Contig<R, S> {
    type RefID = R;
    type Strand = S;
    fn refid(&self) -> &R {
        &self.refid
    }
    fn start(&self) -> isize {
        self.start
    }
    fn length(&self) -> usize {
        self.length
    }
    fn strand(&self) -> S
    where
        S: Copy,
    {
        self.strand
    }

    fn pos_into<T>(&self, pos: &Pos<Self::RefID, T>) -> Option<Pos<(), T>>
    where
        Self::RefID: Eq,
        Self::Strand: Into<ReqStrand> + Copy,
        T: Neg<Output = T> + Copy,
    {
        if self.refid != *pos.refid() {
            None
        } else {
            let offset = pos.pos() - self.start;
            if offset < 0 || offset >= self.length as isize {
                None
            } else {
                Some(match self.strand().into() {
                    ReqStrand::Forward => Pos::new((), offset, pos.strand()),
                    ReqStrand::Reverse => {
                        Pos::new((), self.length as isize - (offset + 1), -pos.strand())
                    }
                })
            }
        }
    }

    fn pos_outof<Q, T>(&self, pos: &Pos<Q, T>) -> Option<Pos<Self::RefID, T>>
    where
        Self::RefID: Clone,
        Self::Strand: Into<ReqStrand> + Copy,
        T: Neg<Output = T> + Copy,
    {
        let offset = match self.strand().into() {
            ReqStrand::Forward => pos.pos(),
            ReqStrand::Reverse => self.length as isize - (pos.pos() + 1),
        };

        if offset >= 0 && offset < self.length as isize {
            Some(Pos::new(
                self.refid.clone(),
                self.start + offset,
                self.strand().into().on_strand(pos.strand()),
            ))
        } else {
            None
        }
    }

    fn contig_intersection<T>(&self, contig: &Contig<Self::RefID, T>) -> Option<Self>
    where
        Self::RefID: PartialEq + Clone,
        Self::Strand: Copy,
    {
        if self.refid() != contig.refid() {
            return None;
        }

        let start = max(self.start, contig.start);
        let end = min(
            self.start + self.length as isize,
            contig.start + contig.length as isize,
        );

        if start <= end {
            Some(Self::new(
                self.refid.clone(),
                start,
                (end - start) as usize,
                self.strand,
            ))
        } else {
            None
        }
    }
}

impl<R, S> Display for Contig<R, S>
where
    R: Display,
    S: Display + Clone + Into<Strand>,
{
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(
            f,
            "{}:{}-{}",
            self.refid,
            self.start,
            self.start + self.length as isize
        )?;
        let strand: Strand = self.strand.clone().into();
        if !strand.is_unknown() {
            write!(f, "({})", strand)?;
        }
        Ok(())
    }
}

impl<R, S> FromStr for Contig<R, S>
where
    R: From<String>,
    S: FromStr<Err = StrandError>,
{
    type Err = ParseAnnotError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        lazy_static! {
            static ref CONTIG_RE: Regex = Regex::new(r"^(.*):(\d+)-(\d+)(\([+-]\))?$").unwrap();
        }

        let cap = CONTIG_RE.captures(s).ok_or(ParseAnnotError::BadAnnot)?;

        let start = cap[2].parse::<isize>().map_err(ParseAnnotError::ParseInt)?;
        let end = cap[3].parse::<isize>().map_err(ParseAnnotError::ParseInt)?;
        let strand = cap
            .get(4)
            .map_or("", |m| m.as_str())
            .parse::<S>()
            .map_err(ParseAnnotError::ParseStrand)?;

        if start <= end {
            Ok(Contig::new(
                R::from(cap[1].to_owned()),
                start,
                (end - start) as usize,
                strand,
            ))
        } else {
            Err(ParseAnnotError::EndBeforeStart)
        }
    }
}

impl<R> From<Contig<R, ReqStrand>> for Contig<R, Strand> {
    fn from(x: Contig<R, ReqStrand>) -> Self {
        Contig {
            refid: x.refid,
            start: x.start,
            length: x.length,
            strand: match x.strand {
                ReqStrand::Forward => Strand::Forward,
                ReqStrand::Reverse => Strand::Reverse,
            },
        }
    }
}

impl<R> From<Contig<R, NoStrand>> for Contig<R, Strand> {
    fn from(x: Contig<R, NoStrand>) -> Self {
        Contig {
            refid: x.refid,
            start: x.start,
            length: x.length,
            strand: Strand::Unknown,
        }
    }
}

impl<R> From<Contig<R, Strand>> for Contig<R, NoStrand> {
    fn from(x: Contig<R, Strand>) -> Self {
        Contig {
            refid: x.refid,
            start: x.start,
            length: x.length,
            strand: NoStrand::Unknown,
        }
    }
}

impl<R> From<Contig<R, ReqStrand>> for Contig<R, NoStrand> {
    fn from(x: Contig<R, ReqStrand>) -> Self {
        Contig {
            refid: x.refid,
            start: x.start,
            length: x.length,
            strand: NoStrand::Unknown,
        }
    }
}

/// Default stranded sequence position on a reference sequence named
/// by a `String`.
pub type SeqContigStranded = Contig<String, ReqStrand>;

/// Default unstranded sequence position on a reference sequence named
/// by a `String`
pub type SeqContigUnstranded = Contig<String, NoStrand>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_and_last() {
        let tma22 = "chrX:461829-462426(+)"
            .parse::<SeqContigStranded>()
            .unwrap();
        let first = tma22.first_pos();
        assert_eq!(first.to_string(), "chrX:461829(+)");
        let last = tma22.last_pos();
        assert_eq!(last.to_string(), "chrX:462425(+)");

        let tma19 = "chrXI:334412-334916(-)"
            .parse::<SeqContigStranded>()
            .unwrap();
        let first = tma19.first_pos();
        assert_eq!(first.to_string(), "chrXI:334915(-)");
        let last = tma19.last_pos();
        assert_eq!(last.to_string(), "chrXI:334412(-)");

        let tma22_first = Pos::new("chrX".to_string(), 461829, ReqStrand::Forward);
        let tma22 = Contig::with_first_length(&tma22_first, 462426 - 461829).unwrap();
        assert_eq!(tma22.to_string(), "chrX:461829-462426(+)");

        let tma19_first = Pos::new("chrXI".to_string(), 335015, ReqStrand::Reverse);
        let tma19 = Contig::with_first_length(&tma19_first, 335016 - 334412).unwrap();
        assert_eq!(tma19.to_string(), "chrXI:334412-335016(-)");
    }

    #[test]
    fn into_outof() {
        let tma22 = "chrX:461829-462426(+)"
            .parse::<SeqContigStranded>()
            .unwrap();
        let p0 = "chrX:461829(+)".parse::<Pos<String, ReqStrand>>().unwrap();
        let p0_into = tma22.pos_into(&p0);
        assert!(Some(Pos::new((), 0, ReqStrand::Forward)).same(&p0_into));
        let p0_outof = tma22.pos_outof(&p0_into.unwrap());
        assert!(Some(p0).same(&p0_outof));

        let p0 = "chrX:461839(-)".parse::<Pos<String, ReqStrand>>().unwrap();
        let p0_into = tma22.pos_into(&p0);
        assert!(Some(Pos::new((), 10, ReqStrand::Reverse)).same(&p0_into));
        let p0_outof = tma22.pos_outof(&p0_into.unwrap());
        assert!(Some(p0).same(&p0_outof));

        let p0 = "chrX:462425(+)".parse::<Pos<String, ReqStrand>>().unwrap();
        let p0_into = tma22.pos_into(&p0);
        assert!(Some(Pos::new((), 596, ReqStrand::Forward)).same(&p0_into));
        let p0_outof = tma22.pos_outof(&p0_into.unwrap());
        assert!(Some(p0).same(&p0_outof));

        let p0 = "chrX:461828(+)".parse::<Pos<String, ReqStrand>>().unwrap();
        let p0_into = tma22.pos_into(&p0);
        assert!(None.same(&p0_into));

        let p0 = "chrV:461829(+)".parse::<Pos<String, ReqStrand>>().unwrap();
        let p0_into = tma22.pos_into(&p0);
        assert!(None.same(&p0_into));

        let p0 = "chrV:462426(+)".parse::<Pos<String, ReqStrand>>().unwrap();
        let p0_into = tma22.pos_into(&p0);
        assert!(None.same(&p0_into));
    }

    fn test_contig_ixn(ca_str: &str, cb_str: &str, cab_str: Option<String>) -> () {
        let ca = ca_str.parse::<SeqContigStranded>().unwrap();
        let cb = cb_str.parse::<SeqContigStranded>().unwrap();
        match ca.contig_intersection(&cb) {
            None => assert_eq!(None, cab_str),
            Some(cab) => assert_eq!(Some(cab.to_string()), cab_str),
        };
    }

    #[test]
    fn test_display_fmt() {
        let tma19 = Contig::new(
            "chrXI".to_owned(),
            334412,
            334916 - 334412,
            ReqStrand::Reverse,
        );
        assert_eq!(format!("{}", tma19), "chrXI:334412-334916(-)");
    }

    #[test]
    fn intersection() {
        test_contig_ixn(
            "chrX:461829-462426(+)",
            "chrX:461800-461900(+)",
            Some("chrX:461829-461900(+)".to_owned()),
        );
        test_contig_ixn(
            "chrX:461829-462426(-)",
            "chrX:461800-461900(+)",
            Some("chrX:461829-461900(-)".to_owned()),
        );
        test_contig_ixn(
            "chrX:461829-462426(+)",
            "chrX:461800-461900(-)",
            Some("chrX:461829-461900(+)".to_owned()),
        );

        test_contig_ixn(
            "chrX:461829-462426(+)",
            "chrX:462000-463000(+)",
            Some("chrX:462000-462426(+)".to_owned()),
        );
        test_contig_ixn(
            "chrX:461829-462426(+)",
            "chrX:461000-463000(+)",
            Some("chrX:461829-462426(+)".to_owned()),
        );
        test_contig_ixn(
            "chrX:461829-462426(+)",
            "chrX:462000-462100(+)",
            Some("chrX:462000-462100(+)".to_owned()),
        );

        test_contig_ixn("chrX:461829-462426(+)", "chrX:461000-461500(+)", None);
        test_contig_ixn("chrX:461829-462426(+)", "chrX:463000-463500(+)", None);
        test_contig_ixn("chrX:461829-462426(+)", "chrV:461000-463000(+)", None);
    }
}
