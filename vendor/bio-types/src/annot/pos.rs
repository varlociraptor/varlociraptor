// Copyright 2017 Nicholas Ingolia
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Positions on a named sequence, e.g., 683,946 on chromosome IV.

use std::convert::Into;
use std::fmt::{self, Display, Formatter};
use std::ops::AddAssign;
use std::ops::Neg;
use std::ops::SubAssign;
use std::str::FromStr;

use regex::Regex;

use crate::annot::contig::Contig;
use crate::annot::loc::Loc;
use crate::annot::*;
use crate::strand::*;

/// Position on a particular, named sequence (e.g. a chromosome).
///
/// Parameterized over the type of the reference sequence identifier
/// and over the strandedness of the position.
///
/// The display format for a `Pos` is _chr:pos(+/-)_. A stranded
/// position must have a _(+)_ or a _(-)_, while an unstranded
/// position does not.
///
/// ```
/// # use bio_types::annot::ParseAnnotError;
/// # fn try_main() -> Result<(), Box<ParseAnnotError>> {
/// use bio_types::annot::pos::Pos;
/// use bio_types::strand::ReqStrand;
/// let start = Pos::new("chrIV".to_owned(), 683946, ReqStrand::Reverse);
/// let start_str = start.to_string();
/// assert_eq!(start_str, "chrIV:683946(-)");
/// let start_str_pos = start_str.parse()?;
/// assert_eq!(start, start_str_pos);
/// # Ok(())
/// # }
/// # fn main() { try_main().unwrap(); }
/// ```
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct Pos<R, S> {
    refid: R,
    pos: isize,
    strand: S,
}

impl<R, S> Pos<R, S> {
    /// Construct a new sequence position
    ///
    /// ```
    /// use std::rc::Rc;
    /// use bio_types::annot::pos::Pos;
    /// use bio_types::strand::ReqStrand;
    /// let chr = Rc::new("chrIV".to_owned());
    /// let start = Pos::new(chr, 683946, ReqStrand::Reverse);
    /// ```
    pub fn new(refid: R, pos: isize, strand: S) -> Self {
        Pos { refid, pos, strand }
    }

    /// Position on the reference sequence (0-based).
    pub fn pos(&self) -> isize {
        self.pos
    }

    /// Convert into a stranded sequence position on the specified strand
    pub fn into_stranded(self, strand: ReqStrand) -> Pos<R, ReqStrand> {
        Pos {
            refid: self.refid,
            pos: self.pos,
            strand,
        }
    }
}

impl<R, T> AddAssign<T> for Pos<R, ReqStrand>
where
    isize: AddAssign<T>,
    isize: SubAssign<T>,
{
    /// Slide the reference position by an offset on the strand of the
    /// annotation.
    ///
    /// # Arguments
    ///
    /// * `dist` specifies the offset for sliding the position. A
    /// positive `dist` will numerically increase the position for
    /// forward-strand features and decrease it for reverse-strand
    /// features.
    ///
    /// ```
    /// use bio_types::annot::pos::Pos;
    /// use bio_types::strand::ReqStrand;
    /// let mut start = Pos::new("chrIV".to_owned(), 683946, ReqStrand::Reverse);
    /// assert_eq!(start.to_string(), "chrIV:683946(-)");
    /// start += 100;
    /// assert_eq!(start.to_string(), "chrIV:683846(-)");
    /// ```
    fn add_assign(&mut self, dist: T) {
        match self.strand {
            ReqStrand::Forward => self.pos += dist,
            ReqStrand::Reverse => self.pos -= dist,
        }
    }
}

impl<R, T> SubAssign<T> for Pos<R, ReqStrand>
where
    isize: AddAssign<T>,
    isize: SubAssign<T>,
{
    /// Slide the reference position by an offset on the strand of the
    /// annotation.
    ///
    /// # Arguments
    ///
    /// * `dist` specifies the offset for sliding the position. A
    /// positive `dist` will numerically decrease the position for
    /// forward-strand features and increase it for reverse-strand
    /// features.
    ///
    /// ```
    /// use bio_types::annot::pos::Pos;
    /// use bio_types::strand::ReqStrand;
    /// let mut start = Pos::new("chrIV".to_owned(), 683946, ReqStrand::Reverse);
    /// assert_eq!(start.to_string(), "chrIV:683946(-)");
    /// start -= 100;
    /// assert_eq!(start.to_string(), "chrIV:684046(-)");
    /// ```
    fn sub_assign(&mut self, dist: T) {
        match self.strand {
            ReqStrand::Forward => self.pos -= dist,
            ReqStrand::Reverse => self.pos += dist,
        }
    }
}

impl<R, S> Loc for Pos<R, S> {
    type RefID = R;
    type Strand = S;
    fn refid(&self) -> &R {
        &self.refid
    }
    fn start(&self) -> isize {
        self.pos
    }
    fn length(&self) -> usize {
        1
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
        if (self.refid != pos.refid) || (self.pos != pos.pos) {
            None
        } else {
            Some(Pos::new(
                (),
                0,
                self.strand().into().on_strand(pos.strand()),
            ))
        }
    }

    fn pos_outof<Q, T>(&self, pos: &Pos<Q, T>) -> Option<Pos<Self::RefID, T>>
    where
        Self::RefID: Clone,
        Self::Strand: Into<ReqStrand> + Copy,
        T: Neg<Output = T> + Copy,
    {
        if pos.pos == 0 {
            Some(Pos::new(
                self.refid.clone(),
                self.pos,
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

        if (self.pos >= contig.start()) && (self.pos < (contig.start() + contig.length() as isize))
        {
            Some(self.clone())
        } else {
            None
        }
    }
}

impl<R, S> Same for Pos<R, S>
where
    R: Eq,
    S: Same,
{
    /// Indicate when two positions are the "same" -- when positions
    /// have unknown/unspecified strands they can be the "same" but
    /// not equal.
    fn same(&self, p: &Self) -> bool {
        self.pos == p.pos && self.refid == p.refid && self.strand.same(&p.strand)
    }
}

impl<R, S> Display for Pos<R, S>
where
    R: Display,
    S: Display + Clone + Into<Strand>,
{
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        let strand: Strand = self.strand.clone().into();
        if strand.is_unknown() {
            write!(f, "{}:{}", self.refid, self.pos)
        } else {
            write!(f, "{}:{}({})", self.refid, self.pos, strand)
        }
    }
}

impl<R, S> FromStr for Pos<R, S>
where
    R: From<String>,
    S: FromStr<Err = StrandError>,
{
    type Err = ParseAnnotError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        lazy_static! {
            static ref POS_RE: Regex = Regex::new(r"^(.*):(\d+)(\([+-]\))?$").unwrap();
        }

        let cap = POS_RE.captures(s).ok_or(ParseAnnotError::BadAnnot)?;

        let strand = cap
            .get(3)
            .map_or("", |m| m.as_str())
            .parse::<S>()
            .map_err(ParseAnnotError::ParseStrand)?;

        Ok(Pos::new(
            R::from(cap[1].to_owned()),
            cap[2].parse::<isize>().map_err(ParseAnnotError::ParseInt)?,
            strand,
        ))
    }
}

impl<R> From<Pos<R, ReqStrand>> for Pos<R, Strand> {
    fn from(x: Pos<R, ReqStrand>) -> Self {
        Pos {
            refid: x.refid,
            pos: x.pos,
            strand: match x.strand {
                ReqStrand::Forward => Strand::Forward,
                ReqStrand::Reverse => Strand::Reverse,
            },
        }
    }
}

impl<R> From<Pos<R, NoStrand>> for Pos<R, Strand> {
    fn from(x: Pos<R, NoStrand>) -> Self {
        Pos {
            refid: x.refid,
            pos: x.pos,
            strand: Strand::Unknown,
        }
    }
}

impl<R> From<Pos<R, Strand>> for Pos<R, NoStrand> {
    fn from(x: Pos<R, Strand>) -> Self {
        Pos {
            refid: x.refid,
            pos: x.pos,
            strand: NoStrand::Unknown,
        }
    }
}

impl<R> From<Pos<R, ReqStrand>> for Pos<R, NoStrand> {
    fn from(x: Pos<R, ReqStrand>) -> Self {
        Pos {
            refid: x.refid,
            pos: x.pos,
            strand: NoStrand::Unknown,
        }
    }
}

/// Default stranded sequence position on a reference sequence named
/// by a `String`.
pub type SeqPosStranded = Pos<String, ReqStrand>;

/// Default unstranded sequence position on a reference sequence named
/// by a `String`
pub type SeqPosUnstranded = Pos<String, NoStrand>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pos_accessors() {
        let start = Pos::new("chrIV".to_owned(), 683946, Strand::Unknown);
        assert_eq!(start.refid(), "chrIV");
        assert_eq!(start.pos(), 683946);
        assert!(start.strand().same(&Strand::Unknown));

        let start = Pos::new("chrIV".to_owned(), 683946, Strand::Reverse);
        assert_eq!(start.refid(), "chrIV");
        assert_eq!(start.pos(), 683946);
        assert!(start.strand().same(&Strand::Reverse));

        let start = Pos::new("chrXV".to_owned(), 493433, Strand::Forward);
        assert_eq!(start.refid(), "chrXV");
        assert_eq!(start.pos(), 493433);
        assert!(start.strand().same(&Strand::Forward));
    }

    #[test]
    fn strand_conversion() {
        let start = "chrIV:683946(-)".parse::<Pos<String, Strand>>().unwrap();
        let start_un: Pos<String, NoStrand> = start.into();
        assert!(start_un.same(&"chrIV:683946".parse::<Pos<String, NoStrand>>().unwrap()));
        let start_re = start_un.into_stranded(ReqStrand::Reverse);
        assert!(start_re.same(&"chrIV:683946(-)".parse::<Pos<String, ReqStrand>>().unwrap()));

        let start = "chrXV:493433(+)".parse::<Pos<String, Strand>>().unwrap();
        let start_un: Pos<String, NoStrand> = start.into();
        assert!(start_un.same(&"chrXV:493433".parse::<Pos<String, NoStrand>>().unwrap()));
        let start_re = start_un.into_stranded(ReqStrand::Forward);
        assert!(start_re.same(&"chrXV:493433(+)".parse::<Pos<String, ReqStrand>>().unwrap()));
    }

    #[test]
    fn string_representation() {
        let start = Pos::new("chrIV".to_owned(), 683946, NoStrand::Unknown);
        assert_eq!(start.to_string(), "chrIV:683946");
        assert!(start.same(&"chrIV:683946".parse::<Pos<String, NoStrand>>().unwrap()));

        let start = Pos::new("chrIV".to_owned(), 683946, Strand::Unknown);
        assert_eq!(start.to_string(), "chrIV:683946");
        assert!(start.same(&"chrIV:683946".parse::<Pos<String, Strand>>().unwrap()));

        let start = Pos::new("chrIV".to_owned(), 683946, Strand::Reverse);
        assert_eq!(start.to_string(), "chrIV:683946(-)");
        assert!(start.same(&"chrIV:683946(-)".parse::<Pos<String, Strand>>().unwrap()));

        let start = Pos::new("chrXV".to_owned(), 493433, Strand::Forward);
        assert_eq!(start.to_string(), "chrXV:493433(+)");
        assert!(start.same(&"chrXV:493433(+)".parse::<Pos<String, Strand>>().unwrap()));

        let start = Pos::new("chrIV".to_owned(), 683946, ReqStrand::Reverse);
        assert_eq!(start.to_string(), "chrIV:683946(-)");
        assert!(start.same(&"chrIV:683946(-)".parse::<Pos<String, ReqStrand>>().unwrap()));

        let start = Pos::new("chrXV".to_owned(), 493433, ReqStrand::Forward);
        assert_eq!(start.to_string(), "chrXV:493433(+)");
        assert!(start.same(&"chrXV:493433(+)".parse::<Pos<String, ReqStrand>>().unwrap()));
    }

    #[test]
    fn loc_impl() {
        let start = Pos::new("chrIV".to_owned(), 683946, ReqStrand::Forward);

        assert_eq!(
            None,
            start.contig_intersection(&Contig::new(
                "chrIV".to_owned(),
                683900,
                40,
                ReqStrand::Forward
            ))
        );
        assert_eq!(
            None,
            start.contig_intersection(&Contig::new(
                "chrV".to_owned(),
                683900,
                100,
                ReqStrand::Forward
            ))
        );
        assert_eq!(
            None,
            start.contig_intersection(&Contig::new(
                "chrIV".to_owned(),
                683950,
                40,
                ReqStrand::Forward
            ))
        );

        assert_eq!(
            Some(start.clone()),
            start.contig_intersection(&Contig::new(
                "chrIV".to_owned(),
                683900,
                100,
                ReqStrand::Forward
            ))
        );
        assert_eq!(
            Some(start.clone()),
            start.contig_intersection(&Contig::new(
                "chrIV".to_owned(),
                683900,
                100,
                ReqStrand::Reverse
            ))
        );

        let rstart = Pos::new("chrIV".to_owned(), 683946, ReqStrand::Reverse);
        assert_eq!(
            Some(rstart.clone()),
            rstart.contig_intersection(&Contig::new(
                "chrIV".to_owned(),
                683900,
                100,
                ReqStrand::Forward
            ))
        );
        assert_eq!(
            Some(rstart.clone()),
            rstart.contig_intersection(&Contig::new(
                "chrIV".to_owned(),
                683900,
                100,
                ReqStrand::Reverse
            ))
        );
    }
}
// chrXV:493433..494470
