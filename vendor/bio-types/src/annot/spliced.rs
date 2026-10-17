// Copyright 2017 Nicholas Ingolia
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Spliced region on a named sequence, e.g., the reverse strand of
//! chromosome V, exon #1 at 166,885 through 166,875 and exon #2 at
//! 166,771 through 166,237.

use std::cmp::{max, min};
use std::convert::Into;
use std::fmt::{self, Display, Formatter};
use std::ops::Neg;
use std::str::FromStr;

use regex::Regex;

use crate::annot::contig::Contig;
use crate::annot::loc::Loc;
use crate::annot::pos::Pos;
use crate::annot::*;
use crate::strand::*;

// The spliced location representation inherently cannot represent
// "bad" splicing structures. Locations comprise a first exon length,
// with a vector of intron-length / exon-length pairs, all type usize
// and hence non-negative.
//
// The InEx data type used to represent intron-exon pairs after the
// first exon enforces an additional strict positivity constraint
// (i.e., >0) on the lengths. This further eliminates degeneracy in
// representations: only one unique set of positive-length exons with
// interleaved positive-length introns can represent a splicing
// structure.
mod inex {
    use std::slice::Iter;

    use super::SplicingError;

    #[derive(Debug, Clone, Hash, PartialEq, Eq)]
    pub struct InEx {
        intron_length: usize,
        exon_length: usize,
    }

    impl InEx {
        pub fn new(intron_length: usize, exon_length: usize) -> Result<Self, SplicingError> {
            if intron_length < 1 {
                Err(SplicingError::IntronLength)
            } else if exon_length < 1 {
                Err(SplicingError::ExonLength)
            } else {
                Ok(InEx {
                    intron_length,
                    exon_length,
                })
            }
        }
        pub fn intron_length(&self) -> usize {
            self.intron_length
        }
        pub fn exon_length(&self) -> usize {
            self.exon_length
        }
        pub fn length(&self) -> usize {
            self.intron_length + self.exon_length
        }
    }

    // Represent just the start (relative to the start of the location) and length of exons
    // Useful for internal coordinate math
    #[derive(Debug, Clone, Hash, PartialEq, Eq)]
    pub struct Ex {
        start: usize,
        length: usize,
    }

    impl Ex {
        pub fn start(&self) -> usize {
            self.start
        }
        pub fn length(&self) -> usize {
            self.length
        }
        pub fn end(&self) -> usize {
            self.start + self.length
        }
    }

    // Iterator over the Ex exons from lowest to highest coordinate
    pub struct Exes<'a> {
        state: ExesState,
        curr_start: usize,
        rest: Iter<'a, InEx>,
    }

    enum ExesState {
        FirstExon(usize),
        LaterExon,
    }

    impl<'a> Exes<'a> {
        pub fn new(exon_0_length: usize, inexes: &'a [InEx]) -> Self {
            Exes {
                state: ExesState::FirstExon(exon_0_length),
                curr_start: 0,
                rest: inexes.iter(),
            }
        }
    }

    impl<'a> Iterator for Exes<'a> {
        type Item = Ex;

        fn next(&mut self) -> Option<Ex> {
            match self.state {
                ExesState::FirstExon(len) => {
                    let ex = Ex {
                        start: self.curr_start,
                        length: len,
                    };
                    self.curr_start += len;
                    self.state = ExesState::LaterExon;
                    Some(ex)
                }
                ExesState::LaterExon => match self.rest.next() {
                    Some(inex) => {
                        let ex = Ex {
                            start: self.curr_start + inex.intron_length(),
                            length: inex.exon_length(),
                        };
                        self.curr_start += inex.length();
                        Some(ex)
                    }
                    None => None,
                },
            }
        }
    }

    // Represent just the start (relative to the start of the location) and length of introns
    // Useful for internal coordinate math
    #[derive(Debug, Clone, Hash, PartialEq, Eq)]
    pub struct In {
        start: usize,
        length: usize,
    }

    #[allow(dead_code)]
    impl In {
        pub fn start(&self) -> usize {
            self.start
        }
        pub fn length(&self) -> usize {
            self.length
        }
        pub fn end(&self) -> usize {
            self.start + self.length
        }
    }

    // Iterator over the Ex introns from lowest to highest coordinate
    pub struct Ins<'a> {
        curr_start: usize,
        rest: Iter<'a, InEx>,
    }

    impl<'a> Ins<'a> {
        #[allow(dead_code)]
        pub fn new(exon_0_length: usize, inexes: &'a [InEx]) -> Self {
            Ins {
                curr_start: exon_0_length,
                rest: inexes.iter(),
            }
        }
    }

    impl<'a> Iterator for Ins<'a> {
        type Item = In;

        fn next(&mut self) -> Option<In> {
            match self.rest.next() {
                Some(inex) => {
                    let intr = In {
                        start: self.curr_start,
                        length: inex.intron_length(),
                    };
                    self.curr_start += inex.length();
                    Some(intr)
                }
                None => None,
            }
        }
    }
}

/// Spliced sequence annotation on a particular, named sequence
/// (e.g. a chromosome).
///
/// Parameterized over the type of the reference sequence identifier
/// and over the strandedness of the position.
///
/// The display format for a `Spliced` is
/// _chr:start_0-end_0;start_1-end_1;...;start_N-end_N(+/-/.)_. The
/// boundaries for each individual exon are given as a half-open
/// 0-based interval, like the Rust `Range` and BED format.
///
/// ```
/// # use bio_types::strand::ReqStrand;
/// # use bio_types::annot::AnnotError;
/// # use bio_types::annot::spliced::{Spliced,SplicingError};
/// # fn try_main() -> Result<(), Box<SplicingError>> {
/// let tad3 = Spliced::with_lengths_starts("chrXII".to_owned(), 765265,
///                                         &vec![808,52,109], &vec![0,864,984],
///                                         ReqStrand::Reverse)?;
/// assert_eq!(tad3.to_string(), "chrXII:765265-766073;766129-766181;766249-766358(-)");
/// let tad3_exons = tad3.exon_contigs();
/// assert_eq!(tad3_exons.len(), 3);
/// assert_eq!(tad3_exons[0].to_string(), "chrXII:766249-766358(-)");
/// assert_eq!(tad3_exons[1].to_string(), "chrXII:766129-766181(-)");
/// assert_eq!(tad3_exons[2].to_string(), "chrXII:765265-766073(-)");
/// # Ok(())
/// # }
/// # fn main() { try_main().unwrap(); }
/// ```
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct Spliced<R, S> {
    refid: R,
    start: isize,
    exon_0_length: usize,
    inexes: Vec<inex::InEx>,
    strand: S,
}

impl<R, S> Spliced<R, S> {
    /// Construct a new, single-exon "spliced" location
    ///
    /// ```
    /// use std::rc::Rc;
    /// use bio_types::annot::spliced::Spliced;
    /// use bio_types::strand::ReqStrand;
    /// let chr = Rc::new("chrX".to_owned());
    /// let tma22 = Spliced::new(chr, 461829, 462426 - 461829, ReqStrand::Forward);
    /// ```
    pub fn new(refid: R, start: isize, exon_0_length: usize, strand: S) -> Self {
        Spliced {
            refid,
            start,
            exon_0_length,
            inexes: Vec::new(),
            strand,
        }
    }

    /// Construct a multi-exon "spliced" location using BED-style exon
    /// starts and lengths.
    pub fn with_lengths_starts(
        refid: R,
        start: isize,
        exon_lengths: &[usize],
        exon_starts: &[usize],
        strand: S,
    ) -> Result<Self, SplicingError> {
        if exon_starts.is_empty() {
            return Err(SplicingError::NoExons);
        } else if exon_starts[0] != 0 {
            return Err(SplicingError::BlockStart);
        } else if exon_starts.len() != exon_lengths.len() {
            return Err(SplicingError::BlockMismatch);
        }

        let exon_0_length = exon_lengths[0];
        let mut intron_start = exon_0_length;
        let mut inexes = Vec::new();
        for exno in 1..exon_starts.len() {
            let exon_start = exon_starts[exno];
            if intron_start >= exon_start {
                return Err(SplicingError::BlockOverlap);
            }
            let intron_length = exon_start - intron_start;
            let exon_length = exon_lengths[exno];
            if exon_length == 0 {
                return Err(SplicingError::ExonLength);
            }
            inexes.push(inex::InEx::new(intron_length, exon_length)?);
            intron_start = exon_start + exon_length;
        }

        Ok(Spliced {
            refid,
            start,
            exon_0_length,
            inexes,
            strand,
        })
    }

    /// Number of exons
    pub fn exon_count(&self) -> usize {
        self.inexes.len() + 1
    }

    // This is the unique representation for zero-length Spliced
    // locations, because InEx pairs have positive lengths for both
    // introns and exons.
    #[allow(dead_code)]
    fn is_zero_length(&self) -> bool {
        self.exon_0_length == 0 && self.inexes.is_empty()
    }

    fn exes(&self) -> inex::Exes {
        inex::Exes::new(self.exon_0_length, &self.inexes)
    }

    /// Vector of exon starting positions, relative to the start of
    /// the location overall.
    ///
    /// These positions run from left to right on the reference
    /// sequence, regardless of the location's strand.
    pub fn exon_starts(&self) -> Vec<usize> {
        let mut starts = vec![0];
        let mut intron_start = self.exon_0_length;
        for inex in self.inexes.iter() {
            starts.push(intron_start + inex.intron_length());
            intron_start += inex.length();
        }
        starts
    }

    /// Vector of exon lengths.
    ///
    /// Exon lengths are given from left to right on the reference
    /// sequence, regardless of the location's strand.
    pub fn exon_lengths(&self) -> Vec<usize> {
        let mut lengths = vec![self.exon_0_length];
        for inex in self.inexes.iter() {
            lengths.push(inex.exon_length());
        }
        lengths
    }

    /// Total length of exons only.
    ///
    /// The `length` method from the `Loc` trait returns the total
    /// length spanned by the annotation, including both introns and
    /// exons.
    pub fn exon_total_length(&self) -> usize {
        self.exes().map(|e| e.length()).sum()
    }

    /// Convert into a stranded sequence location on the specified strand
    pub fn into_stranded(self, strand: ReqStrand) -> Spliced<R, ReqStrand> {
        Spliced {
            refid: self.refid,
            start: self.start,
            exon_0_length: self.exon_0_length,
            inexes: self.inexes,
            strand,
        }
    }

    pub fn contig_cover(self) -> Contig<R, S> {
        let length = self.length();
        Contig::new(self.refid, self.start, length, self.strand)
    }

    // Get exons in reference sequence order.
    fn exon_contigs_vec(&self) -> Vec<Contig<R, S>>
    where
        R: Clone,
        S: Copy,
    {
        let mut exons = Vec::new();

        for ex in self.exes() {
            exons.push(Contig::new(
                self.refid().clone(),
                self.start + ex.start() as isize,
                ex.length(),
                self.strand,
            ));
        }

        exons
    }
}

impl<R> Spliced<R, ReqStrand> {
    pub fn exon_contigs(&self) -> Vec<Contig<R, ReqStrand>>
    where
        R: Clone,
    {
        let mut exons = self.exon_contigs_vec();
        if self.strand == ReqStrand::Reverse {
            exons.reverse()
        }
        exons
    }
}

impl<R, S> Loc for Spliced<R, S> {
    type RefID = R;
    type Strand = S;
    fn refid(&self) -> &R {
        &self.refid
    }
    fn start(&self) -> isize {
        self.start
    }
    fn length(&self) -> usize {
        let mut len = self.exon_0_length;
        for inex in self.inexes.iter() {
            len += inex.length()
        }
        len
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
        if (self.refid != *pos.refid()) || pos.pos() < self.start {
            return None;
        }

        let pos_offset = (pos.pos() - self.start) as usize;

        let mut offset_before = 0;
        for ex in self.exes() {
            if pos_offset >= ex.start() && pos_offset < ex.end() {
                let offset = (offset_before + pos_offset - ex.start()) as isize;
                let into = match self.strand().into() {
                    ReqStrand::Forward => Pos::new((), offset, pos.strand()),
                    ReqStrand::Reverse => Pos::new(
                        (),
                        self.exon_total_length() as isize - (offset + 1),
                        -pos.strand(),
                    ),
                };
                return Some(into);
            }
            offset_before += ex.length();
        }

        None
    }

    fn pos_outof<Q, T>(&self, pos: &Pos<Q, T>) -> Option<Pos<Self::RefID, T>>
    where
        Self::RefID: Clone,
        Self::Strand: Into<ReqStrand> + Copy,
        T: Neg<Output = T> + Copy,
    {
        let mut offset = match self.strand().into() {
            ReqStrand::Forward => pos.pos(),
            ReqStrand::Reverse => self.exon_total_length() as isize - (pos.pos() + 1),
        };

        if offset < 0 {
            return None;
        }

        for ex in self.exes() {
            if offset < ex.length() as isize {
                return Some(Pos::new(
                    self.refid.clone(),
                    self.start + ex.start() as isize + offset,
                    self.strand.into().on_strand(pos.strand()),
                ));
            }
            offset -= ex.length() as isize;
        }

        None
    }

    fn contig_intersection<T>(&self, contig: &Contig<Self::RefID, T>) -> Option<Self>
    where
        Self::RefID: PartialEq + Clone,
        Self::Strand: Copy,
    {
        if self.refid() != contig.refid() {
            return None;
        }

        let contig_rel_start = if contig.start() < self.start {
            0
        } else {
            (contig.start() - self.start) as usize
        };
        let contig_end = contig.start() + contig.length() as isize;
        let contig_rel_end = if contig_end < self.start {
            0
        } else {
            (contig_end - self.start) as usize
        };

        let mut exon_lengths = Vec::new();
        let mut exon_starts = Vec::new();

        for ex in self.exes() {
            let start = max(contig_rel_start, ex.start());
            let end = min(contig_rel_end, ex.end());

            if start < end {
                exon_starts.push(start - contig_rel_start);
                exon_lengths.push(end - start);
            }
        }

        if !exon_starts.is_empty() {
            let first_start = exon_starts[0];
            for start in exon_starts.iter_mut() {
                *start -= first_start;
            }
            let ixn = Self::with_lengths_starts(
                self.refid.clone(),
                max(self.start, contig.start()) + first_start as isize,
                &exon_lengths,
                &exon_starts,
                self.strand,
            )
            .unwrap_or_else(|e| {
                panic!(
                    "Creating intersection spliced: {:?} for {:?} {:?}",
                    e, exon_lengths, exon_starts
                )
            });

            Some(ixn)
        } else {
            None
        }
    }
}

impl<R, S> Display for Spliced<R, S>
where
    R: Display,
    S: Display + Clone + Into<Strand>,
{
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{}:", self.refid)?;

        let mut sep = false;
        for ex in self.exes() {
            let ex_start = self.start + ex.start() as isize;
            write!(
                f,
                "{}{}-{}",
                if sep { ";" } else { "" },
                ex_start,
                ex_start + ex.length() as isize
            )?;
            sep = true;
        }
        let strand: Strand = self.strand.clone().into();
        if !strand.is_unknown() {
            write!(f, "({})", self.strand)?;
        }
        Ok(())
    }
}

impl<R, S> FromStr for Spliced<R, S>
where
    R: From<String>,
    S: FromStr<Err = StrandError>,
{
    type Err = ParseAnnotError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        lazy_static! {
            static ref SPLICED_RE: Regex =
                Regex::new(r"^(.*):(\d+)-(\d+)((?:;\d+-\d+)*)(\([+-]\))?$").unwrap();
            static ref EXON_RE: Regex = Regex::new(r";(\d+)-(\d+)").unwrap();
        }

        let cap = SPLICED_RE.captures(s).ok_or(ParseAnnotError::BadAnnot)?;

        let mut starts = Vec::new();
        let mut lengths = Vec::new();

        let first_start = cap[2].parse::<isize>().map_err(ParseAnnotError::ParseInt)?;
        let first_end = cap[3].parse::<isize>().map_err(ParseAnnotError::ParseInt)?;
        let strand = cap[5].parse::<S>().map_err(ParseAnnotError::ParseStrand)?;

        starts.push(0);
        lengths.push((first_end - first_start) as usize);

        let exon_caps = EXON_RE.captures_iter(&cap[4]);

        for exon_cap in exon_caps {
            let next_start = exon_cap[1]
                .parse::<isize>()
                .map_err(ParseAnnotError::ParseInt)?;
            let next_end = exon_cap[2]
                .parse::<isize>()
                .map_err(ParseAnnotError::ParseInt)?;
            starts.push((next_start - first_start) as usize);
            lengths.push((next_end - next_start) as usize);
        }

        let spliced = Spliced::with_lengths_starts(
            R::from(cap[1].to_owned()),
            first_start,
            lengths.as_slice(),
            starts.as_slice(),
            strand,
        )
        .map_err(ParseAnnotError::Splicing)?;
        Ok(spliced)
    }
}

impl<R> From<Spliced<R, ReqStrand>> for Spliced<R, Strand> {
    fn from(x: Spliced<R, ReqStrand>) -> Self {
        Spliced {
            refid: x.refid,
            start: x.start,
            exon_0_length: x.exon_0_length,
            inexes: x.inexes,
            strand: match x.strand {
                ReqStrand::Forward => Strand::Forward,
                ReqStrand::Reverse => Strand::Reverse,
            },
        }
    }
}

impl<R> From<Spliced<R, NoStrand>> for Spliced<R, Strand> {
    fn from(x: Spliced<R, NoStrand>) -> Self {
        Spliced {
            refid: x.refid,
            start: x.start,
            exon_0_length: x.exon_0_length,
            inexes: x.inexes,
            strand: Strand::Unknown,
        }
    }
}

impl<R> From<Spliced<R, Strand>> for Spliced<R, NoStrand> {
    fn from(x: Spliced<R, Strand>) -> Self {
        Spliced {
            refid: x.refid,
            start: x.start,
            exon_0_length: x.exon_0_length,
            inexes: x.inexes,
            strand: NoStrand::Unknown,
        }
    }
}

impl<R> From<Spliced<R, ReqStrand>> for Spliced<R, NoStrand> {
    fn from(x: Spliced<R, ReqStrand>) -> Self {
        Spliced {
            refid: x.refid,
            start: x.start,
            exon_0_length: x.exon_0_length,
            inexes: x.inexes,
            strand: NoStrand::Unknown,
        }
    }
}

/// Default stranded sequence position on a reference sequence named
/// by a `String`.
pub type SeqSplicedStranded = Spliced<String, ReqStrand>;

/// Default unstranded sequence position on a reference sequence named
/// by a `String`
pub type SeqSplicedUnstranded = Spliced<String, NoStrand>;

#[derive(Error, Debug)]
pub enum SplicingError {
    #[error("Invalid (non-positive) intron length")]
    IntronLength,
    #[error("Invalid (non-positive) exon length")]
    ExonLength,
    #[error("No exons")]
    NoExons,
    #[error("Exons do not start at position 0")]
    BlockStart,
    #[error("Number of exon starts != number of exon lengths")]
    BlockMismatch,
    #[error("Exon blocks overlap")]
    BlockOverlap,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn length_start_to_contig() {
        //chrV    166236  166885  YER007C-A       0       -       166236  166885  0       2       535,11, 0,638,
        let tma20 = Spliced::with_lengths_starts(
            "chrV".to_owned(),
            166236,
            &vec![535, 11],
            &vec![0, 638],
            ReqStrand::Reverse,
        )
        .unwrap();
        assert_eq!(tma20.exon_starts(), vec![0, 638]);
        assert_eq!(tma20.exon_lengths(), vec![535, 11]);
        assert_eq!(tma20.to_string(), "chrV:166236-166771;166874-166885(-)");
        assert_eq!(
            tma20,
            tma20
                .to_string()
                .parse::<Spliced<String, ReqStrand>>()
                .unwrap()
        );
        let tma20_exons = tma20.exon_contigs();
        assert_eq!(tma20_exons.len(), 2);
        assert_eq!(tma20_exons[0].to_string(), "chrV:166874-166885(-)");
        assert_eq!(tma20_exons[1].to_string(), "chrV:166236-166771(-)");

        //chrXVI  173151  174702  YPL198W 0       +       173151  174702  0       3       11,94,630,      0,420,921,
        let rpl7b = Spliced::with_lengths_starts(
            "chrXVI".to_owned(),
            173151,
            &vec![11, 94, 630],
            &vec![0, 420, 921],
            ReqStrand::Forward,
        )
        .unwrap();
        assert_eq!(
            rpl7b.to_string(),
            "chrXVI:173151-173162;173571-173665;174072-174702(+)"
        );
        assert_eq!(
            rpl7b,
            rpl7b
                .to_string()
                .parse::<Spliced<String, ReqStrand>>()
                .unwrap()
        );
        let rpl7b_exons = rpl7b.exon_contigs();
        assert_eq!(rpl7b_exons.len(), 3);
        assert_eq!(rpl7b_exons[0].to_string(), "chrXVI:173151-173162(+)");
        assert_eq!(rpl7b_exons[1].to_string(), "chrXVI:173571-173665(+)");
        assert_eq!(rpl7b_exons[2].to_string(), "chrXVI:174072-174702(+)");

        //chrXII  765265  766358  YLR316C 0       -       765265  766358  0       3       808,52,109,     0,864,984,
        let tad3 = Spliced::with_lengths_starts(
            "chrXII".to_owned(),
            765265,
            &vec![808, 52, 109],
            &vec![0, 864, 984],
            ReqStrand::Reverse,
        )
        .unwrap();
        assert_eq!(
            tad3.to_string(),
            "chrXII:765265-766073;766129-766181;766249-766358(-)"
        );
        assert_eq!(
            tad3,
            tad3.to_string()
                .parse::<Spliced<String, ReqStrand>>()
                .unwrap()
        );
        let tad3_exons = tad3.exon_contigs();
        assert_eq!(tad3_exons.len(), 3);
        assert_eq!(tad3_exons[0].to_string(), "chrXII:766249-766358(-)");
        assert_eq!(tad3_exons[1].to_string(), "chrXII:766129-766181(-)");
        assert_eq!(tad3_exons[2].to_string(), "chrXII:765265-766073(-)");
    }

    fn test_into_outof(
        loc: &Spliced<String, ReqStrand>,
        outstr: &str,
        in_offset: isize,
        in_strand: ReqStrand,
    ) -> () {
        let p0 = outstr.parse::<Pos<String, ReqStrand>>().unwrap();
        let p0_into_expected = Pos::new((), in_offset, in_strand);
        let p0_into_actual = loc.pos_into(&p0);
        let p0_back_out_actual = loc.pos_outof(&p0_into_expected);
        println!(
            "{}\t{}\t{:?}\t{:?}\t{:?}",
            outstr, p0, p0_into_expected, p0_into_actual, p0_back_out_actual
        );
        assert!(Some(p0_into_expected).same(&p0_into_actual));
        assert!(Some(p0).same(&p0_back_out_actual));
    }

    fn test_no_into(loc: &Spliced<String, ReqStrand>, outstr: &str) -> () {
        let p0 = outstr.parse::<Pos<String, ReqStrand>>().unwrap();
        assert!(None.same(&loc.pos_into(&p0)));
    }

    #[test]
    fn into_outof() {
        //chrXVI  173151  174702  YPL198W 0       +       173151  174702  0       3       11,94,630,      0,420,921,
        let rpl7b = Spliced::with_lengths_starts(
            "chrXVI".to_owned(),
            173151,
            &vec![11, 94, 630],
            &vec![0, 420, 921],
            ReqStrand::Forward,
        )
        .unwrap();
        let p0_into = Pos::new((), -1, ReqStrand::Forward);
        assert!(None.same(&rpl7b.pos_outof(&p0_into)));

        test_no_into(&rpl7b, "chrXVI:173150(+)");
        test_into_outof(&rpl7b, "chrXVI:173151(+)", 0, ReqStrand::Forward);
        test_into_outof(&rpl7b, "chrXVI:173152(-)", 1, ReqStrand::Reverse);
        test_into_outof(&rpl7b, "chrXVI:173161(+)", 10, ReqStrand::Forward);
        test_no_into(&rpl7b, "chrXVI:173162(+)");
        test_no_into(&rpl7b, "chrXVI:173570(+)");
        test_into_outof(&rpl7b, "chrXVI:173571(+)", 11, ReqStrand::Forward);
        test_into_outof(&rpl7b, "chrXVI:173664(+)", 104, ReqStrand::Forward);
        test_no_into(&rpl7b, "chrXVI:173665(+)");
        test_no_into(&rpl7b, "chrXVI:174071(+)");
        test_into_outof(&rpl7b, "chrXVI:174072(+)", 105, ReqStrand::Forward);
        test_into_outof(&rpl7b, "chrXVI:174701(+)", 734, ReqStrand::Forward);
        test_no_into(&rpl7b, "chrXVI:174702(+)");

        let p0_into = Pos::new((), 735, ReqStrand::Forward);
        assert!(None.same(&rpl7b.pos_outof(&p0_into)));

        //chrXII  765265  766358  YLR316C 0       -       765265  766358  0       3       808,52,109,     0,864,984,
        let tad3 = Spliced::with_lengths_starts(
            "chrXII".to_owned(),
            765265,
            &vec![808, 52, 109],
            &vec![0, 864, 984],
            ReqStrand::Reverse,
        )
        .unwrap();

        let p0_into = Pos::new((), -1, ReqStrand::Forward);
        assert!(None.same(&tad3.pos_outof(&p0_into)));

        test_no_into(&tad3, "chrXII:765264(-)");
        test_into_outof(&tad3, "chrXII:765265(-)", 968, ReqStrand::Forward);
        test_into_outof(&tad3, "chrXII:765266(+)", 967, ReqStrand::Reverse);
        test_into_outof(&tad3, "chrXII:766072(-)", 161, ReqStrand::Forward);
        test_no_into(&tad3, "chrXII:766073(-)");

        test_no_into(&tad3, "chrXII:766128(-)");
        test_into_outof(&tad3, "chrXII:766129(-)", 160, ReqStrand::Forward);
        test_into_outof(&tad3, "chrXII:766180(-)", 109, ReqStrand::Forward);
        test_no_into(&tad3, "chrXII:766181(-)");

        test_no_into(&tad3, "chrXII:766248(-)");
        test_into_outof(&tad3, "chrXII:766249(-)", 108, ReqStrand::Forward);
        test_into_outof(&tad3, "chrXII:766357(-)", 0, ReqStrand::Forward);
        test_no_into(&tad3, "chrXII:766358(-)");

        let p0_into = Pos::new((), 969, ReqStrand::Forward);
        assert!(None.same(&tad3.pos_outof(&p0_into)));
    }

    fn test_contig_ixn(
        spl: &Spliced<String, ReqStrand>,
        cb_str: &str,
        cab_str: Option<String>,
    ) -> () {
        let cb = cb_str.parse::<Contig<String, ReqStrand>>().unwrap();
        match spl.contig_intersection(&cb) {
            None => assert_eq!(None, cab_str),
            Some(cab) => assert_eq!(Some(cab.to_string()), cab_str),
        };
    }

    #[test]
    fn intersection() {
        //chrXVI  173151  174702  YPL198W 0       +       173151  174702  0       3       11,94,630,      0,420,921,
        let rpl7b = Spliced::with_lengths_starts(
            "chrXVI".to_owned(),
            173151,
            &vec![11, 94, 630],
            &vec![0, 420, 921],
            ReqStrand::Forward,
        )
        .unwrap();

        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-175000(+)",
            Some("chrXVI:173151-173162;173571-173665;174072-174702(+)".to_owned()),
        );

        test_contig_ixn(
            &rpl7b,
            "chrXVI:173150-175000(+)",
            Some("chrXVI:173151-173162;173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173151-175000(+)",
            Some("chrXVI:173151-173162;173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173152-175000(+)",
            Some("chrXVI:173152-173162;173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173155-175000(+)",
            Some("chrXVI:173155-173162;173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173161-175000(+)",
            Some("chrXVI:173161-173162;173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173162-175000(+)",
            Some("chrXVI:173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173500-175000(+)",
            Some("chrXVI:173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173570-175000(+)",
            Some("chrXVI:173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173571-175000(+)",
            Some("chrXVI:173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173572-175000(+)",
            Some("chrXVI:173572-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173600-175000(+)",
            Some("chrXVI:173600-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173664-175000(+)",
            Some("chrXVI:173664-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173665-175000(+)",
            Some("chrXVI:174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:174100-175000(+)",
            Some("chrXVI:174100-174702(+)".to_owned()),
        );
        test_contig_ixn(&rpl7b, "chrXVI:174800-175000(+)", None);

        test_contig_ixn(
            &rpl7b,
            "chrXVI:173150-174703(+)",
            Some("chrXVI:173151-173162;173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173150-174702(+)",
            Some("chrXVI:173151-173162;173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173150-174701(+)",
            Some("chrXVI:173151-173162;173571-173665;174072-174701(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-174500(+)",
            Some("chrXVI:173151-173162;173571-173665;174072-174500(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-174072(+)",
            Some("chrXVI:173151-173162;173571-173665(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-173800(+)",
            Some("chrXVI:173151-173162;173571-173665(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-173666(+)",
            Some("chrXVI:173151-173162;173571-173665(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-173665(+)",
            Some("chrXVI:173151-173162;173571-173665(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-173664(+)",
            Some("chrXVI:173151-173162;173571-173664(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-173600(+)",
            Some("chrXVI:173151-173162;173571-173600(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-173571(+)",
            Some("chrXVI:173151-173162(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-173300(+)",
            Some("chrXVI:173151-173162(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-173155(+)",
            Some("chrXVI:173151-173155(+)".to_owned()),
        );
        test_contig_ixn(&rpl7b, "chrXVI:173000-173100(+)", None);

        test_contig_ixn(
            &rpl7b,
            "chrXVI:173155-174500(+)",
            Some("chrXVI:173155-173162;173571-173665;174072-174500(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173600-174500(+)",
            Some("chrXVI:173600-173665;174072-174500(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173155-173600(+)",
            Some("chrXVI:173155-173162;173571-173600(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173590-173610(+)",
            Some("chrXVI:173590-173610(+)".to_owned()),
        );

        test_contig_ixn(
            &rpl7b,
            "chrXVI:173155-173160(+)",
            Some("chrXVI:173155-173160(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:174400-174500(+)",
            Some("chrXVI:174400-174500(+)".to_owned()),
        );

        test_contig_ixn(&rpl7b, "chrXVI:173200-173300(+)", None);
        test_contig_ixn(&rpl7b, "chrXVI:173800-174000(+)", None);
    }
}
