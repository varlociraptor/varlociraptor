// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cell::RefCell;
use std::cmp;
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::path::Path;
use std::rc::Rc;
use std::str;
use std::sync::Arc;

use anyhow::Result;
use bio::alphabets::dna;
use bio::stats::pairhmm::EmissionParameters;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval, AbstractLocus};
use regex::Regex;
use rust_htslib::bam;
use rust_htslib::bcf::{self, Read};
use vec_map::VecMap;

use crate::errors::Error;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::reference;
use crate::utils;
use crate::variants::evidence::realignment::pairhmm::{ReadEmission, RefBaseEmission};
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::model;
use crate::variants::sampling_bias::{ReadSamplingBias, SamplingBias};
use crate::variants::types::{
    AlleleSupport, MultiLocus, PairedEndEvidence, SingleLocus, SingleLocusBuilder, Variant,
};
use crate::{default_emission, default_ref_base_emission};

const MIN_REF_BASES: u64 = 10;

#[derive(Builder)]
#[builder(build_fn(name = "build_inner"))]
pub(crate) struct BreakendGroup {
    #[builder(default)]
    loci: MultiLocus,
    #[builder(private, default)]
    enclosable_ref_interval: Option<genome::Interval>,
    #[builder(default)]
    // TODO consider making the right side a Vec<Breakend>!
    breakends: BTreeMap<genome::Locus, Breakend>,
    #[builder(default)]
    alt_alleles: RefCell<VecMap<Vec<Arc<AltAllele>>>>,
    #[builder(private)]
    realigner: RefCell<Realigner>,
}

impl BreakendGroupBuilder {
    pub(crate) fn set_realigner(&mut self, realigner: Realigner) -> &mut Self {
        self.realigner = Some(RefCell::new(realigner));

        self
    }

    pub(crate) fn push_breakend(&mut self, breakend: Breakend) -> &mut Self {
        let interval = genome::Interval::new(
            breakend.locus.contig().to_owned(),
            breakend.locus.pos()..breakend.locus.pos() + breakend.ref_allele.len() as u64,
        );

        if self.breakends.is_none() {
            self.breakends = Some(BTreeMap::default());
            self.loci = Some(MultiLocus::default());
        }

        self.breakends
            .as_mut()
            .unwrap()
            .insert(breakend.locus.clone(), breakend);

        self.loci.as_mut().unwrap().push(
            SingleLocusBuilder::default()
                .interval(interval)
                .build()
                .unwrap(),
        );

        self
    }

    pub(crate) fn build(&mut self) -> Result<BreakendGroup, String> {
        // Calculate enclosable reference interval.
        let first = self.breakends.as_ref().unwrap().keys().next().unwrap();
        if self
            .breakends
            .as_ref()
            .unwrap()
            .values()
            .skip(1)
            .all(|bnd| bnd.locus.contig() == first.contig())
        {
            let interval = {
                let last = self.breakends.as_ref().unwrap().values().last().unwrap();
                genome::Interval::new(
                    last.locus.contig().to_owned(),
                    first.pos()
                        ..last.locus.pos()
                            + if !last.is_left_to_right() {
                                last.ref_allele.len() as u64
                            } else {
                                0
                            },
                )
            };
            self.enclosable_ref_interval(Some(interval));
        }

        self.build_inner()
    }
}

impl BreakendGroup {
    pub(crate) fn breakends(&self) -> impl Iterator<Item = &Breakend> {
        self.breakends.values()
    }

    fn upstream_bnd(&self, locus: &genome::Locus) -> Option<&Breakend> {
        for (l, bnd) in self.breakends.range(..locus).rev() {
            if l.contig() == locus.contig() {
                if l.pos() < locus.pos() && !bnd.is_left_to_right() {
                    // Return first locus with smaller position.
                    return Some(bnd);
                }
            } else {
                break;
            }
        }
        None
    }

    fn downstream_bnd(&self, locus: &genome::Locus) -> Option<&Breakend> {
        for (l, bnd) in self.breakends.range(locus..) {
            if l.contig() == locus.contig() {
                if l.pos() > locus.pos() && bnd.is_left_to_right() {
                    // Return first unvisited bnd with larger position.
                    return Some(bnd);
                }
            } else {
                break;
            }
        }
        None
    }

    fn breakend_pair(&self) -> Option<(&Breakend, &Breakend)> {
        if self.breakends.len() == 2 {
            let left = self.breakends.values().next().unwrap();
            let right = self.breakends.values().nth(1).unwrap();
            Some((left, right))
        } else {
            None
        }
    }

    fn is_insertion(&self) -> bool {
        if let Some((left, right)) = self.breakend_pair() {
            if left.locus.pos() + 1 == right.locus.pos()
                && left.locus.contig() == right.locus.contig()
                && !left.emits_revcomp()
                && !right.emits_revcomp()
                && left.is_left_to_right()
                && left.replacement.len() > 1
                && right.replacement[..right.replacement.len() - 1] == left.replacement[1..]
                && !right.is_left_to_right()
            {
                return true;
            }
        }
        false
    }

    fn is_deletion(&self) -> bool {
        if let Some((left, right)) = self.breakend_pair() {
            if left.locus.contig() == right.locus.contig()
                && left.replacement.len() == 1
                && right.replacement.len() == 1
                && left.is_left_to_right()
                && !right.is_left_to_right()
                && !left.emits_revcomp()
                && !right.emits_revcomp()
            {
                return true;
            }
        }
        false
    }
}

impl Variant for BreakendGroup {
    type Evidence = PairedEndEvidence;
    type Loci = MultiLocus;

    fn is_valid_evidence(&self, evidence: &Self::Evidence) -> Option<Vec<usize>> {
        let is_valid_overlap = |locus: &SingleLocus, read| !locus.overlap(read, true).is_none();

        let is_valid_ref_bases = |read: &bam::Record| {
            if let Some(ref interval) = self.enclosable_ref_interval {
                // METHOD: we try to ensure that at least MIN_REF_BASES bases are on the reference
                // By that, we avoid that a read is entirely inside the breakend. Such reads are
                // complicated to map against the alt allele (because e.g. with revcomp that might
                // at a completely different position on the genome). And further, they do not help
                // anyway, because they tend to map equally well to ref like to alt (e.g. as revcomp).
                cmp::max(
                    interval.range().start.saturating_sub(read.pos() as u64),
                    read.cigar_cached()
                        .unwrap()
                        .end_pos()
                        .saturating_sub(interval.range().end as i64) as u64,
                ) > MIN_REF_BASES
            } else {
                true
            }
        };

        let overlapping: Vec<_> = match evidence {
            PairedEndEvidence::SingleEnd(read) => {
                if !is_valid_ref_bases(read) {
                    return None;
                }
                self.loci
                    .iter()
                    .enumerate()
                    .filter_map(|(i, locus)| {
                        if is_valid_overlap(locus, read) {
                            Some(i)
                        } else {
                            None
                        }
                    })
                    .collect()
            }
            PairedEndEvidence::PairedEnd { left, right } => {
                if !is_valid_ref_bases(left) && !is_valid_ref_bases(right) {
                    return None;
                }

                self.loci
                    .iter()
                    .enumerate()
                    .filter_map(|(i, locus)| {
                        if is_valid_overlap(locus, left) || is_valid_overlap(locus, right) {
                            Some(i)
                        } else {
                            None
                        }
                    })
                    .collect()
            }
        };

        if overlapping.is_empty() {
            None
        } else {
            Some(overlapping)
        }
    }

    /// Return variant loci.
    fn loci(&self) -> &Self::Loci {
        &self.loci
    }

    /// Calculate probability for alt and reference allele.
    fn allele_support(
        &self,
        evidence: &Self::Evidence,
        _: &AlignmentProperties,
    ) -> Result<Option<AlleleSupport>> {
        match evidence {
            PairedEndEvidence::SingleEnd(record) => Ok(Some(
                self.realigner
                    .borrow_mut()
                    .allele_support(record, self.loci.iter(), self)?,
            )),
            PairedEndEvidence::PairedEnd { left, right } => {
                let left_support =
                    self.realigner
                        .borrow_mut()
                        .allele_support(left, self.loci.iter(), self)?;
                let right_support =
                    self.realigner
                        .borrow_mut()
                        .allele_support(right, self.loci.iter(), self)?;

                let mut support = left_support;

                support.merge(&right_support);

                Ok(Some(support))
            }
        }
    }

    /// Calculate probability to sample a record length like the given one from the alt allele.
    fn prob_sample_alt(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        match evidence {
            PairedEndEvidence::PairedEnd { left, right } => {
                // METHOD: we do not require the fragment to enclose the breakend group.
                // Hence, we treat both reads independently.
                (self
                    .prob_sample_alt_read(left.seq().len() as u64, alignment_properties)
                    .ln_one_minus_exp()
                    + self
                        .prob_sample_alt_read(right.seq().len() as u64, alignment_properties)
                        .ln_one_minus_exp())
                .ln_one_minus_exp()
            }
            PairedEndEvidence::SingleEnd(read) => {
                self.prob_sample_alt_read(read.seq().len() as u64, alignment_properties)
            }
        }
    }
}

impl SamplingBias for BreakendGroup {
    fn feasible_bases(&self, read_len: u64, alignment_properties: &AlignmentProperties) -> u64 {
        if self.is_deletion() {
            if let Some(len) = self.enclosable_len() {
                if len < (alignment_properties.max_del_cigar_len as u64) {
                    return read_len;
                }
            }
        } else if self.is_insertion() {
            if let Some(len) = self.enclosable_len() {
                if len < (alignment_properties.max_ins_cigar_len as u64) {
                    return read_len;
                }
            }
        }
        (read_len as f64 * alignment_properties.frac_max_softclip) as u64
    }

    fn enclosable_len(&self) -> Option<u64> {
        if self.is_deletion() {
            if let Some((left, right)) = self.breakend_pair() {
                return Some(right.locus.pos() - left.locus.pos());
            }
        } else if self.is_insertion() {
            return Some(self.breakends.values().next().unwrap().replacement.len() as u64 - 1);
        }
        None
    }
}

impl ReadSamplingBias for BreakendGroup {}

impl<'a> Realignable<'a> for BreakendGroup {
    type EmissionParams = BreakendEmissionParams<'a>;

    fn maybe_revcomp(&self) -> bool {
        self.breakends.values().any(|bnd| bnd.emits_revcomp())
    }

    fn alt_emission_params(
        &self,
        read_emission_params: Rc<ReadEmission<'a>>,
        ref_buffer: Arc<reference::Buffer>,
        _: &genome::Interval,
        ref_window: usize,
    ) -> Result<Vec<BreakendEmissionParams<'a>>> {
        // Step 1: fetch contained breakends
        let mut emission_params = Vec::new();

        // METHOD: we consider all breakends, even if they don't overlap.
        // The reason is that the mapper may put reads at the wrong end of a breakend pair.
        for (bnd_idx, first) in self.breakends.values().enumerate() {
            //for (bnd_idx, first) in self.contained_breakends(ref_interval) {
            if !self.alt_alleles.borrow().contains_key(bnd_idx) {
                let mut candidate_alt_alleles = Vec::new();
                // Compute alt allele sequence once.
                let mut alt_allele = AltAllele::default();

                let prefix_range = |bnd: &Breakend| {
                    let start = bnd.locus.pos() as usize;
                    start.saturating_sub(ref_window)..start
                };

                let suffix_range = |bnd: &Breakend, ref_len, inclusive| {
                    let start = bnd.locus.pos() as usize + if inclusive { 0 } else { 1 };
                    start..cmp::min(start + ref_window, ref_len)
                };

                // Decide whether to go from right to left or left to right
                if first.is_left_to_right() {
                    // Add prefix on reference.
                    alt_allele.push_seq(
                        ref_buffer.seq(first.locus.contig())?[prefix_range(first)].iter(),
                        false,
                        false,
                    );
                    //alt_allele.push_seq(b"1".iter(), false); // dbg
                    // Add replacement to alt allele.
                    alt_allele.push_seq(first.replacement().iter(), false, true);
                } else {
                    // Add suffix on reference.
                    let ref_seq = ref_buffer.seq(first.locus.contig())?;
                    alt_allele.push_seq(
                        ref_seq[suffix_range(first, ref_seq.len(), false)].iter(),
                        true,
                        false,
                    );
                    //alt_allele.push_seq(b"1".iter(), true); // dbg
                    // Add replacement to alt allele.
                    alt_allele.push_seq(first.replacement().iter(), true, true);
                }

                let mut revcomp = false;
                let mut next_bnd = Some(first);
                let mut visited = HashSet::new();
                while let Some(current) = next_bnd {
                    if visited.contains(&current.id) {
                        //alt_allele.push_seq(b"2".iter(), !current.is_left_to_right()); // dbg
                        // We are back at a previous breakend, time to stop.
                        let ref_seq = ref_buffer.seq(current.locus.contig())?;
                        if current.is_left_to_right() {
                            // Add suffix on reference.
                            alt_allele.push_seq(
                                ref_seq[suffix_range(current, ref_seq.len(), false)].iter(),
                                false,
                                false,
                            );
                        } else {
                            // Prepend prefix on reference.
                            alt_allele.push_seq(ref_seq[prefix_range(current)].iter(), true, false);
                        }
                        break;
                    }
                    visited.insert(&current.id);

                    let left_to_right = if revcomp {
                        !current.is_left_to_right()
                    } else {
                        current.is_left_to_right()
                    };

                    if let Some(ref join) = current.join() {
                        let ref_seq = ref_buffer.seq(join.locus.contig())?;
                        // Find next breakend from here in order to know how long to extend.
                        let seq = match join.side {
                            Side::LeftOfPos => {
                                next_bnd = self.upstream_bnd(&join.locus);
                                let seq_start =
                                    next_bnd.map_or(0, |bnd| bnd.locus.pos() as usize + 1); // end before the beginning of the next bnd
                                let seq_end = join.locus.pos() as usize + 1; // locus.pos() is meant inclusive, so we need to add 1 here.
                                &ref_seq[seq_start..seq_end]
                            }
                            Side::RightOfPos => {
                                next_bnd = self.downstream_bnd(&join.locus);
                                let seq_start = join.locus.pos() as usize; // locus.pos() is meant inclusive.
                                let seq_end =
                                    next_bnd.map_or(ref_seq.len(), |bnd| bnd.locus.pos() as usize); // end before the start of the next bnd
                                &ref_seq[seq_start..seq_end]
                            }
                        };

                        // Apply replacement operation of next bnd if necessary.
                        let seq: Vec<u8> = match (next_bnd, &join.side) {
                            (Some(next_bnd), Side::RightOfPos) => seq
                                .iter()
                                .chain(next_bnd.replacement().iter())
                                .cloned()
                                .collect(),
                            (Some(next_bnd), Side::LeftOfPos) => next_bnd
                                .replacement()
                                .iter()
                                .chain(seq.iter())
                                .cloned()
                                .collect(),
                            (None, _) => seq.to_owned(),
                        };

                        let extension_modification = if revcomp {
                            join.extension_modification.invert()
                        } else {
                            join.extension_modification
                        };

                        //alt_allele.push_seq(b"|".iter(), !left_to_right); // dbg

                        if next_bnd.is_some() && alt_allele.len() + seq.len() > ref_window {
                            // METHOD: sequence addition will already exceed ref window, hence
                            // we can stop afterwards.
                            next_bnd = None;
                        }

                        // Push sequence to alt allele.
                        match (extension_modification, next_bnd.is_none(), left_to_right) {
                            (ExtensionModification::None, false, is_left_to_right) => {
                                alt_allele.push_seq(seq.iter(), !is_left_to_right, true)
                            }
                            (ExtensionModification::ReverseComplement, false, is_left_to_right) => {
                                alt_allele.push_seq(
                                    dna::revcomp(&seq).iter(),
                                    !is_left_to_right,
                                    true,
                                );
                            }
                            (ExtensionModification::None, true, true) => {
                                alt_allele.push_seq(
                                    seq[..cmp::min(ref_window, seq.len())].iter(),
                                    false,
                                    true,
                                );
                            }
                            (ExtensionModification::ReverseComplement, true, true) => {
                                alt_allele.push_seq(
                                    dna::revcomp(seq)[..ref_window].iter(),
                                    false,
                                    true,
                                );
                            }
                            (ExtensionModification::None, true, false) => {
                                alt_allele.push_seq(
                                    seq[seq.len().saturating_sub(ref_window)..].iter(),
                                    true,
                                    true,
                                );
                            }
                            (ExtensionModification::ReverseComplement, true, false) => {
                                let seq = dna::revcomp(seq);
                                alt_allele.push_seq(
                                    seq[seq.len().saturating_sub(ref_window)..].iter(),
                                    true,
                                    true,
                                );
                            }
                        }

                        // Update revcomp marker for next iteration.
                        revcomp = if let ExtensionModification::ReverseComplement =
                            extension_modification
                        {
                            true
                        } else {
                            false
                        };
                    } else {
                        // Single breakend, assembly stops here.
                        // Nothing else to do, the replacement sequence has already been added in the step before.
                    }
                }

                let alt_allele = Arc::new(alt_allele);

                candidate_alt_alleles.push(alt_allele);

                self.alt_alleles
                    .borrow_mut()
                    .insert(bnd_idx, candidate_alt_alleles);
            }

            for alt_allele in self.alt_alleles.borrow().get(bnd_idx).unwrap() {
                emission_params.push(BreakendEmissionParams {
                    ref_offset: 0,
                    ref_end: alt_allele.len(),
                    alt_allele: Arc::clone(alt_allele),
                    read_emission: Rc::clone(&read_emission_params),
                });
            }
        }

        //dbg!(emission_params.iter().map(|p| std::str::from_utf8(&p.alt_allele.iter().cloned().collect::<Vec<u8>>()).unwrap().to_owned()).collect::<Vec<_>>());

        Ok(emission_params)
    }
}

#[derive(Derefable, Debug, Default)]
pub(crate) struct AltAllele {
    #[deref]
    seq: VecDeque<u8>,
    alt_len: u64,
}

impl AltAllele {
    pub(crate) fn push_seq<'a, S>(&mut self, seq: S, front: bool, is_alt: bool)
    where
        S: Iterator<Item = &'a u8> + DoubleEndedIterator + ExactSizeIterator,
    {
        if is_alt {
            self.alt_len += seq.len() as u64;
        }
        if front {
            for c in seq.rev() {
                self.seq.push_front(*c);
            }
        } else {
            self.seq.extend(seq);
        }
    }
}

pub(crate) struct BreakendEmissionParams<'a> {
    alt_allele: Arc<AltAllele>,
    ref_offset: usize,
    ref_end: usize,
    read_emission: Rc<ReadEmission<'a>>,
}

impl<'a> RefBaseEmission for BreakendEmissionParams<'a> {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        self.alt_allele[i]
    }

    default_ref_base_emission!();
}

impl<'a> EmissionParameters for BreakendEmissionParams<'a> {
    default_emission!();

    #[inline]
    fn len_x(&self) -> usize {
        self.alt_allele.len()
    }
}

/// Modeling of breakends.
#[derive(Getters, CopyGetters, Debug, Clone)]
pub(crate) struct Breakend {
    #[getset(get = "pub(crate)")]
    locus: genome::Locus,
    #[getset(get = "pub(crate)")]
    ref_allele: Vec<u8>,
    #[getset(get = "pub(crate)")]
    replacement: Vec<u8>,
    #[getset(get = "pub(crate)")]
    join: Option<Join>,
    #[getset(get_copy = "pub(crate)")]
    is_left_to_right: bool,
    #[getset(get = "pub(crate)")]
    id: Vec<u8>,
    #[getset(get = "pub(crate)")]
    mateid: Option<Vec<u8>>,
}

impl Breakend {
    pub(crate) fn new(
        locus: genome::Locus,
        ref_allele: &[u8],
        spec: &[u8],
        id: &[u8],
        mateid: Option<Vec<u8>>,
    ) -> Result<Option<Self>> {
        lazy_static! {
            static ref RE: Regex = Regex::new("((?P<replacement>[ACGTN]+)|((?P<bracket1>[\\]\\[])(?P<anglebracket1><)?(?P<contig>[^\\]\\[:>]+)(?P<anglebracket2>>)?(:(?P<pos>[0-9]+))?(?P<bracket2>[\\]\\[])))").unwrap();
        }
        lazy_static! {
            static ref SINGLE_RE: Regex =
                Regex::new("(\\.(?P<from_right>[ACGTN]+))|((?P<from_left>[ACGTN]+)\\.)").unwrap();
        }

        let spec = str::from_utf8(spec).unwrap().to_owned();

        let single_ops: Vec<_> = SINGLE_RE.captures_iter(&spec).collect();
        let ops: Vec<_> = RE.captures_iter(&spec).collect();
        if single_ops.len() == 1 {
            // Parse a single breakend (e.g. .ACTT)

            let caps = &single_ops[0];
            let (is_left_to_right, replacement) = if let Some(repl) = caps.name("from_left") {
                (true, repl)
            } else if let Some(repl) = caps.name("from_right") {
                (false, repl)
            } else {
                return Err(Error::InvalidBNDRecordAlt { spec }.into());
            };

            Ok(Some(Breakend {
                locus,
                ref_allele: ref_allele.to_owned(),
                replacement: replacement.as_str().as_bytes().to_owned(),
                join: None,
                is_left_to_right,
                id: id.to_owned(),
                mateid: None,
            }))
        } else {
            // parse a normal breakend

            if ops.len() != 2 {
                return Err(Error::InvalidBNDRecordAlt { spec }.into());
            }

            let mut replacement = None;
            let mut join = None;
            let mut is_left_to_right = false;
            for caps in ops {
                if let Some(repl) = caps.name("replacement") {
                    let repl = repl.as_str().as_bytes().to_owned();
                    if join.is_none() {
                        is_left_to_right = true;
                    }
                    replacement = Some(repl);
                } else {
                    // Must be in second case.
                    let bracket = caps.name("bracket1").unwrap().as_str();
                    if bracket != caps.name("bracket2").unwrap().as_str() {
                        return Err(Error::InvalidBNDRecordAlt { spec }.into());
                    }
                    // insertion from assembly file
                    match (
                        caps.name("anglebracket1").is_some(),
                        caps.name("anglebracket2").is_some(),
                    ) {
                        (true, true) => {
                            // insertion from assembly file
                            info!(
                                "Skipping BND at {}:{} pointing to assembly file",
                                locus.contig(),
                                locus.pos()
                            );
                            return Ok(None);
                        }
                        (false, false) => (), // no insertion from assembly file,
                        _ => {
                            // angle brackets do not match
                            return Err(Error::InvalidBNDRecordAlt { spec }.into());
                        }
                    }

                    let contig = caps.name("contig").unwrap().as_str().to_owned();
                    let pos = caps.name("pos").unwrap().as_str().parse::<u64>()? - 1; // pos is 1-based, need to convert to 0-based

                    let side = if bracket == "[" {
                        Side::RightOfPos
                    } else {
                        Side::LeftOfPos
                    };
                    let extension_modification = if is_left_to_right {
                        // Replacement sequence is inserted before.
                        // ] means reverse complement.
                        if bracket == "[" {
                            ExtensionModification::None
                        } else {
                            ExtensionModification::ReverseComplement
                        }
                    } else {
                        // Replacement sequence is inserted next.
                        // [ means reverse complement.
                        if bracket == "[" {
                            ExtensionModification::ReverseComplement
                        } else {
                            ExtensionModification::None
                        }
                    };

                    join = Some(Join {
                        locus: genome::Locus::new(contig, pos),
                        side,
                        extension_modification,
                    });
                }
            }

            Ok(Some(Breakend {
                locus,
                ref_allele: ref_allele.to_owned(),
                replacement: replacement.unwrap(),
                join,
                is_left_to_right,
                id: id.to_owned(),
                mateid,
            }))
        }
    }

    pub(crate) fn from_operations(
        locus: genome::Locus,
        ref_allele: &[u8],
        replacement: Vec<u8>,
        join: Join,
        is_left_to_right: bool,
        id: &[u8],
        mateid: &[u8],
    ) -> Self {
        Breakend {
            locus,
            ref_allele: ref_allele.to_owned(),
            replacement,
            join: Some(join),
            is_left_to_right,
            id: id.to_owned(),
            mateid: Some(mateid.to_owned()),
        }
    }

    pub(crate) fn spec(&self) -> Vec<u8> {
        let mut prefix = "";
        let mut suffix = "";

        if let Some(ref join) = self.join {
            let bracket;
            if self.is_left_to_right() {
                prefix = str::from_utf8(&self.replacement).unwrap();
                if self.emits_revcomp() {
                    bracket = ']'
                } else {
                    bracket = '['
                };
            } else {
                if self.emits_revcomp() {
                    bracket = '[';
                } else {
                    bracket = ']';
                };
                suffix = str::from_utf8(&self.replacement).unwrap();
            }

            format!(
                "{prefix}{bracket}{contig}:{pos}{bracket}{suffix}",
                prefix = prefix,
                suffix = suffix,
                bracket = bracket,
                contig = join.locus.contig(),
                pos = join.locus.pos() + 1
            )
            .as_bytes()
            .to_owned()
        } else {
            if self.is_left_to_right() {
                suffix = "."
            } else {
                prefix = "."
            }

            format!(
                "{prefix}{replacement}{suffix}",
                prefix = prefix,
                suffix = suffix,
                replacement = str::from_utf8(&self.replacement).unwrap()
            )
            .as_bytes()
            .to_owned()
        }
    }

    fn emits_revcomp(&self) -> bool {
        if let Some(Join {
            extension_modification: ExtensionModification::ReverseComplement,
            ..
        }) = self.join()
        {
            true
        } else {
            false
        }
    }

    pub(crate) fn to_variant(&self, event: &[u8]) -> model::Variant {
        model::Variant::Breakend {
            ref_allele: self.ref_allele.clone(),
            spec: self.spec(),
            event: event.to_owned(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Side {
    LeftOfPos,
    RightOfPos,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExtensionModification {
    None,
    ReverseComplement,
}

impl ExtensionModification {
    fn invert(&self) -> Self {
        match self {
            ExtensionModification::None => ExtensionModification::ReverseComplement,
            ExtensionModification::ReverseComplement => ExtensionModification::None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, new)]
pub(crate) struct Join {
    locus: genome::Locus,
    side: Side,
    extension_modification: ExtensionModification,
}

#[derive(Default, Debug)]
pub(crate) struct BreakendIndex {
    last_records: HashMap<Vec<u8>, usize>,
}

impl BreakendIndex {
    pub(crate) fn new<P: AsRef<Path>>(inbcf: P) -> Result<Self> {
        let mut bcf_reader = bcf::Reader::from_path(inbcf)?;
        if !utils::is_sv_bcf(&bcf_reader) {
            return Ok(BreakendIndex::default());
        }

        let mut last_records = HashMap::new();

        let mut i = 0;
        loop {
            let mut record = bcf_reader.empty_record();
            if !bcf_reader.read(&mut record)? {
                return Ok(BreakendIndex { last_records });
            }

            if utils::is_bnd(&mut record)? {
                // TODO support records without EVENT tag.
                if let Ok(Some(event)) = record.info(b"EVENT").string() {
                    let event = event[0];
                    last_records.insert(event.to_owned(), i);
                }
            }

            i += 1;
        }
    }

    pub(crate) fn last_record_index(&self, event: &[u8]) -> Option<usize> {
        self.last_records.get(event).cloned()
    }
}

struct LocusPlusOne<'a>(&'a genome::Locus);

impl<'a> AbstractLocus for LocusPlusOne<'a> {
    fn contig(&self) -> &str {
        self.0.contig()
    }

    fn pos(&self) -> u64 {
        self.0.pos() + 1
    }
}

struct LocusMinusOne<'a>(&'a genome::Locus);

impl<'a> AbstractLocus for LocusMinusOne<'a> {
    fn contig(&self) -> &str {
        self.0.contig()
    }

    fn pos(&self) -> u64 {
        self.0.pos() - 1
    }
}
