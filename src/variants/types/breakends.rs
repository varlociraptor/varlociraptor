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
use rust_htslib::bcf::{self, Read};

use crate::errors::Error;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::reference;
use crate::utils;
use crate::variants::evidence::realignment::pairhmm::{ReadEmission, RefBaseEmission};
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::model;
use crate::variants::sampling_bias::{ReadSamplingBias, SamplingBias};
use crate::variants::types::{
    AlleleSupport, MultiLocus, PairedEndEvidence, SingleLocusBuilder, Variant,
};
use crate::{default_emission, default_ref_base_emission};

pub(crate) struct BreakendGroup {
    loci: MultiLocus,
    breakends: BTreeMap<genome::Locus, Breakend>,
    alt_alleles: RefCell<HashMap<Vec<usize>, Vec<Rc<AltAllele>>>>,
    realigner: RefCell<Realigner>,
}

impl BreakendGroup {
    pub(crate) fn new(realigner: Realigner) -> Self {
        BreakendGroup {
            loci: MultiLocus::default(),
            breakends: BTreeMap::new(),
            alt_alleles: RefCell::default(),
            realigner: RefCell::new(realigner),
        }
    }

    pub(crate) fn breakends(&self) -> impl Iterator<Item = &Breakend> {
        self.breakends.values()
    }

    pub(crate) fn push(&mut self, breakend: Breakend) {
        let interval = genome::Interval::new(
            breakend.locus.contig().to_owned(),
            breakend.locus.pos()..breakend.locus.pos() + breakend.ref_allele.len() as u64,
        );

        self.breakends.insert(breakend.locus.clone(), breakend);

        self.loci.push(
            SingleLocusBuilder::default()
                .interval(interval)
                .build()
                .unwrap(),
        );
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

    fn contained_breakend_indices(&self, ref_interval: &genome::Interval) -> Vec<usize> {
        self.breakends
            .keys()
            .enumerate()
            .filter_map(|(i, locus)| {
                // TODO add genome::Interval::contains(genome::Locus) method to genome::Interval in bio-types.
                // Then, simplify this here (see PR https://github.com/rust-bio/rust-bio-types/pull/9).
                if ref_interval.contig() == locus.contig()
                    && locus.pos() >= ref_interval.range().start
                    && locus.pos() < ref_interval.range().end
                {
                    Some(i)
                } else {
                    None
                }
            })
            .collect()
    }
}

impl Variant for BreakendGroup {
    type Evidence = PairedEndEvidence;
    type Loci = MultiLocus;

    fn is_valid_evidence(&self, evidence: &Self::Evidence) -> Option<Vec<usize>> {
        let overlapping: Vec<_> = match evidence {
            PairedEndEvidence::SingleEnd(read) => self
                .loci
                .iter()
                .enumerate()
                .filter_map(|(i, locus)| {
                    if !locus.overlap(read, true).is_none() {
                        Some(i)
                    } else {
                        None
                    }
                })
                .collect(),
            PairedEndEvidence::PairedEnd { left, right } => self
                .loci
                .iter()
                .enumerate()
                .filter_map(|(i, locus)| {
                    if !locus.overlap(left, true).is_none() || !locus.overlap(right, true).is_none()
                    {
                        Some(i)
                    } else {
                        None
                    }
                })
                .collect(),
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
    fn enclosable_len(&self) -> Option<u64> {
        let first = self.breakends.keys().next().unwrap();
        if self
            .breakends
            .values()
            .skip(1)
            .all(|bnd| bnd.locus.contig() == first.contig())
        {
            let last = self.breakends.values().last().unwrap();
            Some(last.locus.pos() + last.ref_allele.len() as u64 - first.pos())
        } else {
            None
        }
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
        ref_interval: &genome::Interval,
        ref_window: usize,
    ) -> Result<Vec<BreakendEmissionParams<'a>>> {
        // Step 1: fetch contained breakends
        let bnds = self.contained_breakend_indices(ref_interval);

        if !self.alt_alleles.borrow().contains_key(&bnds) {
            // Compute alt allele sequence once.

            let mut candidate_alt_alleles = Vec::new();

            for (i, first) in self.breakends.values().enumerate() {
                if !bnds.contains(&i) {
                    continue;
                }

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
                    );
                    //alt_allele.push_seq(b"1".iter(), false); // dbg
                    // Add replacement to alt allele.
                    alt_allele.push_seq(first.replacement().iter(), false);
                } else {
                    // Add suffix on reference.
                    let ref_seq = ref_buffer.seq(first.locus.contig())?;
                    alt_allele.push_seq(
                        ref_seq[suffix_range(first, ref_seq.len(), false)].iter(),
                        true,
                    );
                    //alt_allele.push_seq(b"1".iter(), true); // dbg
                    // Add replacement to alt allele.
                    alt_allele.push_seq(first.replacement().iter(), true);
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
                            );
                        } else {
                            // Prepend prefix on reference.
                            alt_allele.push_seq(ref_seq[prefix_range(current)].iter(), true);
                        }
                        break;
                    }
                    visited.insert(&current.id);

                    if let Operation::Join {
                        ref locus,
                        side,
                        extension_modification,
                    } = current.join() {
                        let ref_seq = ref_buffer.seq(locus.contig())?;
                        // Find next breakend from here in order to know how long to extend.
                        let seq = match side {
                            Side::LeftOfPos => {
                                next_bnd = self.upstream_bnd(locus);
                                let seq_start =
                                    next_bnd.map_or(0, |bnd| bnd.locus.pos() as usize + 1); // end before the beginning of the next bnd
                                let seq_end = locus.pos() as usize + 1; // locus.pos() is meant inclusive, so we need to add 1 here.
                                &ref_seq[seq_start..seq_end]
                            }
                            Side::RightOfPos => {
                                next_bnd = self.downstream_bnd(locus);
                                let seq_start = locus.pos() as usize; // locus.pos() is meant inclusive.
                                let seq_end = next_bnd
                                    .map_or(ref_seq.len(), |bnd| bnd.locus.pos() as usize); // end before the start of the next bnd
                                &ref_seq[seq_start..seq_end]
                            }
                        };

                        // Apply replacement operation of next bnd if necessary.
                        let seq: Vec<u8> = match (next_bnd, side) {
                            (Some(next_bnd), Side::RightOfPos) => seq.iter().chain(next_bnd.replacement().iter()).cloned().collect(),
                            (Some(next_bnd), Side::LeftOfPos) => next_bnd.replacement().iter().chain(seq.iter()).cloned().collect(),
                            (None, _) => seq.to_owned(),
                        };

                        //dbg!((&current, next_bnd, std::str::from_utf8(&seq).unwrap()));

                        let (left_to_right, extension_modification) = if revcomp {
                            (!current.is_left_to_right(), extension_modification.invert())
                        } else {
                            (current.is_left_to_right(), *extension_modification)
                        };

                        //alt_allele.push_seq(b"|".iter(), !left_to_right); // dbg

                        // Push sequence to alt allele.
                        match (
                            extension_modification,
                            next_bnd.is_none(),
                            left_to_right,
                        ) {
                            (ExtensionModification::None, false, is_left_to_right) => {
                                alt_allele.push_seq(seq.iter(), !is_left_to_right)
                            }
                            (
                                ExtensionModification::ReverseComplement,
                                false,
                                is_left_to_right,
                            ) => {
                                alt_allele
                                    .push_seq(dna::revcomp(&seq).iter(), !is_left_to_right);
                            }
                            (ExtensionModification::None, true, true) => {
                                alt_allele.push_seq(seq[..ref_window].iter(), false);
                            }
                            (ExtensionModification::ReverseComplement, true, true) => {
                                alt_allele
                                    .push_seq(dna::revcomp(seq)[..ref_window].iter(), false);
                            }
                            (ExtensionModification::None, true, false) => {
                                alt_allele.push_seq(
                                    seq[seq.len().saturating_sub(ref_window)..].iter(),
                                    true,
                                );
                            }
                            (ExtensionModification::ReverseComplement, true, false) => {
                                let seq = dna::revcomp(seq);
                                alt_allele.push_seq(
                                    seq[seq.len().saturating_sub(ref_window)..].iter(),
                                    true,
                                );
                            }
                        }

                        // Update revcomp marker for next iteration.
                        revcomp = if let ExtensionModification::ReverseComplement = extension_modification {
                            true
                        } else {
                            false
                        };
                    } else {
                        unreachable!();
                    }
                }

                candidate_alt_alleles.push(Rc::new(alt_allele));
            }
            self.alt_alleles.borrow_mut().insert(bnds.clone(), candidate_alt_alleles);
        }

        Ok(self.alt_alleles.borrow().get(&bnds).unwrap().iter().map(|alt_allele| {
            let alt_allele = Rc::clone(alt_allele);

            dbg!(std::str::from_utf8(&Vec::from((**alt_allele.as_ref()).clone())).unwrap());

            BreakendEmissionParams {
                ref_offset: 0,
                ref_end: alt_allele.len(),
                alt_allele,
                read_emission: Rc::clone(&read_emission_params),
            }
        }).collect())
    }
}

#[derive(Derefable, Debug, Default)]
pub(crate) struct AltAllele {
    #[deref]
    seq: VecDeque<u8>,
}

impl AltAllele {
    pub(crate) fn push_seq<'a, S>(&mut self, seq: S, front: bool)
    where
        S: Iterator<Item = &'a u8> + DoubleEndedIterator,
    {
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
    alt_allele: Rc<AltAllele>,
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
#[derive(Getters, Debug)]
pub(crate) struct Breakend {
    #[getset(get = "pub")]
    locus: genome::Locus,
    ref_allele: Vec<u8>,
    operations: [Operation; 2],
    spec: Option<Vec<u8>>,
    #[getset(get = "pub")]
    id: Vec<u8>,
    #[getset(get = "pub")]
    mateid: Vec<u8>,
}

impl Breakend {
    pub(crate) fn new(
        locus: genome::Locus,
        ref_allele: &[u8],
        spec: &[u8],
        id: &[u8],
        mateid: &[u8],
    ) -> Result<Option<Self>> {
        lazy_static! {
            static ref RE: Regex = Regex::new("((?P<replacement>[ACGTN]+)|((?P<bracket1>[\\]\\[])(?P<anglebracket1><)?(?P<contig>[^\\]\\[:>]+)(?P<anglebracket2>>)?(:(?P<pos>[0-9]+))?(?P<bracket2>[\\]\\[])))").unwrap();
        }

        let spec = str::from_utf8(spec).unwrap().to_owned();

        let ops: Vec<_> = RE.captures_iter(&spec).collect();
        if ops.len() != 2 {
            return Err(Error::InvalidBNDRecordAlt { spec }.into());
        }

        let mut operations = Vec::new();
        for caps in ops {
            if let Some(replacement) = caps.name("replacement") {
                let replacement = replacement.as_str().as_bytes().to_owned();
                operations.push(Operation::Replacement(replacement));
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
                let pos: u64 = caps.name("pos").unwrap().as_str().parse()?;

                let side = if bracket == "[" {
                    Side::RightOfPos
                } else {
                    Side::LeftOfPos
                };
                let extension_modification = if !operations.is_empty() {
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

                operations.push(Operation::Join {
                    locus: genome::Locus::new(contig, pos),
                    side,
                    extension_modification,
                });
            }
        }

        Ok(Some(Breakend {
            locus,
            ref_allele: ref_allele.to_owned(),
            operations: [operations[0].clone(), operations[1].clone()],
            spec: Some(spec.as_bytes().to_owned()),
            id: id.to_owned(),
            mateid: mateid.to_owned(),
        }))
    }

    pub(crate) fn from_operations(
        locus: genome::Locus,
        ref_allele: &[u8],
        operations: [Operation; 2],
        id: &[u8],
        mateid: &[u8],
    ) -> Self {
        Breakend {
            locus,
            ref_allele: ref_allele.to_owned(),
            operations,
            spec: None,
            id: id.to_owned(),
            mateid: mateid.to_owned(),
        }
    }

    fn is_left_to_right(&self) -> bool {
        if let Operation::Replacement(_) = self.operations[0] {
            true
        } else {
            false
        }
    }

    fn replacement(&self) -> &[u8] {
        for op in &self.operations {
            if let Operation::Replacement(seq) = op {
                return &seq
            }
        }
        unreachable!();
    }

    fn join(&self) -> &Operation {
        for op in &self.operations {
            if let Operation::Join{ .. } = op {
                return op
            }
        }
        unreachable!();
    }

    fn emits_revcomp(&self) -> bool {
        for op in &self.operations {
            if let Operation::Join {
                extension_modification: ExtensionModification::ReverseComplement,
                ..
            } = op
            {
                return true;
            }
        }

        false
    }

    pub(crate) fn to_variant(&self, event: &[u8]) -> Option<model::Variant> {
        if let Some(ref spec) = self.spec {
            Some(model::Variant::Breakend {
                ref_allele: self.ref_allele.clone(),
                spec: spec.to_owned(),
                event: event.to_owned(),
            })
        } else {
            None
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Operation {
    Replacement(Vec<u8>),
    Join {
        locus: genome::Locus,
        side: Side,
        extension_modification: ExtensionModification,
    },
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