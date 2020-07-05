// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cell::RefCell;
use std::collections::{BTreeMap, HashMap};
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
use crate::variants::types::{AlleleSupport, MultiLocus, PairedEndEvidence, SingleLocus, Variant};
use crate::{default_emission, default_ref_base_emission};

pub(crate) struct BreakendGroup {
    loci: MultiLocus,
    breakends: BTreeMap<genome::Locus, Breakend>,
    alt_alleles: RefCell<HashMap<Vec<usize>, Rc<Vec<u8>>>>,
    realigner: RefCell<Realigner>,
}

impl BreakendGroup {
    pub(crate) fn new(realigner: Realigner) -> Self {
        BreakendGroup {
            loci: MultiLocus::default(),
            breakends: BTreeMap::default(),
            alt_alleles: RefCell::default(),
            realigner: RefCell::new(realigner),
        }
    }

    pub(crate) fn breakends(&self) -> impl Iterator<Item = &Breakend> {
        self.breakends.values()
    }

    pub(crate) fn push(
        &mut self,
        locus: genome::Locus,
        ref_allele: &[u8],
        spec: &[u8],
        id: Vec<u8>,
        mateid: Vec<u8>,
    ) -> Result<()> {
        if let Some(breakend) = Breakend::new(locus, ref_allele, spec, id, mateid)? {
            let interval = genome::Interval::new(
                breakend.locus.contig().to_owned(),
                breakend.locus.pos()..breakend.locus.pos() + ref_allele.len() as u64,
            );

            self.breakends.insert(breakend.locus.clone(), breakend);

            self.loci.push(SingleLocus::new(interval));
        }

        Ok(())
    }

    fn upstream_bnd(&self, locus: &genome::Locus) -> Option<&Breakend> {
        self.breakends.range(..locus).last().and_then(|(_, bnd)| {
            if bnd.locus.contig() == locus.contig() {
                Some(bnd)
            } else {
                None
            }
        })
    }

    fn downstream_bnd(&self, locus: &genome::Locus) -> Option<&Breakend> {
        self.breakends.range(locus..).nth(1).and_then(|(_, bnd)| {
            if bnd.locus.contig() == locus.contig() {
                Some(bnd)
            } else {
                None
            }
        })
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

    fn alt_emission_params(
        &self,
        read_emission_params: Rc<ReadEmission<'a>>,
        ref_buffer: Arc<reference::Buffer>,
        ref_interval: &genome::Interval,
        ref_window: usize,
    ) -> Result<BreakendEmissionParams<'a>> {
        // Step 1: fetch contained breakends
        let bnds: Vec<_> = self
            .breakends
            .keys()
            .enumerate()
            .filter_map(|(i, locus)| {
                if ref_interval.contig() == locus.contig()
                    && locus.pos() >= ref_interval.range().start
                    && locus.pos() < ref_interval.range().end
                {
                    Some(i)
                } else {
                    None
                }
            })
            .collect();

        if !self.alt_alleles.borrow().contains_key(&bnds) {
            // Compute alt allele sequence once.
            let first = self.breakends.values().nth(*bnds.first().unwrap()).unwrap();

            let mut alt_allele = Vec::new();

            // prefix on reference
            let start = first.locus.pos() as usize;
            alt_allele.extend(
                &ref_buffer.seq(first.locus.contig())?[start.saturating_sub(ref_window)..start],
            );

            // breakend operations in between
            let mut next_bnd = Some(first);
            while let Some(current) = next_bnd {
                // apply operation
                for op in &current.operations {
                    match op {
                        Operation::Replacement(seq) => alt_allele.extend(seq),
                        Operation::Join {
                            ref locus,
                            side,
                            extension_modification,
                        } => {
                            let ref_seq = ref_buffer.seq(locus.contig())?;
                            // Find next breakend from here in order to know how long to extend and which to evaluate next.
                            let seq = match side {
                                Side::LeftOfPos => {
                                    next_bnd = self.upstream_bnd(locus);
                                    let seq_start =
                                        next_bnd.map_or(0, |bnd| bnd.locus.pos() as usize + 1);
                                    let seq_end = locus.pos() as usize;

                                    &ref_seq[seq_start..seq_end]
                                }
                                Side::RightOfPos => {
                                    next_bnd = self.downstream_bnd(locus);
                                    let seq_start = locus.pos() as usize;
                                    let seq_end = next_bnd
                                        .map_or(ref_seq.len(), |bnd| bnd.locus.pos() as usize);

                                    &ref_seq[seq_start..seq_end]
                                }
                            };

                            match (extension_modification, next_bnd.is_none()) {
                                (ExtensionModification::None, false) => alt_allele.extend(seq),
                                (ExtensionModification::None, true) => {
                                    alt_allele.extend(&seq[..ref_window])
                                }
                                (ExtensionModification::ReverseComplement, false) => {
                                    alt_allele.extend(dna::revcomp(seq))
                                }
                                (ExtensionModification::ReverseComplement, true) => {
                                    alt_allele.extend(&dna::revcomp(seq)[..ref_window])
                                }
                            }
                        }
                    }
                }
            }
            self.alt_alleles
                .borrow_mut()
                .insert(bnds.clone(), Rc::new(alt_allele));
        }

        let alt_allele = Rc::clone(self.alt_alleles.borrow().get(&bnds).unwrap());

        Ok(BreakendEmissionParams {
            ref_offset: 0,
            ref_end: alt_allele.len(),
            alt_allele,
            read_emission: read_emission_params,
        })
    }
}

#[derive(Debug)]
pub(crate) struct BreakendEmissionParams<'a> {
    alt_allele: Rc<Vec<u8>>,
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
#[derive(Getters)]
pub(crate) struct Breakend {
    #[getset(get = "pub")]
    locus: genome::Locus,
    ref_allele: Vec<u8>,
    operations: [Operation; 2],
    spec: Vec<u8>,
    #[getset(get = "pub")]
    id: Vec<u8>,
    #[getset(get = "pub")]
    mateid: Vec<u8>,
}

impl Breakend {
    fn new(
        locus: genome::Locus,
        ref_allele: &[u8],
        spec: &[u8],
        id: Vec<u8>,
        mateid: Vec<u8>,
    ) -> Result<Option<Self>> {
        lazy_static! {
            static ref RE: Regex = Regex::new("((?P<replacement>[ACGTN]+)|((?P<bracket1>[\\]\\[])(?P<anglebracket1><)?(?P<contig>[^\\]\\[])(?P<anglebracket2>>)?(:(?P<pos>[^\\]\\[]))?(?P<bracket2>[\\]\\[])))").unwrap();
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
            spec: spec.as_bytes().to_owned(),
            id,
            mateid,
        }))
    }

    pub(crate) fn to_variant(&self, event: &[u8]) -> model::Variant {
        model::Variant::Breakend {
            ref_allele: self.ref_allele.clone(),
            spec: self.spec.clone(),
            event: event.to_owned(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Side {
    LeftOfPos,
    RightOfPos,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum ExtensionModification {
    None,
    ReverseComplement,
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
                if let Some(event) = record.info(b"EVENT").string()? {
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
