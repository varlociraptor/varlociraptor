// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cell::RefCell;
use std::cmp;
use std::collections::{BTreeMap, HashMap};
use std::rc::Rc;
use std::str;
use std::sync::Arc;

use anyhow::Result;
use bio::alphabets::dna;
use bio::stats::pairhmm::EmissionParameters;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval, AbstractLocus};
use regex::Regex;

use crate::errors::Error;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::reference;
use crate::variants::evidence::realignment::pairhmm::{ReadEmission, RefBaseEmission};
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::sampling_bias::{ReadSamplingBias, SamplingBias};
use crate::variants::types::{AlleleSupport, MultiLocus, PairedEndEvidence, SingleLocus, Variant};
use crate::{default_emission, default_ref_base_emission};

#[derive(new)]
pub(crate) struct BreakendGroup {
    #[new(default)]
    loci: MultiLocus,
    #[new(default)]
    breakends: BTreeMap<genome::Locus, Breakend>,
    #[new(default)]
    alt_alleles: RefCell<HashMap<Vec<usize>, Rc<Vec<u8>>>>,
    realigner: RefCell<Realigner>,
}

impl BreakendGroup {
    pub(crate) fn push(
        &mut self,
        locus: genome::Locus,
        ref_allele: &[u8],
        spec: &[u8],
    ) -> Result<()> {
        if let Some(breakend) = Breakend::new(locus, ref_allele, spec)? {
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
        self.breakends
            .range(locus..)
            .skip(1)
            .next()
            .and_then(|(_, bnd)| {
                if bnd.locus.contig() == locus.contig() {
                    Some(bnd)
                } else {
                    None
                }
            })
    }
}

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
            let first = self
                .breakends
                .values()
                .skip(*bnds.first().unwrap())
                .next()
                .unwrap();
            let last = self
                .breakends
                .values()
                .skip(*bnds.last().unwrap())
                .next()
                .unwrap();

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
pub(crate) struct Breakend {
    locus: genome::Locus,
    ref_allele: Vec<u8>,
    operations: [Operation; 2],
}

impl Breakend {
    fn new(locus: genome::Locus, ref_allele: &[u8], spec: &[u8]) -> Result<Option<Self>> {
        lazy_static! {
            static ref RE: Regex = Regex::new("((?P<replacement>[ACGTN]+)|((?P<bracket1>[\\]\\[])(?P<anglebracket1><)?(?P<contig>[^\\]\\[])(?P<anglebracket2>>)?(:(?P<pos>[^\\]\\[]))?(?P<bracket2>[\\]\\[])))").unwrap();
        }

        let spec = str::from_utf8(spec).unwrap().to_owned();

        let ops: Vec<_> = RE.captures_iter(&spec).collect();
        if ops.len() != 2 {
            return Err(Error::InvalidBNDRecord { spec }.into());
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
                    return Err(Error::InvalidBNDRecord { spec }.into());
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
                        return Err(Error::InvalidBNDRecord { spec }.into());
                    }
                }

                let contig = caps.name("contig").unwrap().as_str().to_owned();
                let pos: u64 = caps.name("pos").unwrap().as_str().parse()?;

                let side = if bracket == "[" {
                    Side::RightOfPos
                } else {
                    Side::LeftOfPos
                };
                let extension_modification = if operations.len() > 0 {
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
        }))
    }

    fn mate_locus(&self) -> &genome::Locus {
        for op in &self.operations {
            if let Operation::Join { ref locus, .. } = op {
                return locus;
            }
        }
        unreachable!();
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
