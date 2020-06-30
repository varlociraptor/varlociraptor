// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cell::RefCell;
use std::str;
use std::sync::Arc;
use std::rc::Rc;
use std::cmp;
use std::collections::BTreeMap;

use anyhow::Result;
use regex::Regex;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractLocus};
use bio::stats::pairhmm::EmissionParameters;

use crate::errors::Error;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::evidence::realignment::pairhmm::{ReadEmission, RefBaseEmission};
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::sampling_bias::{ReadSamplingBias, SamplingBias};
use crate::variants::types::{AlleleSupport, MultiLocus, PairedEndEvidence, SingleLocus, Variant};
use crate::{default_emission, default_ref_base_emission};
use crate::reference;

#[derive(new)]
pub(crate) struct BreakendGroup<'a> {
    #[new(default)]
    loci: MultiLocus,
    #[new(default)]
    breakends: BTreeMap<genome::Locus, Breakend<'a>>,
    realigner: RefCell<Realigner>,
}

impl<'a> BreakendGroup<'a> {
    pub(crate) fn push(&'a mut self, locus: genome::Locus, ref_allele: &[u8], spec: &[u8]) -> Result<()> {
        self.loci.push(
            SingleLocus::new(
                genome::Interval::new(
                    locus.contig().to_owned(),
                    locus.pos()..locus.pos() + ref_allele.len() as u64
                )
            )
        );

        if let Some(breakend) = Breakend::new(locus, ref_allele, spec, self)? {
            self.breakends.push(breakend);
        }

        Ok(())
    }
}

impl<'a> Realignable<'a> for BreakendGroup {
    type EmissionParams = BreakendGroupEmissionParams<'a>;

    fn alt_emission_params(
        &self,
        read_emission_params: Rc<ReadEmission<'a>>,
        ref_buffer: Arc<reference::Buffer>,
        ref_window: usize,
    ) -> Result<BreakendGroupEmissionParams<'a>> {
        if !self.alt_alleles.contains_key(origin_contig) {
            // Step 1: find breakend that starts close to


            // Compute alt allele sequence once.
            let first = self.breakends.first().unwrap();
            
            let mut alt_allele = Rc::new(Vec::new());
            
            // prefix on reference
            let start = first.locus.pos() as usize;
            alt_allele.extend(&ref_buffer.seq(first.locus.contig())?[start.saturating_sub(ref_window)..start]);

            // breakend operations in between
            for breakend in self.breakends {
                breakend.
            }

            // suffix on reference
            let last = self.breakends.last().unwrap();
            let last_ref_seq = ref_buffer.seq(last.locus.contig())?;
            let end = last.locus.pos() as usize + last.ref_allele.len();
            alt_allele.extend(&last_ref_seq[end..cmp::min(end + ref_window, last_ref_seq.len()]);

            self.alt_allele = Some(alt_allele);
        }

        Ok(BreakendGroupEmissionParams {
            alt_allele: Rc::clone(self.alt_alleles.get(origin_contig).unwrap()),
        })
    }
}

#[derive(Debug)]
pub(crate) struct BreakendEmissionParams<'a> {
    alt_allele: Rc<Vec<u8>>,
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
pub(crate) struct Breakend<'a> {
    locus: genome::Locus,
    ref_allele: Vec<u8>,
    operations: [Operation; 2],
    // lazy bucket for alt allele sequence (used for realignment)
    alt_allele: Option<Rc<Vec<u8>>>,
    group: &'a BreakendGroup<'a>,
}

impl<'a> Breakend<'a> {
    fn new(locus: genome::Locus, ref_allele: &[u8], spec: &[u8], group: &'a BreakendGroup) -> Result<Option<Self>> {
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
                match (caps.name("anglebracket1").is_some(), caps.name("anglebracket2").is_some()) {
                    (true, true) => {
                        // insertion from assembly file
                        info!("Skipping BND at {}:{} pointing to assembly file", locus.contig(), locus.pos());
                        return Ok(None);
                    },
                    (false, false) => (), // no insertion from assembly file,
                    _ => {
                        // angle brackets do not match
                        return Err(Error::InvalidBNDRecord { spec }.into());
                    }
                }

                let contig = caps.name("contig").unwrap().as_str().to_owned();
                let pos: u64 = caps.name("pos").unwrap().as_str().parse()?;

                let side = if bracket == "[" { Side::RightOfPos } else { Side::LeftOfPos };
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
            alt_allele: None,
            group,
        }))
    }
}

impl<'a, 'b> Realignable<'a> for Breakend<'b> {
    type EmissionParams = BreakendEmissionParams<'a>;

    fn alt_emission_params(
        &self,
        read_emission_params: Rc<ReadEmission<'a>>,
        ref_buffer: Arc<reference::Buffer>,
        ref_window: usize,
    ) -> Result<BreakendEmissionParams<'a>> {

        if self.alt_allele.is_none() {
            // Step 1: backtrack from here to find first breakend on reference.
            let backtrack = |breakend: Breakend| {
                let mate_locus = breakend.mate_locus();
                let mate self.group.breakends.get(mate_locus).unwrap();
                if mate.locus > breakend.locus {
                    // TODO stop, we are already at the leftmost breakend.
                } else {
                    backtrack(&mate);
                }
            }

            // Step 2: obtain prefix on reference.

            // Step 3: apply breakend operations.

            // Step 4: obtain suffix on reference.

            // Compute alt allele sequence once.
            let first = self.group.breakends.first().unwrap();
            
            let mut alt_allele = Rc::new(Vec::new());
            
            // prefix on reference
            let start = first.locus.pos() as usize;
            alt_allele.extend(&ref_buffer.seq(first.locus.contig())?[start.saturating_sub(ref_window)..start]);

            // breakend operations in between
            for breakend in self.breakends {
                breakend.
            }

            // suffix on reference
            let last = self.breakends.last().unwrap();
            let last_ref_seq = ref_buffer.seq(last.locus.contig())?;
            let end = last.locus.pos() as usize + last.ref_allele.len();
            alt_allele.extend(&last_ref_seq[end..cmp::min(end + ref_window, last_ref_seq.len()]);

            self.alt_allele = Some(alt_allele);
        }

        Ok(BreakendEmissionParams {
            alt_allele: Rc::clone(self.alt_alleles.get(origin_contig).unwrap()),
        })
    }

    fn mate_locus(&self) -> &genome::Locus {
        for op in self.operations {
            if let Join { ref locus, ..} = op {
                return locus;
            }
        }
        unreachable!;
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

