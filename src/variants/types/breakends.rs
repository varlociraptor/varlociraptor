// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cell::RefCell;
use std::str;

use anyhow::Result;
use regex::Regex;
use bio_types::genome::{self, AbstractLocus};

use crate::errors::Error;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::evidence::realignment::pairhmm::{ReadEmission, RefBaseEmission};
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::sampling_bias::{ReadSamplingBias, SamplingBias};
use crate::variants::types::{AlleleSupport, MultiLocus, PairedEndEvidence, SingleLocus, Variant};
use crate::{default_emission, default_ref_base_emission};

#[derive(Debug, new)]
pub(crate) struct Breakends {
    #[new(default)]
    loci: MultiLocus,
    #[new(default)]
    breakends: Vec<Breakend>,
    realigner: RefCell<Realigner>,
}

impl Breakends {
    pub(crate) fn push(&mut self, breakend: Breakend) {
        self.loci.push(
            SingleLocus::new(
                genome::Interval::new(
                    breakend.locus.contig().to_owned(),
                    breakend.locus.pos()..breakend.locus.pos() + breakend.ref_allele.len() as u64
                )
            )
        );

        self.breakends.push(breakend);
    }
}


/// Modeling of breakends.
#[derive(Debug)]
pub(crate) struct Breakend {
    locus: genome::Locus,
    ref_allele: Vec<u8>,
    operations: [Operation; 2],
}

impl Breakend {
    pub(crate) fn new(locus: genome::Locus, ref_allele: &[u8], spec: &[u8]) -> Result<Option<Self>> {
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
        }))
    }
}


#[derive(Debug, Clone)]
pub(crate) enum Side {
    LeftOfPos,
    RightOfPos,
}

#[derive(Debug, Clone)]
pub(crate) enum ExtensionModification {
    None,
    ReverseComplement,
}

#[derive(Debug, Clone)]
pub(crate) enum Operation {
    Replacement(Vec<u8>),
    Join { 
        locus: genome::Locus,
        side: Side,
        extension_modification: ExtensionModification,
    },
}

