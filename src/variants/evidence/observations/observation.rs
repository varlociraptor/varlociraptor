// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::char;
use std::hash::{Hash, Hasher};
use std::ops;
use std::rc::Rc;
use std::str;

use anyhow::Result;
use bio_types::genome::{self, AbstractLocus};
use bio_types::sequence::SequenceReadPairOrientation;
use counter::Counter;

use rust_htslib::bam;

use serde::Serialize;
// use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::bayesian::BayesFactor;
use itertools::Itertools;

use crate::errors::{self, Error};
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::utils::{self};
use crate::utils::bayes_factor_to_letter;



const INVALID_XA_FORMAT_MSG: &str = "XA tag of bam records in unexpected format. Expecting string (type Z) in bwa format (chr,pos,CIGAR,NM;).";

/// Strand support for observation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Strand {
    Forward,
    Reverse,
    Both,
    #[default]
    None,
}

impl Strand {
    pub(crate) fn from_record_and_pos(record: &bam::Record, pos: usize) -> Result<Self> {
        Ok(
            if let Some(strand_info) = utils::aux_tag_strand_info(record) {
                if let Some(s) = strand_info.get(pos) {
                    Self::from_aux_item(*s)?
                } else {
                    return Err(Error::ReadPosOutOfBounds.into());
                }
            } else {
                Self::from_record(record)
            },
        )
    }

    pub(crate) fn from_record(record: &bam::Record) -> Self {
        // just record global strand information of record
        let reverse_strand = utils::is_reverse_strand(record);
        if reverse_strand {
            Strand::Reverse
        } else {
            Strand::Forward
        }
    }

    pub(crate) fn from_aux(strand_info_aux: &[u8]) -> Result<Self> {
        let mut strand = Strand::None;
        for s in strand_info_aux {
            strand |= Self::from_aux_item(*s)?;
        }
        Ok(strand)
    }

    pub(crate) fn from_aux_item(item: u8) -> Result<Self> {
        Ok(match item {
            b'+' => Strand::Forward,
            b'-' => Strand::Reverse,
            b'*' => Strand::Both,
            b'.' => Strand::None,
            _ => {
                return Err(Error::InvalidStrandInfo {
                    value: char::from_u32(item as u32).unwrap(),
                }
                .into())
            }
        })
    }

    pub(crate) fn no_strand_info() -> Self {
        Self::None
    }
}

impl ops::BitOrAssign for Strand {
    fn bitor_assign(&mut self, rhs: Self) {
        if let Strand::None = self {
            *self = rhs;
        } else if let Strand::None = rhs {
            // no new information
        } else if *self != rhs {
            *self = Strand::Both;
        }
    }
}

pub(crate) fn read_orientation(record: &bam::Record) -> Result<SequenceReadPairOrientation> {
    if let Ok(bam::record::Aux::String(ro)) = record.aux(b"RO") {
        let orientations = ro.as_bytes().split(|e| *e == b',').collect_vec();
        Ok(if orientations.len() != 1 {
            // more than one orientation, return None
            SequenceReadPairOrientation::None
        } else {
            match orientations[0] {
                b"F1R2" => SequenceReadPairOrientation::F1R2,
                b"F2R1" => SequenceReadPairOrientation::F2R1,
                b"F1F2" => SequenceReadPairOrientation::F1F2,
                b"F2F1" => SequenceReadPairOrientation::F2F1,
                b"R1R2" => SequenceReadPairOrientation::R1R2,
                b"R2R1" => SequenceReadPairOrientation::R2R1,
                b"R1F2" => SequenceReadPairOrientation::R1F2,
                b"R2F1" => SequenceReadPairOrientation::R2F1,
                b"None" => SequenceReadPairOrientation::None,
                _ => {
                    return Err(errors::Error::InvalidReadOrientationInfo {
                        value: str::from_utf8(orientations[0]).unwrap().to_owned(),
                    }
                    .into())
                }
            }
        })
    } else {
        Ok(record.read_pair_orientation())
    }
}

#[derive(Debug, Clone, Derefable, Default)]
pub struct ExactAltLoci {
    #[deref]
    inner: Vec<genome::Locus>,
}

impl<'a> From<&'a bam::Record> for ExactAltLoci {
    fn from(record: &'a bam::Record) -> Self {
        match record.aux(b"XA") {
            Ok(bam::record::Aux::String(xa)) => {
                ExactAltLoci {
                    inner: xa
                        .split(';')
                        .filter_map(|xa| {
                            if xa.is_empty() {
                                // last semicolon passed
                                None
                            } else {
                                let items: Vec<_> = xa.split(',').collect();
                                if items.len() == 4 {
                                    let contig = items[0];
                                    let mut pos = items[1];
                                    if pos.starts_with('-') || pos.starts_with('-') {
                                        pos = &pos[1..];
                                    }
                                    if let Ok(pos) = pos.parse() {
                                        Some(genome::Locus::new(contig.to_owned(), pos))
                                    } else {
                                        None
                                    }
                                } else {
                                    warn!("{}", INVALID_XA_FORMAT_MSG);
                                    None
                                }
                            }
                        })
                        .collect(),
                }
            }
            Ok(_tag) => {
                warn!("{}", INVALID_XA_FORMAT_MSG);
                ExactAltLoci::default()
            }
            Err(_e) => {
                // no XA tag found, return empty.
                ExactAltLoci::default()
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) enum AltLocus {
    Major,
    Some,
    None,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum ReadPosition {
    Major,
    #[default]
    Some,
}

pub(crate) enum MaxBayesFactor {
    Alt(BayesFactor),
    Ref(BayesFactor),
    Equal,
}

impl std::fmt::Display for MaxBayesFactor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                MaxBayesFactor::Alt(bf) => format!("A{}", bayes_factor_to_letter(*bf)),
                MaxBayesFactor::Ref(bf) => format!("R{}", bayes_factor_to_letter(*bf)),
                MaxBayesFactor::Equal => "E".to_string(),
            }
        )
    }
}

pub(crate) fn calc_major_feature<F, I>(feature: I) -> Option<F>
where
    I: Iterator<Item = F>,
    F: Clone + Eq + Hash,
{
    let counter: Counter<_> = feature.collect();
    let most_common = counter.most_common();
    if most_common.is_empty() {
        None
    } else {
        let (ref most_common_feature, most_common_count) = most_common[0];
        if most_common_count == 1 {
            // all reads exhibit different features
            None
        } else if most_common.len() == 1 {
            // all the same
            Some(most_common_feature.clone())
        } else {
            let (ref _second_most_common_feature, second_most_common_count) = most_common[1];
            if most_common_count > second_most_common_count {
                // clear winner
                Some(most_common_feature.clone())
            } else {
                // tie
                None
            }
        }
    }
}

pub(crate) fn locus_to_bucket(
    locus: &genome::Locus,
    alignment_properties: &AlignmentProperties,
) -> genome::Locus {
    // METHOD: map each locus to the nearest multiple of the read len from the left.
    // This way, varying reads become comparable

    let coeff = alignment_properties.max_read_len as u64 * 10;

    let bucket = genome::Locus::new(locus.contig().to_owned(), (locus.pos() / coeff) * coeff);

    bucket
}

#[derive(Clone, Eq, Debug)]
pub(crate) enum Evidence {
    SingleEndSequencingRead(Rc<bam::Record>),
    PairedEndSequencingRead {
        left: Rc<bam::Record>,
        right: Rc<bam::Record>,
    },
}

impl Evidence {
    pub(crate) fn read_orientation(&self) -> Result<SequenceReadPairOrientation> {
        match self {
            Evidence::SingleEndSequencingRead(read) => read_orientation(read.as_ref()),
            Evidence::PairedEndSequencingRead { left, right } => {
                let left_orient = read_orientation(left.as_ref())?;
                let right_orient = read_orientation(right.as_ref())?;
                if left_orient != right_orient {
                    warn!(
                        "Discordant read orientation in read pair {}, ignoring \
                        orientation information for this read. This can happen if the \
                        read mapper does not annotate mate positions of supplementary \
                        alignments as expected. If you believe this is a bug \
                        we would be grateful for an issue report with the problematic \
                        reads at https://github.com/varlociraptor/varlociraptor.",
                        std::str::from_utf8(left.qname()).unwrap(),
                    );
                    Ok(SequenceReadPairOrientation::None)
                } else {
                    read_orientation(left.as_ref())
                }
            }
        }
    }

    pub(crate) fn is_paired(&self) -> bool {
        match self {
            Evidence::SingleEndSequencingRead(read) => read.is_paired(),
            Evidence::PairedEndSequencingRead { left, .. } => left.is_paired(),
        }
    }

    pub(crate) fn softclipped(&self) -> bool {
        match self {
            Evidence::SingleEndSequencingRead(rec) => {
                let cigar = rec.cigar_cached().unwrap();
                cigar.leading_softclips() > 0 || cigar.trailing_softclips() > 0
            }
            Evidence::PairedEndSequencingRead { left, right } => {
                let cigar_left = left.cigar_cached().unwrap();
                let cigar_right = right.cigar_cached().unwrap();
                cigar_left.leading_softclips() > 0
                    || cigar_left.trailing_softclips() > 0
                    || cigar_right.leading_softclips() > 0
                    || cigar_right.trailing_softclips() > 0
            }
        }
    }

    pub(crate) fn len(&self) -> usize {
        match self {
            Evidence::SingleEndSequencingRead(rec) => rec.seq_len(),
            Evidence::PairedEndSequencingRead { left, right } => left.seq_len() + right.seq_len(),
        }
    }

    pub(crate) fn id(&self) -> EvidenceIdentifier {
        match self {
            Evidence::PairedEndSequencingRead { left, .. } => {
                EvidenceIdentifier::Bytes(left.qname().to_owned())
            }
            Evidence::SingleEndSequencingRead(rec) => {
                EvidenceIdentifier::Bytes(rec.qname().to_owned())
            }
        }
    }

    pub(crate) fn alt_loci(&self) -> ExactAltLoci {
        match self {
            Evidence::SingleEndSequencingRead(rec) => ExactAltLoci::from(rec.as_ref()),
            Evidence::PairedEndSequencingRead { left, right } => {
                let mut left = ExactAltLoci::from(left.as_ref());
                left.inner.extend(ExactAltLoci::from(right.as_ref()).inner);
                left
            }
        }
    }
}

impl PartialEq for Evidence {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Evidence::SingleEndSequencingRead(a), Evidence::SingleEndSequencingRead(b)) => {
                a.qname() == b.qname()
            }
            (
                Evidence::PairedEndSequencingRead { left: a, .. },
                Evidence::PairedEndSequencingRead { left: b, .. },
            ) => a.qname() == b.qname(),
            _ => false,
        }
    }
}

impl Hash for Evidence {
    fn hash<H: Hasher>(&self, state: &mut H) {
        match self {
            Evidence::SingleEndSequencingRead(a) => a.qname().hash(state),
            Evidence::PairedEndSequencingRead { left: a, .. } => a.qname().hash(state),
        }
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub(crate) enum EvidenceIdentifier {
    Bytes(Vec<u8>),
    Integer(u32),
}

impl std::fmt::Display for EvidenceIdentifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EvidenceIdentifier::Bytes(id) => write!(f, "{}", str::from_utf8(id).unwrap()),
            EvidenceIdentifier::Integer(id) => write!(f, "{}", id),
        }
    }
}

impl Hash for EvidenceIdentifier {
    fn hash<H: Hasher>(&self, state: &mut H) {
        match self {
            EvidenceIdentifier::Bytes(id) => id.hash(state),
            EvidenceIdentifier::Integer(id) => id.hash(state),
        }
    }
}
