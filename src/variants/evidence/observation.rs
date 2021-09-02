// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::char;
use std::hash::{Hash, Hasher};
use std::ops;
use std::ops::Deref;
use std::rc::Rc;
use std::str;

use anyhow::Result;
use bio::stats::LogProb;
use bio_types::sequence::SequenceReadPairOrientation;
use counter::Counter;
use rust_htslib::bam;
use serde::ser::{SerializeStruct, Serializer};
use serde::Serialize;
// use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::bayesian::BayesFactor;
use itertools::Itertools;

use crate::errors::{self, Error};
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::utils;
use crate::variants::sample;
use crate::variants::types::Variant;

/// Calculate expected value of sequencing depth, considering mapping quality.
pub(crate) fn expected_depth(obs: &[Observation<ReadPosition>]) -> u32 {
    LogProb::ln_sum_exp(&obs.iter().map(|o| o.prob_mapping).collect_vec())
        .exp()
        .round() as u32
}

/// Strand support for observation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) enum Strand {
    Forward,
    Reverse,
    Both,
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

impl Default for Strand {
    fn default() -> Self {
        Strand::None
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) enum ReadPosition {
    Major,
    Some,
}

impl Default for ReadPosition {
    fn default() -> Self {
        ReadPosition::Some
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

/// An observation for or against a variant.
#[derive(Clone, Debug, Builder, Default)]
pub(crate) struct Observation<P = Option<u32>>
where
    P: Clone,
{
    /// Posterior probability that the read/read-pair has been mapped correctly (1 - MAPQ).
    prob_mapping: LogProb,
    /// Posterior probability that the read/read-pair has been mapped incorrectly (MAPQ).
    prob_mismapping: LogProb,
    /// Posterior probability that the read/read-pair has been mapped correctly (1 - MAPQ), adjusted form.
    #[builder(private, default = "None")]
    prob_mapping_adj: Option<LogProb>,
    /// Posterior probability that the read/read-pair has been mapped incorrectly (MAPQ), adjusted form.
    #[builder(private, default = "None")]
    prob_mismapping_adj: Option<LogProb>,
    /// Probability that the read/read-pair comes from the alternative allele.
    pub(crate) prob_alt: LogProb,
    /// Probability that the read/read-pair comes from the reference allele.
    pub(crate) prob_ref: LogProb,
    /// Probability that the read/read-pair comes from an unknown allele at an unknown true
    /// locus (in case it is mismapped). This should usually be set as the product of the maxima
    /// of prob_ref and prob_alt per read.
    pub(crate) prob_missed_allele: LogProb,
    /// Probability to sample the alt allele
    pub(crate) prob_sample_alt: LogProb,
    /// Probability to overlap with both strands
    #[builder(private)]
    pub(crate) prob_double_overlap: LogProb,
    /// Probability to overlap with one strand only (1-prob_double_overlap)
    #[builder(private)]
    pub(crate) prob_single_overlap: LogProb,
    pub(crate) prob_hit_base: LogProb,
    /// Strand evidence this observation relies on
    pub(crate) strand: Strand,
    /// Read orientation support this observation relies on
    pub(crate) read_orientation: SequenceReadPairOrientation,
    /// True if obervation contains s
    pub(crate) softclipped: bool,
    pub(crate) paired: bool,
    /// Read position of the variant in the read (for SNV and MNV)
    pub(crate) read_position: P,
    /// Whether the read contains indel operations agains the alt allele
    pub(crate) has_alt_indel_operations: bool,
}

impl<P: Clone> ObservationBuilder<P> {
    pub(crate) fn prob_mapping_mismapping(&mut self, prob_mapping: LogProb) -> &mut Self {
        self.prob_mapping(prob_mapping)
            .prob_mismapping(prob_mapping.ln_one_minus_exp())
    }

    pub(crate) fn prob_overlap(&mut self, prob_double_overlap: LogProb) -> &mut Self {
        self.prob_double_overlap(prob_double_overlap)
            .prob_single_overlap(prob_double_overlap.ln_one_minus_exp())
    }
}

impl Observation<Option<u32>> {
    pub(crate) fn process(&self, major_read_position: Option<u32>) -> Observation<ReadPosition> {
        Observation {
            prob_mapping: self.prob_mapping,
            prob_mismapping: self.prob_mismapping,
            prob_mapping_adj: self.prob_mapping_adj,
            prob_mismapping_adj: self.prob_mismapping_adj,
            prob_alt: self.prob_alt,
            prob_ref: self.prob_ref,
            prob_missed_allele: self.prob_missed_allele,
            prob_sample_alt: self.prob_sample_alt,
            prob_double_overlap: self.prob_double_overlap,
            prob_single_overlap: self.prob_single_overlap,
            prob_hit_base: self.prob_hit_base,
            strand: self.strand,
            read_orientation: self.read_orientation,
            softclipped: self.softclipped,
            paired: self.paired,
            read_position: self.read_position.map_or(ReadPosition::Some, |pos| {
                if let Some(major_pos) = major_read_position {
                    if pos == major_pos {
                        ReadPosition::Major
                    } else {
                        ReadPosition::Some
                    }
                } else {
                    ReadPosition::Some
                }
            }),
            has_alt_indel_operations: self.has_alt_indel_operations,
        }
    }
}

impl<P: Clone> Observation<P> {
    pub(crate) fn bayes_factor_alt(&self) -> BayesFactor {
        BayesFactor::new(self.prob_alt, self.prob_ref)
    }

    pub(crate) fn prob_mapping_orig(&self) -> LogProb {
        self.prob_mapping
    }

    pub(crate) fn prob_mapping(&self) -> LogProb {
        self.prob_mapping_adj.unwrap_or(self.prob_mapping)
    }

    pub(crate) fn prob_mismapping(&self) -> LogProb {
        self.prob_mismapping_adj.unwrap_or(self.prob_mismapping)
    }

    /// Remove all non-standard alignments from pileup (softclipped observations, non-standard read orientations).
    pub(crate) fn remove_nonstandard_alignments(
        pileup: Vec<Self>,
        omit_read_orientation_bias: bool,
    ) -> Vec<Self> {
        // METHOD: this can be helpful to get cleaner SNV and MNV calls. Support for those should be
        // solely driven by standard alignments, that are in expected orientation.
        // Otherwise called SNVs can be artifacts of near SVs.
        pileup
            .into_iter()
            .filter(|obs| {
                omit_read_orientation_bias
                    || (obs.read_orientation == SequenceReadPairOrientation::F1R2
                        || obs.read_orientation == SequenceReadPairOrientation::F2R1
                        || obs.read_orientation == SequenceReadPairOrientation::None)
            })
            .collect()
    }
}

pub(crate) fn major_read_position(pileup: &[Observation<Option<u32>>]) -> Option<u32> {
    let counter: Counter<_> = pileup.iter().filter_map(|obs| obs.read_position).collect();
    let most_common = counter.most_common();
    if most_common.is_empty() {
        None
    } else {
        Some(most_common[0].0)
    }
}

impl Serialize for Observation {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut s = serializer.serialize_struct("Observation", 3)?;
        s.serialize_field("prob_mapping", &self.prob_mapping)?;
        s.serialize_field("prob_mismapping", &self.prob_mismapping)?;
        s.serialize_field("prob_alt", &self.prob_alt)?;
        s.serialize_field("prob_ref", &self.prob_ref)?;
        s.serialize_field("prob_sample_alt", &self.prob_sample_alt)?;
        s.end()
    }
}

/// Something that can be converted into observations.
pub(crate) trait Observable<E>: Variant<Evidence = E>
where
    E: Evidence + Eq + Hash,
{
    fn extract_observations(
        &self,
        buffer: &mut sample::RecordBuffer,
        alignment_properties: &mut AlignmentProperties,
        max_depth: usize,
    ) -> Result<Vec<Observation>>;

    /// Convert MAPQ (from read mapper) to LogProb for the event that the read maps
    /// correctly.
    fn prob_mapping(&self, evidence: &E) -> LogProb;

    /// Calculate an observation from the given evidence.
    fn evidence_to_observation(
        &self,
        evidence: &E,
        alignment_properties: &AlignmentProperties,
    ) -> Result<Option<Observation>> {
        Ok(match self.allele_support(evidence, alignment_properties)? {
            // METHOD: only consider allele support if it comes either from forward or reverse strand.
            // Unstranded observations (e.g. only insert size), are too unreliable, or do not contain
            // any information (e.g. no overlap).
            Some(allele_support) if allele_support.strand() != Strand::None => {
                let obs = ObservationBuilder::default()
                    .prob_mapping_mismapping(self.prob_mapping(evidence))
                    .prob_alt(allele_support.prob_alt_allele())
                    .prob_ref(allele_support.prob_ref_allele())
                    .prob_sample_alt(self.prob_sample_alt(evidence, alignment_properties))
                    .prob_missed_allele(allele_support.prob_missed_allele())
                    .prob_overlap(if allele_support.strand() == Strand::Both {
                        LogProb::ln_one()
                    } else {
                        LogProb::ln_zero()
                    })
                    .strand(allele_support.strand())
                    .read_orientation(evidence.read_orientation()?)
                    .softclipped(evidence.softclipped())
                    .has_alt_indel_operations(if self.report_indel_operations() {
                        allele_support.has_alt_indel_operations()
                    } else {
                        // METHOD: do not report any operations if the variant chooses to not report them.
                        false
                    })
                    .read_position(allele_support.read_position())
                    .paired(evidence.is_paired())
                    .prob_hit_base(LogProb::ln_one() - LogProb((evidence.len() as f64).ln()))
                    .build()
                    .unwrap();
                Some(obs)
            }
            _ => None,
        })
    }
}

pub(crate) trait Evidence {
    fn read_orientation(&self) -> Result<SequenceReadPairOrientation>;

    fn softclipped(&self) -> bool;

    fn is_paired(&self) -> bool;

    fn len(&self) -> usize;
}

#[derive(new, Clone, Eq, Debug)]
pub(crate) struct SingleEndEvidence {
    inner: Rc<bam::Record>,
}

impl Deref for SingleEndEvidence {
    type Target = bam::Record;

    fn deref(&self) -> &bam::Record {
        self.inner.as_ref()
    }
}

impl Evidence for SingleEndEvidence {
    fn read_orientation(&self) -> Result<SequenceReadPairOrientation> {
        // Single end evidence can just mean that we only need to consider each read alone,
        // although they are paired. Hence we can still check for read orientation.
        read_orientation(self.inner.as_ref())
    }

    fn is_paired(&self) -> bool {
        self.inner.is_paired()
    }

    fn softclipped(&self) -> bool {
        let cigar = self.cigar_cached().unwrap();
        cigar.leading_softclips() > 0 || cigar.trailing_softclips() > 0
    }

    fn len(&self) -> usize {
        self.inner.seq_len()
    }
}

impl PartialEq for SingleEndEvidence {
    fn eq(&self, other: &Self) -> bool {
        self.qname() == other.qname()
    }
}

impl Hash for SingleEndEvidence {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.qname().hash(state);
    }
}

#[derive(Clone, Eq, Debug)]
pub(crate) enum PairedEndEvidence {
    SingleEnd(Rc<bam::Record>),
    PairedEnd {
        left: Rc<bam::Record>,
        right: Rc<bam::Record>,
    },
}

impl Evidence for PairedEndEvidence {
    fn read_orientation(&self) -> Result<SequenceReadPairOrientation> {
        match self {
            PairedEndEvidence::SingleEnd(read) => read_orientation(read.as_ref()),
            PairedEndEvidence::PairedEnd { left, .. } => read_orientation(left.as_ref()),
        }
    }

    fn is_paired(&self) -> bool {
        match self {
            PairedEndEvidence::SingleEnd(read) => read.is_paired(),
            PairedEndEvidence::PairedEnd { left, .. } => left.is_paired(),
        }
    }

    fn softclipped(&self) -> bool {
        match self {
            PairedEndEvidence::SingleEnd(rec) => {
                let cigar = rec.cigar_cached().unwrap();
                cigar.leading_softclips() > 0 || cigar.trailing_softclips() > 0
            }
            PairedEndEvidence::PairedEnd { left, right } => {
                let cigar_left = left.cigar_cached().unwrap();
                let cigar_right = right.cigar_cached().unwrap();
                cigar_left.leading_softclips() > 0
                    || cigar_left.trailing_softclips() > 0
                    || cigar_right.leading_softclips() > 0
                    || cigar_right.trailing_softclips() > 0
            }
        }
    }

    fn len(&self) -> usize {
        match self {
            PairedEndEvidence::SingleEnd(rec) => rec.seq_len(),
            PairedEndEvidence::PairedEnd { left, right } => left.seq_len() + right.seq_len(),
        }
    }
}

impl PartialEq for PairedEndEvidence {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (PairedEndEvidence::SingleEnd(a), PairedEndEvidence::SingleEnd(b)) => {
                a.qname() == b.qname()
            }
            (
                PairedEndEvidence::PairedEnd { left: a, .. },
                PairedEndEvidence::PairedEnd { left: b, .. },
            ) => a.qname() == b.qname(),
            _ => false,
        }
    }
}

impl Hash for PairedEndEvidence {
    fn hash<H: Hasher>(&self, state: &mut H) {
        match self {
            PairedEndEvidence::SingleEnd(a) => a.qname().hash(state),
            PairedEndEvidence::PairedEnd { left: a, .. } => a.qname().hash(state),
        }
    }
}
