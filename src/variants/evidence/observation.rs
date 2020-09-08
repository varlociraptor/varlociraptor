// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::hash::{Hash, Hasher};
use std::ops::Deref;
use std::rc::Rc;

use anyhow::Result;
use bio::stats::{LogProb, Prob};
use rust_htslib::bam;
use serde::ser::{SerializeStruct, Serializer};
use serde::Serialize;
// use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::bayesian::BayesFactor;
use itertools::Itertools;

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::sample;
use crate::variants::types::Variant;

/// Calculate expected value of sequencing depth, considering mapping quality.
pub(crate) fn expected_depth(obs: &[Observation]) -> u32 {
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

impl Default for Strand {
    fn default() -> Self {
        Strand::None
    }
}

/// Read orientation support for observation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) enum ReadOrientation {
    F1R2,
    F2R1,
    R1F2,
    R2F1,
    F1F2,
    R1R2,
    F2F1,
    R2R1,
    None,
}

impl ReadOrientation {
    pub(crate) fn new(record: &bam::Record) -> Self {
        if record.is_paired() && record.is_proper_pair() && record.tid() == record.mtid() {
            let (is_reverse, is_first_in_template, is_mate_reverse) =
                if record.pos() < record.mpos() {
                    // given record is the left one
                    (
                        record.is_reverse(),
                        record.is_first_in_template(),
                        record.is_mate_reverse(),
                    )
                } else {
                    // given record is the right one
                    (
                        record.is_mate_reverse(),
                        record.is_last_in_template(),
                        record.is_reverse(),
                    )
                };
            match (is_reverse, is_first_in_template, is_mate_reverse) {
                (false, false, false) => ReadOrientation::F2F1,
                (false, false, true) => ReadOrientation::F2R1,
                (false, true, false) => ReadOrientation::F1F2,
                (true, false, false) => ReadOrientation::R2F1,
                (false, true, true) => ReadOrientation::F1R2,
                (true, false, true) => ReadOrientation::R2R1,
                (true, true, false) => ReadOrientation::R1F2,
                (true, true, true) => ReadOrientation::R1R2,
            }
        } else {
            ReadOrientation::None
        }
    }
}

impl Default for ReadOrientation {
    fn default() -> Self {
        ReadOrientation::None
    }
}

/// An observation for or against a variant.
#[derive(Clone, Debug, Builder, Default)]
pub(crate) struct Observation {
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
    /// Strand evidence this observation relies on
    pub(crate) strand: Strand,
    /// Read orientation support this observation relies on
    pub(crate) read_orientation: ReadOrientation,
    /// True if obervation contains softclips
    pub(crate) softclipped: bool,
}

impl ObservationBuilder {
    pub(crate) fn prob_mapping_mismapping(&mut self, prob_mapping: LogProb) -> &mut Self {
        self.prob_mapping(prob_mapping)
            .prob_mismapping(prob_mapping.ln_one_minus_exp())
    }

    pub(crate) fn prob_overlap(&mut self, prob_double_overlap: LogProb) -> &mut Self {
        self.prob_double_overlap(prob_double_overlap)
            .prob_single_overlap(prob_double_overlap.ln_one_minus_exp())
    }
}

impl Observation {
    pub(crate) fn bayes_factor_alt(&self) -> BayesFactor {
        BayesFactor::new(self.prob_alt, self.prob_ref)
    }

    pub(crate) fn is_paired(&self) -> bool {
        self.read_orientation != ReadOrientation::None
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

    /// Adjust prob_mapping in the given pileup to its average (arithmetic mean of the regular probabilities).
    pub(crate) fn adjust_prob_mapping(pileup: &mut [Self]) {
        // METHOD: a pileup can, although consisting largely of uncertain mappers, contain
        // some reads that by accident are certain mappers (although they don't belong here). If those
        // support the variant, they can lead to an artifact.
        // By taking the average MAPQ over the pileup, we make a conservative choice, justified by the fact
        // that MAPQ choice of the mapper is influenced by (a) the locus ambiguity and (b) stochastic noise 
        // driven by sequencing errors and variants at homologous loci. By averaging, we eliminate
        // the stochastic noise. These assumptions are only valid for SNV and MNV loci.
        let adj_mapq = LogProb((LogProb::ln_sum_exp(&pileup.iter().map(|obs| obs.prob_mapping).collect::<Vec<_>>()).exp() / pileup.len() as f64).ln());
        for obs in pileup {
            if obs.prob_mapping > adj_mapq {
                obs.prob_mapping_adj = Some(adj_mapq);
                obs.prob_mismapping_adj = Some(adj_mapq.ln_one_minus_exp());
            }
        }
    }

    /// Remove all non-standard alignments from pileup (softclipped observations, non-standard read orientations).
    pub(crate) fn remove_nonstandard_alignments(pileup: Vec<Self>) -> Vec<Self> {
        /// METHOD: this can be helpful to get cleaner SNV and MNV calls. Support for those should be 
        /// solely driven by standard alignments, that are not clipped and in expected orientation.
        /// Otherwise called SNVs can be artifacts of near SVs.
        pileup.into_iter().filter(|obs| !obs.softclipped && (obs.read_orientation == ReadOrientation::F1R2 || obs.read_orientation == ReadOrientation::F2R1 || obs.read_orientation == ReadOrientation::None)).collect()
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
                    .prob_overlap(LogProb::ln_zero()) // no double overlap possible (TODO: check this!)
                    .strand(allele_support.strand())
                    .read_orientation(evidence.read_orientation())
                    .softclipped(evidence.softclipped())
                    .build()
                    .unwrap();
                Some(obs)
            }
            _ => None,
        })
    }
}

pub(crate) trait Evidence {
    fn read_orientation(&self) -> ReadOrientation;

    fn softclipped(&self) -> bool;
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
    fn read_orientation(&self) -> ReadOrientation {
        // Single end evidence can just mean that we only need to consider each read alone,
        // although they are paired. Hence we can still check for read orientation.
        ReadOrientation::new(self.inner.as_ref())
    }

    fn softclipped(&self) -> bool {
        let cigar = self.cigar_cached().unwrap();
        cigar.leading_softclips() > 0 || cigar.trailing_softclips() > 0
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
    fn read_orientation(&self) -> ReadOrientation {
        match self {
            PairedEndEvidence::SingleEnd(_) => ReadOrientation::None,
            PairedEndEvidence::PairedEnd { left, .. } => ReadOrientation::new(left.as_ref()),
        }
    }

    fn softclipped(&self) -> bool {
        match self {
            PairedEndEvidence::SingleEnd(rec) => {
                let cigar = rec.cigar_cached().unwrap();
                cigar.leading_softclips() > 0 || cigar.trailing_softclips() > 0
            },
            PairedEndEvidence::PairedEnd { left, right } => {
                let cigar_left = left.cigar_cached().unwrap();
                let cigar_right = right.cigar_cached().unwrap();
                cigar_left.leading_softclips() > 0 || cigar_left.trailing_softclips() > 0 || cigar_right.leading_softclips() > 0 || cigar_right.trailing_softclips() > 0
            },
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
