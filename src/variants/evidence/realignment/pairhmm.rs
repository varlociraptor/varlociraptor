// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;
use std::fmt::Debug;
use std::ops::Range;

use std::sync::Arc;

use bio::stats::pairhmm;
use bio::stats::{LogProb, Prob};
use num_traits::Zero;
use rust_htslib::bam;

use crate::variants::evidence::bases::prob_read_base_miscall;
use crate::variants::evidence::realignment::edit_distance::EditDistanceHit;

/// Width of band around alignment with optimal edit distance.
pub(crate) const EDIT_BAND: usize = 4;

lazy_static! {
    static ref PROB_CONFUSION: LogProb = LogProb::from(Prob(0.3333));
}

pub(crate) trait RefBaseEmission {
    fn ref_base(&self, i: usize) -> u8;

    fn ref_offset(&self) -> usize;

    fn ref_end(&self) -> usize;

    fn set_ref_offset(&mut self, value: usize);

    fn set_ref_end(&mut self, value: usize);

    /// Reference area that is altered by the variant.
    /// Can return None if not applicable (default).
    fn variant_homopolymer_ref_range(&self) -> Option<Range<u64>>;

    fn shrink_to_hit(&mut self, hit: &EditDistanceHit) {
        self.set_ref_end(cmp::min(
            self.ref_offset() + hit.end() + EDIT_BAND,
            self.ref_end(),
        ));
        self.set_ref_offset(self.ref_offset() + hit.start().saturating_sub(EDIT_BAND));
    }

    fn len_x(&self) -> usize;
}

pub(crate) trait VariantEmission {
    fn is_homopolymer_indel(&self) -> bool;
}

#[macro_export]
macro_rules! default_ref_base_emission {
    () => {
        fn ref_offset(&self) -> usize {
            self.ref_offset
        }

        fn ref_end(&self) -> usize {
            self.ref_end
        }

        fn set_ref_offset(&mut self, value: usize) {
            self.ref_offset = value;
        }

        fn set_ref_end(&mut self, value: usize) {
            self.ref_end = value;
        }
    };
}

/// Gap parameters for PairHMM.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct GapParams {
    #[serde(deserialize_with = "parse_float_or_null")]
    pub(crate) prob_insertion_artifact: LogProb,
    #[serde(deserialize_with = "parse_float_or_null")]
    pub(crate) prob_deletion_artifact: LogProb,
    #[serde(deserialize_with = "parse_float_or_null")]
    pub(crate) prob_insertion_extend_artifact: LogProb,
    #[serde(deserialize_with = "parse_float_or_null")]
    pub(crate) prob_deletion_extend_artifact: LogProb,
}

impl Default for GapParams {
    fn default() -> Self {
        Self {
            prob_insertion_artifact: LogProb::from(Prob(2.8e-6)),
            prob_deletion_artifact: LogProb::from(Prob(5.1e-6)),
            prob_insertion_extend_artifact: LogProb::zero(),
            prob_deletion_extend_artifact: LogProb::zero(),
        }
    }
}

fn parse_float_or_null<'de, D>(d: D) -> Result<LogProb, D::Error>
where
    D: serde::Deserializer<'de>,
{
    serde::Deserialize::deserialize(d)
        .map(|x: Option<f64>| x.map(|v| LogProb(v)).unwrap_or(LogProb::zero()))
}

fn parse_vec_of_float_or_null<'de, D>(d: D) -> Result<Vec<LogProb>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    serde::Deserialize::deserialize(d).map(|values: Vec<Option<f64>>| {
        values
            .into_iter()
            .map(|x| x.map(|v| LogProb(v)).unwrap_or(LogProb::zero()))
            .collect()
    })
}

impl pairhmm::GapParameters for GapParams {
    #[inline]
    fn prob_gap_x(&self) -> LogProb {
        self.prob_insertion_artifact
    }

    #[inline]
    fn prob_gap_y(&self) -> LogProb {
        self.prob_deletion_artifact
    }

    #[inline]
    fn prob_gap_x_extend(&self) -> LogProb {
        self.prob_insertion_extend_artifact
    }

    #[inline]
    fn prob_gap_y_extend(&self) -> LogProb {
        self.prob_deletion_extend_artifact
    }
}

impl pairhmm::StartEndGapParameters for GapParams {
    /// Semiglobal alignment: return true.
    #[inline]
    fn free_start_gap_x(&self) -> bool {
        true
    }

    /// Semiglobal alignment: return true.
    #[inline]
    fn free_end_gap_x(&self) -> bool {
        true
    }

    /// Semiglobal alignment: return 1.0.
    #[inline]
    fn prob_start_gap_x(&self, _: usize) -> LogProb {
        LogProb::ln_one()
    }
}

/// Hop parameters for HomopolyPairHMM.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct HopParams {
    #[serde(deserialize_with = "parse_vec_of_float_or_null")]
    pub(crate) prob_seq_homopolymer: Vec<LogProb>,
    #[serde(deserialize_with = "parse_vec_of_float_or_null")]
    pub(crate) prob_ref_homopolymer: Vec<LogProb>,
    #[serde(deserialize_with = "parse_vec_of_float_or_null")]
    pub(crate) prob_seq_extend_homopolymer: Vec<LogProb>,
    #[serde(deserialize_with = "parse_vec_of_float_or_null")]
    pub(crate) prob_ref_extend_homopolymer: Vec<LogProb>,
}

impl Default for HopParams {
    fn default() -> Self {
        Self {
            prob_seq_homopolymer: vec![
                LogProb::zero(),
                LogProb::zero(),
                LogProb::zero(),
                LogProb::zero(),
            ],
            prob_ref_homopolymer: vec![
                LogProb::zero(),
                LogProb::zero(),
                LogProb::zero(),
                LogProb::zero(),
            ],
            prob_seq_extend_homopolymer: vec![
                LogProb::zero(),
                LogProb::zero(),
                LogProb::zero(),
                LogProb::zero(),
            ],
            prob_ref_extend_homopolymer: vec![
                LogProb::zero(),
                LogProb::zero(),
                LogProb::zero(),
                LogProb::zero(),
            ],
        }
    }
}

impl pairhmm::BaseSpecificHopParameters for HopParams {
    #[inline]
    fn prob_hop_x_with_base(&self, base: u8) -> LogProb {
        match base.to_ascii_uppercase() {
            b'A' => self.prob_seq_homopolymer[0],
            b'C' => self.prob_seq_homopolymer[1],
            b'G' => self.prob_seq_homopolymer[2],
            b'T' => self.prob_seq_homopolymer[3],
            _ => unreachable!(),
        }
    }

    #[inline]
    fn prob_hop_y_with_base(&self, base: u8) -> LogProb {
        match base.to_ascii_uppercase() {
            b'A' => self.prob_ref_homopolymer[0],
            b'C' => self.prob_ref_homopolymer[1],
            b'G' => self.prob_ref_homopolymer[2],
            b'T' => self.prob_ref_homopolymer[3],
            _ => unreachable!(),
        }
    }

    #[inline]
    fn prob_hop_x_extend_with_base(&self, base: u8) -> LogProb {
        match base.to_ascii_uppercase() {
            b'A' => self.prob_seq_extend_homopolymer[0],
            b'C' => self.prob_seq_extend_homopolymer[1],
            b'G' => self.prob_seq_extend_homopolymer[2],
            b'T' => self.prob_seq_extend_homopolymer[3],
            _ => unreachable!(),
        }
    }

    #[inline]
    fn prob_hop_y_extend_with_base(&self, base: u8) -> LogProb {
        match base.to_ascii_uppercase() {
            b'A' => self.prob_ref_extend_homopolymer[0],
            b'C' => self.prob_ref_extend_homopolymer[1],
            b'G' => self.prob_ref_extend_homopolymer[2],
            b'T' => self.prob_ref_extend_homopolymer[3],
            _ => unreachable!(),
        }
    }
}

#[derive(Getters, new)]
pub(crate) struct ReadVsAlleleEmission<'a> {
    #[getset(get = "pub(crate)")]
    read_emission: &'a ReadEmission<'a>,
    allele_emission: Box<dyn RefBaseVariantEmission>,
}

impl<'a> RefBaseEmission for ReadVsAlleleEmission<'a> {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        self.allele_emission.ref_base(i)
    }

    #[inline]
    fn ref_offset(&self) -> usize {
        self.allele_emission.ref_offset()
    }

    #[inline]
    fn ref_end(&self) -> usize {
        self.allele_emission.ref_end()
    }

    #[inline]
    fn set_ref_offset(&mut self, value: usize) {
        self.allele_emission.set_ref_offset(value)
    }

    #[inline]
    fn set_ref_end(&mut self, value: usize) {
        self.allele_emission.set_ref_end(value)
    }

    #[inline]
    fn variant_homopolymer_ref_range(&self) -> Option<Range<u64>> {
        self.allele_emission.variant_homopolymer_ref_range()
    }

    #[inline]
    fn len_x(&self) -> usize {
        self.allele_emission.len_x()
    }
}

impl<'a> VariantEmission for ReadVsAlleleEmission<'a> {
    #[inline]
    fn is_homopolymer_indel(&self) -> bool {
        self.allele_emission.is_homopolymer_indel()
    }
}

impl<'a> pairhmm::EmissionParameters for ReadVsAlleleEmission<'a> {
    #[inline]
    fn prob_emit_xy(&self, i: usize, j: usize) -> bio::stats::pairhmm::XYEmission {
        let r = self.allele_emission.ref_base(i);
        self.read_emission.prob_match_mismatch(j, r)
    }

    #[inline]
    fn prob_emit_x(&self, _: usize) -> LogProb {
        LogProb::ln_one()
    }

    #[inline]
    fn prob_emit_y(&self, j: usize) -> LogProb {
        self.read_emission.prob_insertion(j)
    }

    #[inline]
    fn len_y(&self) -> usize {
        self.read_emission.read_end() - self.read_emission.read_offset()
    }

    #[inline]
    fn len_x(&self) -> usize {
        self.allele_emission.len_x()
    }
}

impl<'a> bio::stats::pairhmm::Emission for ReadVsAlleleEmission<'a> {
    fn emission_x(&self, i: usize) -> u8 {
        self.allele_emission.ref_base(i)
    }

    fn emission_y(&self, j: usize) -> u8 {
        unsafe {
            self.read_emission
                .read_seq
                .decoded_base_unchecked(self.read_emission.project_j(j))
        }
    }
}

#[derive(Getters)]
#[getset(get = "pub")]
pub(crate) struct ReadEmission<'a> {
    pub(crate) read_seq: bam::record::Seq<'a>,
    any_miscall: Vec<LogProb>,
    no_miscall: Vec<LogProb>,
    read_offset: usize,
    read_end: usize,
}

impl<'a> ReadEmission<'a> {
    pub(crate) fn new(
        read_seq: bam::record::Seq<'a>,
        qual: &[u8],
        read_offset: usize,
        read_end: usize,
    ) -> Self {
        let mut any_miscall = vec![LogProb::ln_zero(); read_end - read_offset];
        let mut no_miscall = any_miscall.clone();
        for (j, j_) in (read_offset..read_end).enumerate() {
            let prob_miscall = prob_read_base_miscall(*unsafe { qual.get_unchecked(j_) });
            any_miscall[j] = prob_miscall;
            no_miscall[j] = prob_miscall.ln_one_minus_exp();
        }
        ReadEmission {
            read_seq,
            any_miscall,
            no_miscall,
            read_offset,
            read_end,
        }
    }

    fn particular_miscall(&self, j: usize) -> LogProb {
        (unsafe { self.any_miscall.get_unchecked(j) }) + *PROB_CONFUSION
    }

    /// Calculate probability of read_base given ref_base.
    pub(crate) fn prob_match_mismatch(&self, j: usize, ref_base: u8) -> pairhmm::XYEmission {
        let read_base = unsafe { self.read_seq.decoded_base_unchecked(self.project_j(j)) };

        if read_base == ref_base.to_ascii_uppercase() {
            pairhmm::XYEmission::Match(*unsafe { self.no_miscall.get_unchecked(j) })
        } else {
            // TODO replace the second term with technology specific confusion matrix
            pairhmm::XYEmission::Mismatch(self.particular_miscall(j))
        }
    }

    pub(crate) fn prob_insertion(&self, j: usize) -> LogProb {
        *unsafe { self.any_miscall.get_unchecked(j) }
    }

    #[inline]
    pub(crate) fn project_j(&self, j: usize) -> usize {
        j + self.read_offset
    }

    /// Calculate probability that none of the bases is miscalled.
    pub(crate) fn certainty_est(&self) -> LogProb {
        self.no_miscall.iter().sum()
    }
}

/// Emission parameters for PairHMM over reference allele.
pub(crate) struct ReferenceEmissionParams {
    pub(crate) ref_seq: Arc<Vec<u8>>,
    pub(crate) ref_offset: usize,
    pub(crate) ref_end: usize,
}

impl RefBaseEmission for ReferenceEmissionParams {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        self.ref_seq[i + self.ref_offset]
    }

    fn variant_homopolymer_ref_range(&self) -> Option<Range<u64>> {
        None
    }

    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset
    }

    default_ref_base_emission!();
}

impl VariantEmission for ReferenceEmissionParams {
    fn is_homopolymer_indel(&self) -> bool {
        false
    }
}

pub(crate) trait RefBaseVariantEmission: RefBaseEmission + VariantEmission {}

impl<T: RefBaseEmission + VariantEmission> RefBaseVariantEmission for T {}
