// Copyright 2016-2022 Johannes Köster, David Lähnemann, Till Hartmann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use rust_htslib::bam::{FetchDefinition, Read};
use std::cmp;
use std::collections::HashMap;
use std::convert::TryFrom;
use std::f64;
use std::ops::AddAssign;
use std::path::Path;
use std::str;
use std::u32;

use anyhow::{anyhow, Result};
use bio::stats::pairhmm::HomopolyPairHMM;
use bio::stats::{LogProb, Prob};
use counter::Counter;
use itertools::Itertools;
use num_traits::Zero;
use ordered_float::NotNan;
use rust_htslib::bam::{self, record::Cigar};
use statrs::distribution::ContinuousCDF;
use statrs::statistics::{Data, Distribution, OrderStatistics};

use crate::reference;
use crate::utils::aux_tag_is_entire_fragment;
use crate::utils::bam_utils::idxstats;
use crate::utils::homopolymers::is_homopolymer_seq;
use crate::utils::homopolymers::{extend_homopolymer_stretch, is_homopolymer_iter};
use crate::utils::SimpleCounter;

pub(crate) const MIN_HOMOPOLYMER_LEN: usize = 2;

use crate::variants::evidence::realignment::pairhmm::{GapParams, HopParams};

pub(crate) const NUM_FRAGMENTS: usize = 1_000_000;

struct BackwardsCompatibility;

impl BackwardsCompatibility {
    fn default_homopolymer_error_model() -> HashMap<i16, f64> {
        let mut model = HashMap::new();
        model.insert(0, 0.9975414130829068);
        model.insert(1, 0.0010076175889726332);
        model.insert(-1, 0.0010076175889726332);
        model.insert(-2, 0.00020152351779452663);
        model.insert(2, 0.00010076175889726332);
        model.insert(3, 5.038_087_944_863_166e-5);
        model.insert(-3, 9.068_558_300_753_699e-5);

        model
    }

    fn default_max_mapq() -> u8 {
        60
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub(crate) struct AlignmentProperties {
    pub(crate) insert_size: Option<InsertSize>,
    pub(crate) max_del_cigar_len: Option<u32>,
    pub(crate) max_ins_cigar_len: Option<u32>,
    pub(crate) frac_max_softclip: Option<f64>,
    pub(crate) max_read_len: u32,
    #[serde(default = "BackwardsCompatibility::default_max_mapq")]
    pub(crate) max_mapq: u8,
    #[serde(default)]
    pub(crate) cigar_counts: Option<CigarStats>,
    #[serde(default)]
    pub(crate) transition_counts: Option<TransitionCounts>,
    #[serde(default)]
    pub(crate) gap_params: GapParams,
    #[serde(default)]
    pub(crate) hop_params: HopParams,
    #[serde(default = "BackwardsCompatibility::default_homopolymer_error_model")]
    pub(crate) wildtype_homopolymer_error_model: HashMap<i16, f64>,
    #[serde(default)]
    initial: bool,
    #[serde(skip, default)]
    epsilon_gap: f64,
}

impl AlignmentProperties {
    /// Update maximum observed cigar operation lengths. Return whether any D, I, S, or H operation
    /// was found in the cigar string.
    /// The argument `update_unknown` denotes whether unknown properties shall be updated as well.
    /// This is only desired during initial estimation.
    pub(crate) fn update_max_cigar_ops_len(
        &mut self,
        record: &bam::Record,
        allow_hardclips: bool,
    ) -> (bool, bool) {
        let norm = |j| NotNan::new(j as f64 / record.seq().len() as f64).unwrap();

        let mut is_regular = true;
        let mut has_soft_clip = false;
        let iter = if let Some(cigar) = record.cigar_cached() {
            Box::new(cigar.iter().copied()) as Box<dyn Iterator<Item = Cigar>>
        } else {
            Box::new(iter_cigar(record))
        };
        for c in iter {
            match c {
                Cigar::SoftClip(l) => {
                    let s = norm(l);
                    if let Some(ref mut maxclipfrac) = self.frac_max_softclip {
                        *maxclipfrac = *cmp::max(s, NotNan::new(*maxclipfrac).unwrap())
                    } else if self.initial {
                        self.frac_max_softclip = Some(*s);
                    }
                    is_regular = false;
                    has_soft_clip = true;
                }
                Cigar::Del(l) => {
                    if let Some(ref mut maxlen) = self.max_del_cigar_len {
                        *maxlen = cmp::max(l, *maxlen);
                    } else if self.initial {
                        self.max_del_cigar_len = Some(l);
                    }
                    is_regular = false;
                }
                Cigar::Ins(l) => {
                    if let Some(ref mut maxlen) = self.max_ins_cigar_len {
                        *maxlen = cmp::max(l, *maxlen);
                    } else if self.initial {
                        self.max_ins_cigar_len = Some(l);
                    }
                    is_regular = false;
                }
                Cigar::HardClip(_) if !allow_hardclips => {
                    is_regular = false;
                }
                _ => continue,
            }
        }

        (is_regular, has_soft_clip)
    }

    /// Estimate `AlignmentProperties` from first NUM_FRAGMENTS fragments of bam file.
    /// Only reads that are mapped, not duplicates and where quality checks passed are taken.
    pub(crate) fn estimate(
        path: impl AsRef<Path>,
        omit_insert_size: bool,
        reference_buffer: &mut reference::Buffer,
        num_records: Option<usize>,
        epsilon_gap: f64,
    ) -> Result<Self> {
        // If we do not consider insert size, it is safe to also process hardclipped reads.
        let allow_hardclips = omit_insert_size;

        let mut properties = AlignmentProperties {
            insert_size: None,
            max_del_cigar_len: None,
            max_ins_cigar_len: None,
            frac_max_softclip: None,
            max_read_len: 0,
            max_mapq: 0,
            cigar_counts: Default::default(),
            transition_counts: Default::default(),
            wildtype_homopolymer_error_model: HashMap::new(),
            initial: true,
            gap_params: Default::default(),
            hop_params: Default::default(),
            epsilon_gap,
        };

        #[derive(Debug)]
        struct RecordFlagStats {
            inner: Counter<(&'static str, bool), usize>,
        }
        impl RecordFlagStats {
            fn new() -> Self {
                Self {
                    inner: Counter::new(),
                }
            }
            fn update(&mut self, record: &bam::Record) {
                self.inner[&("mapq_is_zero", record.mapq() == 0)] += 1;
                self.inner[&("is_duplicate", record.is_duplicate())] += 1;
                self.inner[&("is_unmapped", record.is_unmapped())] += 1;
                self.inner[&("is_quality_check_failed", record.is_quality_check_failed())] += 1;
                self.inner[&("is_paired", record.is_paired())] += 1;
                self.inner[&("is_first_in_template", record.is_first_in_template())] += 1;
                self.inner[&("is_mate_unmapped", record.is_mate_unmapped())] += 1;
                self.inner[&("mate_ids_match", (record.tid() == record.mtid()))] += 1;
            }
        }

        struct RecordStats {
            mapq: u8,
            read_len: u32,
            cigar_counts: CigarStats,
            transition_counts: TransitionCounts,
            insert_size: Option<f64>,
        }

        #[derive(Default, Debug)]
        struct AlignmentStats {
            n_reads: usize,
            max_mapq: u8,
            max_read_len: u32,
            n_softclips: u32,
            n_not_usable: u32,
            frac_max_softclip: Option<f64>,
            max_del: Option<u32>,
            max_ins: Option<u32>,
            cigar_counts: CigarStats,
            transition_counts: TransitionCounts,
            tlens: Vec<f64>,
        }

        let mut bam = bam::IndexedReader::from_path(path.as_ref())?;
        bam.fetch(FetchDefinition::All)?;

        // Retrieve number of alignments in the bam file.
        // Use this to estimate the number of alignments needed to estimate the
        // HPHMM's transition probabilities to a certain precision.
        let stats = idxstats(path.as_ref())?;
        let num_alignments = stats
            .iter()
            .map(|(_, (_, n_mapped, _))| n_mapped)
            .sum::<u64>();

        let min_num_alignments_needed =
            Self::estimate_number_of_alignments_for_hphmm_mle_param_estimation(num_alignments);

        // Given the number of alignments needed, calculate the step size for the record iterator.
        // The reason to simply use a fixed step size is that
        // - because it is unknown how alignments are distributed across the genome -
        // it is difficult to define an unbiased sampling strategy + it requires random access.
        let step = (num_alignments as f64 / min_num_alignments_needed as f64)
            .floor()
            .max(1.) as usize;

        let header = bam.header().clone();
        let mut n_records_analysed = 0;
        let mut n_records_skipped = 0;
        let mut record_flag_stats = RecordFlagStats::new();
        let records = bam.records();
        let records = records
            .map(|record| record.unwrap())
            .filter(|record| {
                let skip = record.mapq() == 0
                    || record.is_duplicate()
                    || record.is_quality_check_failed()
                    || record.is_unmapped();
                if skip {
                    n_records_skipped += 1;
                }
                !skip
            })
            .step_by(step)
            .take(min_num_alignments_needed);
        let records = if let Some(num_records) = num_records {
            Box::new(records.take(num_records)) as Box<dyn Iterator<Item = bam::Record>>
        } else {
            Box::new(records)
        };
        let all_stats = records
            .map(|mut record| {
                record_flag_stats.update(&record);
                n_records_analysed += 1;
                // don't cache cigar for "ultralong" reads
                if record.seq_len() <= 50000 {
                    record.cache_cigar();
                }

                let chrom = std::str::from_utf8(header.tid2name(record.tid() as u32)).unwrap();
                let (cigar_counts, transition_counts) = cigar_stats(
                    &record,
                    &reference_buffer.seq(chrom).unwrap(),
                    allow_hardclips,
                );

                let insert_size = {
                    if !cigar_counts.is_not_regular && !omit_insert_size {
                        if record.is_paired()
                            && record.is_first_in_template()
                            && record.tid() == record.mtid()
                            && !record.is_mate_unmapped()
                        {
                            Some(record.insert_size().abs() as f64)
                        } else if !record.is_paired() && aux_tag_is_entire_fragment(&record) {
                            // METHOD: if the read is not paired, but the aux EF tag indicates that it covers the entire
                            // fragment, we can use the length of the read as the insert size.
                            Some((record.cigar_cached().unwrap().end_pos() - record.pos()) as f64)
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                };

                RecordStats {
                    mapq: record.mapq(),
                    read_len: record.seq().len() as u32,
                    cigar_counts,
                    transition_counts,
                    insert_size,
                }
            })
            .fold(AlignmentStats::default(), |mut acc, rs| {
                acc.n_reads += 1;
                acc.max_mapq = acc.max_mapq.max(rs.mapq);
                acc.max_read_len = acc.max_read_len.max(rs.read_len);
                acc.n_softclips += rs.cigar_counts.has_soft_clip as u32;
                acc.n_not_usable += rs.cigar_counts.is_not_regular as u32;
                acc.frac_max_softclip =
                    acc.frac_max_softclip.max(rs.cigar_counts.frac_max_softclip);
                acc.max_ins = OptionMax::max(acc.max_ins, rs.cigar_counts.max_ins);
                acc.max_del = OptionMax::max(acc.max_del, rs.cigar_counts.max_del);
                acc.cigar_counts += rs.cigar_counts;
                acc.transition_counts += rs.transition_counts;
                if let Some(insert_size) = rs.insert_size {
                    acc.tlens.push(insert_size);
                }
                acc
            });

        properties.cigar_counts = Some(all_stats.cigar_counts.clone());
        properties.transition_counts = Some(all_stats.transition_counts.clone());

        properties.wildtype_homopolymer_error_model = properties.wildtype_homopolymer_error_model();

        properties.gap_params = properties.estimate_gap_params().unwrap_or_default();
        properties.hop_params = properties.estimate_hop_params().unwrap_or_default();
        properties.max_read_len = all_stats.max_read_len;
        properties.max_del_cigar_len = all_stats.max_del;
        properties.max_ins_cigar_len = all_stats.max_ins;
        properties.frac_max_softclip = all_stats.frac_max_softclip;
        properties.max_mapq = all_stats.max_mapq;
        properties.max_read_len = all_stats.max_read_len;

        let s = if let Some(n) = num_records {
            format!("in {} alignments", n)
        } else {
            "in any alignments".into()
        };
        if properties.max_del_cigar_len.is_none() {
            warn!(
                "No deletion CIGAR operations found {}. \
                Varlociraptor will be unable to estimate the sampling bias for deletions.",
                s,
            );
        }
        if properties.max_ins_cigar_len.is_none() {
            warn!(
                "No deletion CIGAR operations found {}. \
                Varlociraptor will be unable to estimate the sampling bias for insertions.",
                s
            );
        }
        if properties.frac_max_softclip.is_none() {
            warn!(
                "No softclip CIGAR operations found {}. \
                Varlociraptor will be unable to estimate the sampling bias for larger indels.",
                s,
            )
        }

        // Mark initial estimation as done.
        properties.initial = false;

        if all_stats.tlens.is_empty() {
            warn!(
                "\nFound no records to use for estimating the insert size. Will assume\n\
                single end sequencing data and calculate deletion probabilities without\n\
                considering the insert size.\n\
                \n\
                If your data should be paired end, please consider manually providing\n\
                --alignment-properties, e.g. computed with `samtools stats`. Also,\n\
                the following counts of unusable records might indicate a source of\n\
                this problem:\n\n\
                - I, D, S or H CIGAR operation: {nu}\n\
                - S CIGAR (soft clip, e.g. due to UMIs or adapters): {sc}\n\
                \n\
                In addition, {nr} records were skipped in the estimation for one\n\
                of the following reasons:\n\
                - not paired: {not_paired}\n\
                - not the first segment with regard to the template sequence: {not_first_in_template}\n\
                - mapping quality of 0: {mapq_zero}\n\
                - marked as a duplicate: {duplicate}\n\
                - mate mapped to different template (e.g. different chromosome): {mate_ids_dont_match}\n\
                - failed some quality check according to the 512 SAM flag: {quality_fail}\n\
                - mate unmapped: {mate_unmapped}\n\
                - record unmapped: {record_unmapped}\n",
                nu = all_stats.n_not_usable,
                sc = all_stats.n_softclips,
                nr = n_records_skipped,
                not_paired = record_flag_stats.inner.get(&("is_paired", false)).unwrap_or(&0),
                not_first_in_template = record_flag_stats.inner.get(&("is_first_in_template", false)).unwrap_or(&0),
                mapq_zero = record_flag_stats.inner.get(&("mapq_is_zero", true)).unwrap_or(&0),
                duplicate = record_flag_stats.inner.get(&("is_duplicate", true)).unwrap_or(&0),
                mate_ids_dont_match = record_flag_stats.inner.get(&("mate_ids_match", false)).unwrap_or(&0),
                quality_fail = record_flag_stats.inner.get(&("is_quality_check_failed", true)).unwrap_or(&0),
                mate_unmapped = record_flag_stats.inner.get(&("is_mate_unmapped", true)).unwrap_or(&0),
                record_unmapped = record_flag_stats.inner.get(&("is_unmapped", true)).unwrap_or(&0),
            );
            properties.insert_size = None;
            Ok(properties)
        } else {
            let mut tlens = Data::new(all_stats.tlens);
            let upper = tlens.percentile(95);
            let lower = tlens.percentile(5);
            let valid = Data::new(
                tlens
                    .iter()
                    .cloned()
                    .filter(|l| *l <= upper && *l >= lower)
                    .collect_vec(),
            );

            properties.insert_size = Some(InsertSize {
                mean: valid.iter().sum::<f64>() / valid.len() as f64,
                sd: valid.std_dev().unwrap(),
            });
            Ok(properties)
        }
    }

    fn estimate_number_of_alignments_for_hphmm_mle_param_estimation(num_alignments: u64) -> usize {
        // TODO: Get rough estimate of read length, assume 100 for illumina for now
        // Most reads in the wild are at least 100bp long, so this is a conservative estimate.
        let transitions_per_alignment = 100;

        // The estimation is in terms of the number of transitions,
        // which is roughly the number of alignments times the average read length.
        let num_transitions = num_alignments as f64 * transitions_per_alignment as f64;

        // The target precision for the transition probabilities.
        let precision = 1e-1;

        // The confidence level used during estimation.
        let confidence_level = 0.1;

        // The number of valid transitions (i.e. nonzero transition probability) in the HPHMM.
        let number_of_valid_transitions = 82;

        // Upper confidence_level / number_of_valid_transitions quantile
        // of the chi-squared distribution with 1 degree of freedom.
        let b = statrs::distribution::ChiSquared::new(1.0)
            .expect("Failed constructing ChiSquared distribution with 1 degree of freedom")
            .inverse_cdf(1. - confidence_level / number_of_valid_transitions as f64);

        // Given expected transition probabilities, calculate the number of alignments needed
        // to estimate the transition probabilities *for this specific set of alignments* to a
        // certain precision and confidence level.
        // Generally, we expect rather low transition probabilities for homopolymer errors,
        // and roughly 1/4 for match -> match transitions.
        let num_alignments_needed = [0.25, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5].map(|p| {
            let p_rel = precision * p;
            ((b * num_transitions * p * (1. - p)
                / (p_rel.powi(2) * (num_transitions - 1.) + b * p * (1. - p)))
                / transitions_per_alignment as f64)
                .ceil() as usize
        });
        // The minimum number of alignments that need to be examined to estimate the transition
        // probabilities is the maximum of the estimates for each transition probability above.
        let min_num_alignments_needed = *num_alignments_needed.iter().max().unwrap();
        min_num_alignments_needed
    }
}

#[repr(usize)]
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum State {
    MatchA = 0,
    MatchC = 1,
    MatchG = 2,
    MatchT = 3,
    GapX = 4,
    GapY = 5,
    HopAX = 6,
    HopAY = 7,
    HopCX = 8,
    HopCY = 9,
    HopGX = 10,
    HopGY = 11,
    HopTX = 12,
    HopTY = 13,
    Other = 14,
}

const STATES: [State; 14] = [
    State::MatchA,
    State::MatchC,
    State::MatchG,
    State::MatchT,
    State::GapX,
    State::GapY,
    State::HopAX,
    State::HopAY,
    State::HopCX,
    State::HopCY,
    State::HopGX,
    State::HopGY,
    State::HopTX,
    State::HopTY,
];

impl State {
    fn match_(c: u8) -> Self {
        use State::*;
        match c.to_ascii_uppercase() {
            b'A' => MatchA,
            b'C' => MatchC,
            b'G' => MatchG,
            b'T' => MatchT,
            _ => Other,
        }
    }

    fn hop_x(c: u8) -> Self {
        use State::*;
        match c.to_ascii_uppercase() {
            b'A' => HopAX,
            b'C' => HopCX,
            b'G' => HopGX,
            b'T' => HopTX,
            _ => Other,
        }
    }

    fn hop_y(c: u8) -> Self {
        use State::*;
        match c.to_ascii_uppercase() {
            b'A' => HopAY,
            b'C' => HopCY,
            b'G' => HopGY,
            b'T' => HopTY,
            _ => Other,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct TransitionCounts {
    // a State x State matrix with each cell containing the number of transitions
    // from state x to state y
    inner: ndarray::Array2<usize>,
}

impl Default for TransitionCounts {
    fn default() -> Self {
        Self {
            inner: ndarray::Array2::default((16, 16)),
        }
    }
}

impl AddAssign for TransitionCounts {
    fn add_assign(&mut self, rhs: Self) {
        self.inner += &rhs.inner;
    }
}

impl TransitionCounts {
    pub(crate) fn increment(&mut self, from: State, to: State) {
        self.increment_by(from, to, 1);
    }

    pub(crate) fn increment_by(&mut self, from: State, to: State, value: usize) {
        self.inner[(from as usize, to as usize)] += value;
    }

    pub(crate) fn count(&self, from: State, to: State) -> usize {
        self.inner[(from as usize, to as usize)]
    }
}

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub(crate) struct CigarStats {
    pub(crate) gap_counts: SimpleCounter<isize>,
    pub(crate) hop_counts: HashMap<u8, SimpleCounter<(usize, usize)>>,
    pub(crate) match_counts: SimpleCounter<u32>,
    pub(crate) num_match_bases: u64,
    pub(crate) num_ins_bases: u64,
    pub(crate) num_del_bases: u64,
    pub(crate) is_not_regular: bool,
    pub(crate) has_soft_clip: bool,
    pub(crate) frac_max_softclip: Option<f64>,
    pub(crate) max_del: Option<u32>,
    pub(crate) max_ins: Option<u32>,
}

impl AddAssign for CigarStats {
    fn add_assign(&mut self, rhs: Self) {
        self.num_match_bases += rhs.num_match_bases;
        self.num_ins_bases += rhs.num_ins_bases;
        self.num_del_bases += rhs.num_del_bases;
        self.gap_counts += rhs.gap_counts;
        for (base, counts) in rhs.hop_counts {
            *self.hop_counts.entry(base).or_insert_with(Default::default) += counts;
        }
        self.match_counts += rhs.match_counts;
    }
}

fn iter_cigar(record: &bam::Record) -> impl Iterator<Item = Cigar> + '_ {
    record.raw_cigar().iter().map(|&c| {
        let len = c >> 4;
        match c & 0b1111 {
            0 => Cigar::Match(len),
            1 => Cigar::Ins(len),
            2 => Cigar::Del(len),
            3 => Cigar::RefSkip(len),
            4 => Cigar::SoftClip(len),
            5 => Cigar::HardClip(len),
            6 => Cigar::Pad(len),
            7 => Cigar::Equal(len),
            8 => Cigar::Diff(len),
            _ => panic!("Unexpected cigar operation"),
        }
    })
}

/// Count the following operations:
/// - gaps, both insertions and deletions, by length
/// - homopolymer runs, both insertions and deletions, by length and base
/// - matches, independent of the fact whether they are exact or contain substitutions
fn cigar_stats(
    record: &bam::Record,
    refseq: &[u8],
    allow_hardclips: bool,
) -> (CigarStats, TransitionCounts) {
    let norm = |j| NotNan::new(j as f64 / record.seq().len() as f64).unwrap();

    let mut cigar_stats = CigarStats::default();
    let mut transition_stats = TransitionCounts::default();
    let qseq = record.seq();
    let mut qpos = 0usize;
    let mut rpos = record.pos() as usize;

    let iter = if let Some(cigar) = record.cigar_cached() {
        Box::new(cigar.iter().copied()) as Box<dyn Iterator<Item = Cigar>>
    } else {
        Box::new(iter_cigar(record))
    };

    use State::{GapX, GapY};
    for c in iter {
        match c {
            Cigar::Del(l) => {
                cigar_stats.max_del = OptionMax::max(cigar_stats.max_del, Some(l));
                cigar_stats.is_not_regular = true;

                cigar_stats.num_del_bases += l as u64;
                let l = l as usize;
                if l < i16::MAX as usize {
                    let base = refseq[rpos];
                    let is_homopolymer = is_homopolymer_seq(&refseq[rpos..rpos + l]);
                    if is_homopolymer {
                        let mut len = l;
                        if rpos + l < refseq.len() {
                            len += extend_homopolymer_stretch(
                                refseq[rpos],
                                &mut refseq[rpos + l..].iter(),
                            )
                        }
                        if rpos > 1 {
                            len += extend_homopolymer_stretch(
                                refseq[rpos],
                                &mut refseq[..rpos - 1].iter().rev(),
                            )
                        }
                        if len >= MIN_HOMOPOLYMER_LEN {
                            let ms = State::match_(base);
                            let hs = State::hop_x(base);
                            transition_stats.increment_by(ms, ms, l);
                            transition_stats.increment_by(ms, hs, 1);
                            let num_hops = len.saturating_sub(l.saturating_sub(2));
                            transition_stats.increment_by(hs, hs, num_hops);
                            if rpos + len + 1 < refseq.len() {
                                transition_stats
                                    .increment(hs, State::match_(refseq[rpos + len + 1]));
                            }
                            cigar_stats
                                .hop_counts
                                .entry(base)
                                .or_insert_with(SimpleCounter::default)
                                .incr((len, len - l));
                        }
                    }
                    if !is_homopolymer || l == 1 {
                        transition_stats.increment_by(State::match_(base), GapX, 1);
                        transition_stats.increment_by(GapX, GapX, l.saturating_sub(2));
                        if rpos + l + 1 < refseq.len() {
                            transition_stats.increment(GapX, State::match_(refseq[rpos + l + 1]));
                        }
                        cigar_stats.gap_counts.incr(-(isize::try_from(l).unwrap()));
                    }
                }
                rpos += l as usize;
            }
            Cigar::Ins(l) => {
                cigar_stats.max_ins = OptionMax::max(cigar_stats.max_ins, Some(l));
                cigar_stats.is_not_regular = true;

                cigar_stats.num_ins_bases += l as u64;
                let l = l as usize;
                if l < i16::MAX as usize {
                    let base = if refseq[rpos].to_ascii_uppercase() == qseq[qpos] {
                        refseq[rpos]
                    } else {
                        qseq[qpos]
                    };
                    let is_homopolymer = is_homopolymer_iter((qpos..qpos + l).map(|i| qseq[i]));
                    if is_homopolymer {
                        let mut len =
                            l + extend_homopolymer_stretch(qseq[qpos], &mut refseq[rpos..].iter());
                        if rpos > 0 {
                            len += extend_homopolymer_stretch(
                                qseq[qpos],
                                &mut refseq[..rpos].iter().rev(),
                            );
                        }
                        if len >= MIN_HOMOPOLYMER_LEN {
                            let ms = State::match_(base);
                            let hs = State::hop_y(base);
                            transition_stats.increment_by(ms, ms, l);
                            transition_stats.increment_by(ms, hs, 1);
                            let num_hops = len.saturating_sub(l.saturating_sub(2));
                            transition_stats.increment_by(hs, hs, num_hops);
                            if rpos + 1 < refseq.len() {
                                transition_stats.increment(hs, State::match_(refseq[rpos + 1]));
                            }
                            cigar_stats
                                .hop_counts
                                .entry(base)
                                .or_insert_with(SimpleCounter::default)
                                .incr((len - l, l));
                        }
                    }
                    if !is_homopolymer || l == 1 {
                        transition_stats.increment_by(State::match_(base), GapY, 1);
                        transition_stats.increment_by(GapY, GapY, l.saturating_sub(2));
                        if rpos + l + 1 < refseq.len() {
                            transition_stats.increment(GapY, State::match_(refseq[rpos + l + 1]));
                        }
                        cigar_stats.gap_counts.incr(isize::try_from(l).unwrap());
                    }
                }
                qpos += l as usize;
            }
            Cigar::Match(l) | Cigar::Diff(l) | Cigar::Equal(l) => {
                cigar_stats.num_match_bases += l as u64;
                cigar_stats.match_counts.incr(l);
                let l = l as usize;
                for ((rbase, qbase), stretch) in &refseq[rpos..rpos + l]
                    .iter()
                    .zip((qpos..qpos + l).map(|i| qseq[i]))
                    .group_by(|(c_r, c_q)| (**c_r, *c_q))
                {
                    if rbase.to_ascii_uppercase() == qbase {
                        let len = stretch.count();
                        if len >= MIN_HOMOPOLYMER_LEN {
                            cigar_stats
                                .hop_counts
                                .entry(rbase)
                                .or_insert_with(SimpleCounter::default)
                                .incr((len, len));
                        }
                    }
                }
                for (r1, r2) in refseq[rpos..rpos + l].iter().tuple_windows() {
                    transition_stats.increment(State::match_(*r1), State::match_(*r2));
                }
                qpos += l as usize;
                rpos += l as usize;
            }
            Cigar::SoftClip(l) => {
                let s = norm(l);
                cigar_stats.frac_max_softclip =
                    OptionMaxFloat::max(cigar_stats.frac_max_softclip, Some(*s));
                cigar_stats.is_not_regular = true;
                cigar_stats.has_soft_clip = true;
                qpos += l as usize;
            }
            Cigar::RefSkip(l) => {
                rpos += l as usize;
            }
            Cigar::HardClip(_) if !allow_hardclips => {
                cigar_stats.is_not_regular = true;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => continue,
        }
    }
    (cigar_stats, transition_stats)
}

impl AlignmentProperties {
    pub(crate) fn estimate_gap_params(&self) -> Result<GapParams> {
        if let Some(cigar_counts) = &self.cigar_counts {
            let gaps = &cigar_counts.gap_counts;
            let gap_counts_with_length = |length: isize| {
                gaps.iter()
                    .filter(|(len, _)| **len == length)
                    .map(|(_, c)| c)
                    .sum::<usize>()
            };
            if ![2, -2].iter().any(|i| gap_counts_with_length(*i) > 0) || gaps.total_count() < 1000
            {
                warn!("Insufficient observations for gap parameter estimation, falling back to default gap parameters");
                return Err(anyhow!(
                    "Insufficient observations for gap parameter estimation"
                ));
            }

            let [(del_open, del_extend), (ins_open, ins_extend)] = [-1, 1].map(|sign| {
                let num_gap1 = gap_counts_with_length(sign);
                let num_gap2 = gap_counts_with_length(2 * sign);

                let gap_open = (num_gap1 + num_gap2) as f64
                    / (cigar_counts.num_match_bases + cigar_counts.num_ins_bases) as f64;
                let gap_extend = num_gap2 as f64 / (num_gap1 as f64 + num_gap2 as f64);
                let gap_extend = if gap_extend < self.epsilon_gap {
                    0.
                } else {
                    gap_extend
                };
                (gap_open, gap_extend)
            });

            Ok(GapParams {
                prob_insertion_artifact: LogProb::from(
                    Prob::checked(ins_open).unwrap_or_else(|_| Prob::zero()),
                ),
                prob_deletion_artifact: LogProb::from(
                    Prob::checked(del_open).unwrap_or_else(|_| Prob::zero()),
                ),
                prob_insertion_extend_artifact: LogProb::from(
                    Prob::checked(ins_extend).unwrap_or_else(|_| Prob::zero()),
                ),
                prob_deletion_extend_artifact: LogProb::from(
                    Prob::checked(del_extend).unwrap_or_else(|_| Prob::zero()),
                ),
            })
        } else {
            warn!("No cigar operations counted, falling back to default gap parameters.");
            Err(anyhow!(
                "Insufficient observations for gap parameter estimation"
            ))
        }
    }

    pub(crate) fn estimate_hop_params(&self) -> Result<HopParams> {
        if let Some(transition_counts) = &self.transition_counts {
            let mut insufficient_counts = false;
            let (
                (prob_seq_homopolymer, prob_ref_homopolymer),
                (prob_seq_extend_homopolymer, prob_ref_extend_homopolymer),
            ): ((Vec<_>, Vec<_>), (Vec<_>, Vec<_>)) = [b'A', b'C', b'G', b'T']
                .iter()
                .map(|&base| {
                    // METHOD: only look at homopolymer events with a length difference of 1 (either ins or del) since these are very likely always artifacts
                    let match_ = State::match_(base);
                    let hop_x = State::hop_x(base);
                    let hop_y = State::hop_y(base);

                    let [prob_ins, prob_del] = [hop_x, hop_y].map(|hop| {
                        let start_hop = transition_counts.count(match_, hop);
                        let extend_hop = transition_counts.count(hop, hop);
                        if start_hop + extend_hop < 100 {
                            insufficient_counts |= true;
                        }
                        let from_previous = STATES
                            .map(|s| transition_counts.count(match_, s))
                            .iter()
                            .sum::<usize>();
                        let prob = (start_hop + extend_hop) as f64 / from_previous as f64;
                        LogProb::from(Prob::checked(prob).unwrap_or_else(|_| Prob::zero()))
                    });

                    ((prob_ins, prob_del), (prob_ins, prob_del))
                })
                .unzip();

            if insufficient_counts {
                warn!("Insufficient observations for hop parameter estimation, falling back to default hop parameters");
                return Err(anyhow!(
                    "Insufficient observations for hop parameter estimation"
                ));
            }
            Ok(HopParams {
                prob_seq_homopolymer,
                prob_ref_homopolymer,
                prob_seq_extend_homopolymer,
                prob_ref_extend_homopolymer,
            })
        } else {
            warn!("Insufficient observations for hop parameter estimation, falling back to default hop parameters");
            Err(anyhow!(
                "No cigar operations counted, falling back to default hop parameters"
            ))
        }
    }

    pub(crate) fn wildtype_homopolymer_error_model(&self) -> HashMap<i16, f64> {
        if let Some(cigar_counts) = &self.cigar_counts {
            let n = cigar_counts
                .hop_counts
                .values()
                .flat_map(|counter| counter.values())
                .filter(|count| **count >= 10)
                .sum::<usize>() as f64;
            // group by length, i.e. discard base information
            let counts = cigar_counts
                .hop_counts
                .values()
                .flat_map(|counter| counter.iter())
                .map(|(k, v)| {
                    (
                        (isize::try_from(k.0).unwrap() - isize::try_from(k.1).unwrap())
                            .clamp(i16::MIN.into(), i16::MAX.into()) as i16,
                        *v,
                    )
                })
                .collect_vec();
            let grouped = counts
                .iter()
                .sorted_unstable_by_key(|(length, _)| length)
                .group_by(|(length, _)| length);
            grouped
                .into_iter()
                .map(|(length, group)| (length, group.map(|(_, count)| count).sum::<usize>()))
                .map(|(length, count)| (*length, count as f64 / n))
                .collect()
        } else {
            self.wildtype_homopolymer_error_model.clone()
        }
    }
}

/// Expected insert size in terms of mean and standard deviation.
/// This should be estimated from unsorted(!) bam files to avoid positional biases.
#[derive(Copy, Clone, Debug, Default, Serialize, Deserialize)]
pub(crate) struct InsertSize {
    pub(crate) mean: f64,
    pub(crate) sd: f64,
}

trait OptionMaxFloat<F> {
    fn max(self, other: Option<F>) -> Option<F>;
}

trait OptionMax<V> {
    fn max(self, other: Option<V>) -> Option<V>;
}

impl OptionMaxFloat<f64> for Option<f64> {
    fn max(self, other: Option<f64>) -> Option<f64> {
        match (self, other) {
            (None, None) => None,
            (Some(a), None) => Some(a),
            (None, Some(b)) => Some(b),
            (Some(a), Some(b)) => Some(a.max(b)),
        }
    }
}

impl<V: Ord> OptionMax<V> for Option<V> {
    fn max(self, other: Option<V>) -> Option<V> {
        match (self, other) {
            (None, None) => None,
            (Some(a), None) => Some(a),
            (None, Some(b)) => Some(b),
            (Some(a), Some(b)) => Some(a.max(b)),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::utils::bam_utils::idxstats;
    use bio::io::fasta;

    use super::*;

    fn reference_buffer() -> reference::Buffer {
        reference::Buffer::new(
            fasta::IndexedReader::from_file(&"tests/resources/chr10.fa").unwrap(),
            1,
        )
    }

    #[test]
    fn test_estimate() {
        let path = "tests/resources/tumor-first30000.bam";
        let mut reference_buffer = reference_buffer();

        let props = AlignmentProperties::estimate(
            path,
            false,
            &mut reference_buffer,
            Some(NUM_FRAGMENTS),
            0.,
        )
        .unwrap();
        println!("{:?}", props);

        if let Some(isize) = props.insert_size {
            assert_relative_eq!(isize.mean.round(), 312.0);
            assert_relative_eq!(isize.sd.round(), 12.0);
        } else {
            panic!("test_estimate(): props.insert_size was None. Something is wrong, here.");
        }
        assert_eq!(props.max_del_cigar_len, Some(30));
        assert_eq!(props.max_ins_cigar_len, Some(12));
        assert_eq!(props.frac_max_softclip, Some(0.63));
    }

    #[test]
    fn test_estimate_all_reads_have_short_clips() {
        let path = "tests/resources/tumor-first30000.reads_with_soft_clips.bam";
        let mut reference_buffer = reference_buffer();

        let props = AlignmentProperties::estimate(
            path,
            false,
            &mut reference_buffer,
            Some(NUM_FRAGMENTS),
            0.,
        )
        .unwrap();
        println!("{:?}", props);

        assert!(props.insert_size.is_none());
        assert_eq!(props.max_del_cigar_len, Some(2));
        assert_eq!(props.max_ins_cigar_len, Some(4));
        assert_eq!(props.frac_max_softclip, Some(0.63));
    }

    #[test]
    fn test_estimate_all_reads_single_end() {
        // this file contains only single-ended reads (artificially made single-ended with awk)
        let path = "tests/resources/tumor-first30000.bunch_of_reads_made_single_ended.bam";
        let mut reference_buffer = reference_buffer();

        let props = AlignmentProperties::estimate(
            path,
            false,
            &mut reference_buffer,
            Some(NUM_FRAGMENTS),
            0.,
        )
        .unwrap();
        println!("{:?}", props);

        assert!(props.insert_size.is_none());
        assert_eq!(props.max_del_cigar_len, None);
        assert_eq!(props.max_ins_cigar_len, None);
        assert_eq!(props.frac_max_softclip, Some(0.03));
    }
}
