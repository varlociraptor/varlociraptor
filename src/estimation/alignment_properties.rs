// Copyright 2016-2022 Johannes Köster, David Lähnemann, Till Hartmann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;
use std::collections::HashMap;
use std::convert::TryFrom;
use std::f64;
use std::ops::AddAssign;
use std::str;
use std::u32;

use crate::estimation::alignment_properties::State::HopAX;
use anyhow::{anyhow, Result};
use bio::stats::pairhmm::HomopolyPairHMM;
use bio::stats::{LogProb, Prob};
use counter::Counter;
use itertools::Itertools;
use num_traits::Zero;
use ordered_float::NotNan;
use rust_htslib::bam::{self, record::Cigar};
use statrs::statistics::{Data, Distribution, OrderStatistics};

use crate::reference;
use crate::utils::aux_tag_is_entire_fragment;
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
    pub(crate) fn estimate<R: bam::Read>(
        bam: &mut R,
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
            // TODO expose step parameter to cli
            .step_by(1);
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
                let cigar_counts = cigar_stats(
                    &record,
                    &reference_buffer.seq(chrom).unwrap(),
                    allow_hardclips,
                );
                let transition_counts = transition_stats(
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

        for i in 0..16 {
            for j in 0..16 {
                eprint!("{}\t", all_stats.transition_counts.inner[(i, j)]);
            }
            eprintln!();
        }
        properties.cigar_counts = Some(all_stats.cigar_counts.clone());

        properties.wildtype_homopolymer_error_model = properties.wildtype_homopolymer_error_model();

        properties.gap_params = properties.estimate_gap_params().unwrap_or_default();
        properties.hop_params = properties.estimate_hop_params().unwrap_or_default();
        #[derive(Serialize)]
        struct Props {
            gap_counts: Vec<(isize, usize)>,
            hop_counts: Vec<(u8, Vec<((usize, usize), usize)>)>,
            transition_counts: TransitionCounts,
        }
        let g_counts = properties.cigar_counts.as_ref().unwrap().gap_counts.clone();
        let h_counts = properties.cigar_counts.as_ref().unwrap().hop_counts.clone();
        serde_json::to_writer_pretty(
            &mut std::io::stdout(),
            &Props {
                gap_counts: g_counts.iter().map(|(x, y)| (*x, *y)).collect(),
                hop_counts: h_counts
                    .into_iter()
                    .map(|(b, v)| (b, v.iter().map(|(x, y)| (*x, *y)).collect::<Vec<_>>()))
                    .collect(),
                transition_counts: all_stats.transition_counts,
            },
        )?;
        std::process::exit(0);

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
    Start = 15,
}

impl State {
    fn supports(&self, c: u8) -> bool {
        use State::*;
        match self {
            MatchA | HopAX | HopAY if c.eq_ignore_ascii_case(&b'A') => true,
            MatchC | HopCX | HopCY if c.eq_ignore_ascii_case(&b'C') => true,
            MatchG | HopGX | HopGY if c.eq_ignore_ascii_case(&b'G') => true,
            MatchT | HopTX | HopTY if c.eq_ignore_ascii_case(&b'T') => true,
            _ => false,
        }
    }

    fn is_match(&self) -> bool {
        use State::*;
        match self {
            MatchA | MatchC | MatchG | MatchT => true,
            _ => false,
        }
    }

    fn is_hop(&self) -> bool {
        use State::*;
        match self {
            HopAX | HopAY | HopCX | HopCY | HopGX | HopGY | HopTX | HopTY => true,
            _ => false,
        }
    }

    fn is_gap(&self) -> bool {
        use State::*;
        match self {
            GapX | GapY => true,
            _ => false,
        }
    }

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
    inner: ndarray::Array2<usize>, // enum State num members/variants == 14, + 1 for start state?
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
}

fn advance_positions(c: Cigar, rpos: &mut usize, qpos: &mut usize) {
    match c {
        Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
            *rpos += n as usize;
            *qpos += n as usize;
        }
        Cigar::Ins(n) => {
            *qpos += n as usize;
        }
        Cigar::Del(n) => {
            *rpos += n as usize;
        }
        Cigar::RefSkip(n) => {
            *rpos += n as usize;
        }
        Cigar::SoftClip(n) => {
            *qpos += n as usize;
        }
        Cigar::HardClip(_) => {}
        Cigar::Pad(_) => {}
    }
}

fn transition_stats(record: &bam::Record, rseq: &[u8], allow_hardclips: bool) -> TransitionCounts {
    let mut transition_stats = TransitionCounts::default();
    let qseq = record.seq();
    let mut qpos = 0usize;
    let mut rpos = record.pos() as usize;
    let iter = iter_cigar(record);
    let mut prev_state = State::Start;
    for (r1, q1) in iter
        .flat_map(|c| {
            advance_positions(c, &mut rpos, &mut qpos);
            match c {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    let n = n as usize;
                    (0..n)
                        .map(|i| (rseq[rpos + i - n], qseq[qpos + i - n]))
                        .collect_vec()
                }
                Cigar::Ins(n) => {
                    let n = n as usize;
                    (0..n).map(|i| (b'-', qseq[qpos + i - n])).collect_vec()
                }
                Cigar::Del(n) => {
                    let n = n as usize;
                    (0..n).map(|i| (rseq[rpos + i - n], b'-')).collect_vec()
                }
                _ => {
                    vec![(b'-', b'-')]
                }
            }
        })
    {
        use State::*;
        match (r1, q1) {
            (b'-', b'-') => {
                prev_state = Start;
            }
            (b'-', q) => match prev_state {
                s if s.supports(q) && (s.is_hop() || s.is_match()) => {
                    prev_state = State::hop_x(q);
                    transition_stats.increment(s, prev_state)
                }
                s => {
                    prev_state = GapX;
                    transition_stats.increment(s, prev_state);
                }
            },
            (r, b'-') => match prev_state {
                s if s.supports(r) && (s.is_hop() || s.is_match()) => {
                    prev_state = State::hop_y(r);
                    transition_stats.increment(s, prev_state)
                }
                s => {
                    prev_state = GapY;
                    transition_stats.increment(s, prev_state);
                }
            },
            (r, _q) => {
                let s = prev_state.clone();
                prev_state = State::match_(r);
                transition_stats.increment(s, prev_state);
            }
            _ => unreachable!(),
        }
    }

    transition_stats
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
fn cigar_stats(record: &bam::Record, refseq: &[u8], allow_hardclips: bool) -> CigarStats {
    let norm = |j| NotNan::new(j as f64 / record.seq().len() as f64).unwrap();

    let mut cigar_stats = CigarStats::default();
    let qseq = record.seq();
    let mut qpos = 0usize;
    let mut rpos = record.pos() as usize;

    let iter = if let Some(cigar) = record.cigar_cached() {
        Box::new(cigar.iter().copied()) as Box<dyn Iterator<Item = Cigar>>
    } else {
        Box::new(iter_cigar(record))
    };
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
                            cigar_stats
                                .hop_counts
                                .entry(base)
                                .or_insert_with(SimpleCounter::default)
                                .incr((len, len - l));
                        }
                    }
                    if !is_homopolymer || l == 1 {
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
                            cigar_stats
                                .hop_counts
                                .entry(base)
                                .or_insert_with(SimpleCounter::default)
                                .incr((len - l, l));
                        }
                    }
                    if !is_homopolymer || l == 1 {
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
    cigar_stats
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
        if let Some(cigar_counts) = &self.cigar_counts {
            let mut insufficient_counts = false;
            let empty = SimpleCounter::default();
            let counts = |base| {
                cigar_counts
                    .hop_counts
                    .get(&base)
                    .unwrap_or(&empty)
                    .iter()
                    .map(|(lengths, count)| (*lengths, *count))
            };
            let (
                (prob_seq_homopolymer, prob_ref_homopolymer),
                (prob_seq_extend_homopolymer, prob_ref_extend_homopolymer),
            ): ((Vec<_>, Vec<_>), (Vec<_>, Vec<_>)) = [b'A', b'C', b'G', b'T']
                .iter()
                .map(|base| {
                    let events = |predicate: fn((usize, usize)) -> bool| {
                        counts(*base).filter(move |(l, _)| predicate(*l))
                    };
                    // METHOD: only look at homopolymer events with a length difference of 1 (either ins or del) since these are very likely always artifacts
                    let dels = events(|(r, q)| r == q + 1).collect_vec();
                    let num_diff_1_del_events = dels.iter().map(|(_, c)| c).sum::<usize>();
                    let num_diff_1_del_bases = dels.iter().map(|((r, _), c)| r * c).sum::<usize>();

                    let ins = events(|(r, q)| q == r + 1).collect_vec();
                    let num_diff_1_ins_events = ins.iter().map(|(_, c)| c).sum::<usize>();
                    let num_diff_1_ins_bases = ins.iter().map(|((_, q), c)| q * c).sum::<usize>();

                    let prob_del_hop_start_or_extend =
                        num_diff_1_del_events as f64 / num_diff_1_del_bases as f64;
                    let prob_ins_hop_start_or_extend =
                        num_diff_1_ins_events as f64 / num_diff_1_ins_bases as f64;

                    if num_diff_1_del_events < 100 || num_diff_1_ins_events < 100 {
                        insufficient_counts |= true;
                    }

                    let prob_ins = LogProb::from(
                        Prob::checked(prob_ins_hop_start_or_extend)
                            .unwrap_or_else(|_| Prob::zero()),
                    );
                    let prob_del = LogProb::from(
                        Prob::checked(prob_del_hop_start_or_extend)
                            .unwrap_or_else(|_| Prob::zero()),
                    );
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
        let mut bam = bam::Reader::from_path("tests/resources/tumor-first30000.bam").unwrap();
        let mut reference_buffer = reference_buffer();

        let props = AlignmentProperties::estimate(
            &mut bam,
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
        let mut bam =
            bam::Reader::from_path("tests/resources/tumor-first30000.reads_with_soft_clips.bam")
                .unwrap();
        let mut reference_buffer = reference_buffer();

        let props = AlignmentProperties::estimate(
            &mut bam,
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
        let mut bam = bam::Reader::from_path(
            "tests/resources/tumor-first30000.bunch_of_reads_made_single_ended.bam",
        )
        .unwrap();
        let mut reference_buffer = reference_buffer();

        let props = AlignmentProperties::estimate(
            &mut bam,
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
