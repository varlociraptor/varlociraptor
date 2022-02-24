// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::f64;
use std::ops::Deref;
use std::str;
use std::u32;

use anyhow::Result;
use bio::stats::{LogProb, Prob};
use counter::Counter;
use itertools::Itertools;
use ordered_float::NotNan;
use rust_htslib::bam::{self, record::Cigar};
use statrs::statistics::{Data, Distribution, OrderStatistics};

use crate::reference;
use crate::utils::homopolymers::is_homopolymer_seq;
use crate::utils::homopolymers::{extend_homopolymer_stretch, is_homopolymer_iter};
use crate::utils::SimpleCounter;

pub(crate) const MIN_HOMOPOLYMER_LEN: usize = 4;

use crate::variants::evidence::realignment::pairhmm::{GapParams, HopParams};
use rayon::prelude::*;

const NUM_FRAGMENTS: usize = 1000000;

fn default_homopolymer_error_model() -> HashMap<i16, f64> {
    let mut model = HashMap::new();
    model.insert(0, 0.9975414130829068);
    model.insert(1, 0.0010076175889726332);
    model.insert(-1, 0.0010076175889726332);
    model.insert(-2, 0.00020152351779452663);
    model.insert(2, 0.00010076175889726332);
    model.insert(3, 5.038087944863166e-05);
    model.insert(-3, 9.068558300753699e-05);

    model
}

fn default_max_mapq() -> u8 {
    60
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub(crate) struct AlignmentProperties {
    pub(crate) insert_size: Option<InsertSize>,
    pub(crate) max_del_cigar_len: Option<u32>,
    pub(crate) max_ins_cigar_len: Option<u32>,
    pub(crate) frac_max_softclip: Option<f64>,
    pub(crate) max_read_len: u32,
    #[serde(default = "default_max_mapq")]
    pub(crate) max_mapq: u8,
    pub(crate) homopolymer_counts: HashMap<(u8, i16), usize>,
    pub(crate) gap_counts: HashMap<i16, usize>,
    #[serde(default = "default_homopolymer_error_model")]
    pub(crate) wildtype_homopolymer_error_model: HashMap<i16, f64>,
    #[serde(default)]
    initial: bool,
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

fn gap_and_hop_counts(
    record: &bam::Record,
    refseq: &[u8],
) -> (SimpleCounter<i16>, SimpleCounter<(u8, i16)>) {
    let mut hop_counts = SimpleCounter::default();
    let mut gap_counts = SimpleCounter::default();
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
                let l = l as usize;
                if l < i16::MAX as usize {
                    let base = refseq[rpos];
                    if is_homopolymer_seq(&refseq[rpos..rpos + l]) {
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
                            hop_counts.incr((base, -(l as i16)));
                        }
                    } else {
                        gap_counts.incr(-(l as i16));
                    }
                }
                rpos += l as usize;
            }
            Cigar::Ins(l) => {
                let l = l as usize;
                if l < i16::MAX as usize {
                    let base = qseq[qpos];
                    if is_homopolymer_iter((qpos..qpos + l).map(|i| qseq[i])) {
                        let mut len =
                            l + extend_homopolymer_stretch(qseq[qpos], &mut refseq[rpos..].iter());
                        if rpos > 0 {
                            len += extend_homopolymer_stretch(
                                qseq[qpos],
                                &mut refseq[..rpos].iter().rev(),
                            );
                        }
                        if len >= MIN_HOMOPOLYMER_LEN {
                            hop_counts.incr((base, l as i16));
                        }
                    } else {
                        gap_counts.incr(l as i16);
                    }
                }
                qpos += l as usize;
            }
            Cigar::Match(l) | Cigar::Diff(l) | Cigar::Equal(l) => {
                let l = l as usize;
                for (base, stretch) in &refseq[rpos..rpos + l].iter().group_by(|c| **c) {
                    if stretch.count() >= MIN_HOMOPOLYMER_LEN {
                        hop_counts.incr((base, 0));
                    }
                }
                qpos += l as usize;
                rpos += l as usize;
            }
            Cigar::SoftClip(l) => {
                qpos += l as usize;
            }
            Cigar::RefSkip(l) => {
                rpos += l as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => continue,
        }
    }
    (gap_counts, hop_counts)
}

struct CigarStats {
    is_regular: bool,
    has_soft_clip: bool,
    frac_max_softclip: Option<f64>,
    max_del: Option<u32>,
    max_ins: Option<u32>,
}

fn cigar_stats(record: &bam::Record, allow_hardclips: bool) -> CigarStats {
    let norm = |j| NotNan::new(j as f64 / record.seq().len() as f64).unwrap();

    let mut is_regular = true;
    let mut has_soft_clip = false;
    let mut frac_max_softclip = None;
    let mut max_del_cigar_len = None;
    let mut max_ins_cigar_len = None;
    let iter = if let Some(cigar) = record.cigar_cached() {
        Box::new(cigar.iter().copied()) as Box<dyn Iterator<Item = Cigar>>
    } else {
        Box::new(iter_cigar(record))
    };
    for c in iter {
        match c {
            Cigar::SoftClip(l) => {
                let s = norm(l);
                if let Some(ref mut maxclipfrac) = frac_max_softclip {
                    *maxclipfrac = *cmp::max(s, NotNan::new(*maxclipfrac).unwrap())
                } else {
                    frac_max_softclip = Some(*s);
                }
                is_regular = false;
                has_soft_clip = true;
            }
            Cigar::Del(l) => {
                if let Some(ref mut maxlen) = max_del_cigar_len {
                    *maxlen = cmp::max(l, *maxlen);
                } else {
                    max_del_cigar_len = Some(l);
                }
                is_regular = false;
            }
            Cigar::Ins(l) => {
                if let Some(ref mut maxlen) = max_ins_cigar_len {
                    *maxlen = cmp::max(l, *maxlen);
                } else {
                    max_ins_cigar_len = Some(l);
                }
                is_regular = false;
            }
            Cigar::HardClip(_) if !allow_hardclips => {
                is_regular = false;
            }
            _ => continue,
        }
    }

    CigarStats {
        is_regular,
        has_soft_clip,
        frac_max_softclip,
        max_del: max_del_cigar_len,
        max_ins: max_ins_cigar_len,
    }
}

trait OptionMax {
    fn max(self, other: Option<f64>) -> Option<f64>;
}

impl OptionMax for Option<f64> {
    fn max(self, other: Option<f64>) -> Option<f64> {
        match (self, other) {
            (None, None) => None,
            (Some(a), None) => Some(a),
            (None, Some(b)) => Some(b),
            (Some(a), Some(b)) => Some(a.max(b)),
        }
    }
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
            homopolymer_counts: Default::default(),
            gap_counts: Default::default(),
            wildtype_homopolymer_error_model: HashMap::new(),
            initial: true,
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
            cigar_stats: CigarStats,
            homopolymer_counts: SimpleCounter<(u8, i16)>,
            gap_counts: SimpleCounter<i16>,
            insert_size: Option<f64>,
        }

        #[derive(Default)]
        struct AlignmentStats {
            n_reads: usize,
            max_mapq: u8,
            max_read_len: u32,
            n_softclips: u32,
            n_not_usable: u32,
            frac_max_softclip: Option<f64>,
            max_del: Option<u32>,
            max_ins: Option<u32>,
            homopolymer_counts: SimpleCounter<(u8, i16)>,
            gap_counts: SimpleCounter<i16>,
            tlens: Vec<f64>,
        }

        impl AlignmentStats {
            fn update(&mut self, other: Self) {
                self.n_reads += other.n_reads;
                self.max_mapq = self.max_mapq.max(other.max_mapq);
                self.max_read_len = self.max_read_len.max(other.max_read_len);
                self.n_softclips += other.n_softclips;
                self.n_not_usable += other.n_not_usable;
                self.homopolymer_counts += other.homopolymer_counts;
                self.gap_counts += other.gap_counts;
                self.frac_max_softclip = self.frac_max_softclip.max(other.frac_max_softclip);
                self.max_ins = self.max_ins.max(other.max_ins);
                self.max_del = self.max_del.max(other.max_del);
            }
        }

        let header = bam.header().clone();
        let tid_to_tname: HashMap<_, _> = header
            .target_names()
            .iter()
            .map(|tname| {
                let tid = header.tid(tname).unwrap();
                (tid, *tname)
            })
            .collect();
        let reference_seqs: HashMap<&[u8], Vec<u8>> = tid_to_tname
            .values()
            .map(|tname| {
                (
                    *tname,
                    reference_buffer
                        .seq(std::str::from_utf8(tname).unwrap())
                        .unwrap()
                        .deref()
                        .clone(),
                )
            })
            .collect();

        let buf_size = rayon::current_num_threads();
        let mut n_records_read = 0;
        let mut n_records_skipped = 0;
        let mut all_stats = AlignmentStats::default();
        let mut record_flag_stats = RecordFlagStats::new();
        let mut eof = false;
        let mut records = bam.records();
        while !eof && n_records_read < NUM_FRAGMENTS {
            let mut record_buffer = Vec::with_capacity(buf_size);
            while record_buffer.len() < buf_size && !eof {
                match records.next() {
                    None => {
                        eof = true;
                        break;
                    }
                    Some(res) => {
                        let mut record = res?;
                        record_flag_stats.update(&record);
                        n_records_read += 1;
                        // don't cache cigar for ultralong reads
                        if record.seq_len() <= 32000 {
                            record.cache_cigar();
                        }
                        record_buffer.push(record);
                    }
                }
            }
            let batch_stats = record_buffer
                .into_par_iter()
                .enumerate()
                .filter(|(_i, record)| {
                    !(record.mapq() == 0
                        || record.is_duplicate()
                        || record.is_quality_check_failed()
                        || record.is_unmapped())
                })
                .map(|(_i, record)| {
                    let cigar_stats = cigar_stats(&record, allow_hardclips);

                    let (gap_counts, homopolymer_counts) = gap_and_hop_counts(
                        &record,
                        &reference_seqs[tid_to_tname[&(record.tid() as u32)]],
                    );

                    let insert_size = if !cigar_stats.is_regular {
                        None
                    } else {
                        // record insert size
                        Some(record.insert_size().abs() as f64)
                    };
                    RecordStats {
                        mapq: record.mapq(),
                        read_len: record.seq().len() as u32,
                        cigar_stats,
                        homopolymer_counts,
                        gap_counts,
                        insert_size,
                    }
                })
                .fold(AlignmentStats::default, |mut acc, rs| {
                    acc.n_reads += 1;
                    acc.max_mapq = acc.max_mapq.max(rs.mapq);
                    acc.max_read_len = acc.max_read_len.max(rs.read_len);
                    acc.n_softclips += rs.cigar_stats.has_soft_clip as u32;
                    acc.n_not_usable += (!rs.cigar_stats.is_regular) as u32;
                    acc.homopolymer_counts += rs.homopolymer_counts;
                    acc.gap_counts += rs.gap_counts;
                    if let Some(insert_size) = rs.insert_size {
                        acc.tlens.push(insert_size);
                    }
                    acc.frac_max_softclip =
                        acc.frac_max_softclip.max(rs.cigar_stats.frac_max_softclip);
                    acc.max_ins = acc.max_ins.max(rs.cigar_stats.max_ins);
                    acc.max_del = acc.max_del.max(rs.cigar_stats.max_del);
                    acc
                })
                .reduce(AlignmentStats::default, |mut a, b| {
                    a.update(b);
                    a
                });
            n_records_skipped += buf_size - batch_stats.n_reads;
            all_stats.update(batch_stats);
        }

        properties.homopolymer_counts = all_stats.homopolymer_counts.clone();
        properties.gap_counts = all_stats.gap_counts.clone();

        properties.wildtype_homopolymer_error_model = {
            let n = all_stats
                .homopolymer_counts
                .values()
                .filter(|count| **count >= 10)
                .sum::<usize>() as f64;
            // group by length, i.e. discard base information
            let grouped = all_stats
                .homopolymer_counts
                .iter()
                .sorted_unstable_by_key(|((_, length), _)| length)
                .group_by(|((_, length), _)| length);
            grouped
                .into_iter()
                .map(|(length, group)| (length, group.map(|(_, count)| count).sum::<usize>()))
                .map(|(length, count)| (*length, count as f64 / n))
                .collect()
        };

        properties.max_read_len = all_stats.max_read_len;
        properties.max_del_cigar_len = all_stats.max_del;
        properties.max_ins_cigar_len = all_stats.max_ins;
        properties.frac_max_softclip = all_stats.frac_max_softclip;
        properties.max_mapq = all_stats.max_mapq;
        properties.max_read_len = all_stats.max_read_len;

        if properties.max_del_cigar_len.is_none() {
            warn!(
                "No deletion CIGAR operations found in first 10000 alignments. \
                Varlociraptor will be unable to estimate the sampling bias for deletions."
            );
        }
        if properties.max_ins_cigar_len.is_none() {
            warn!(
                "No deletion CIGAR operations found in first 10000 alignments. \
                Varlociraptor will be unable to estimate the sampling bias for insertions."
            );
        }
        if properties.frac_max_softclip.is_none() {
            warn!(
                "No softclip CIGAR operations found in the first 10000 alignments. \
                Varlociraptor will be unable to estimate the sampling bias for larger indels."
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

fn exponential_mle<V: Into<usize>>(value_counts: impl Iterator<Item = (V, usize)>) -> f64 {
    let (sum, count) = value_counts.fold((0usize, 0usize), |(sum, counts), (value, count)| {
        (sum + value.into() * count, counts + count)
    });
    // the MLE of the exponential distribution's lambda parameter is simply 1 / sample_mean
    let lambda = count as f64 / sum as f64;
    // … the estimator is slightly biased, which can be corrected for:
    lambda - (lambda / count as f64)
}

impl AlignmentProperties {
    fn gap_params(&self) -> GapParams {
        let counts = |length_predicate: fn(i16) -> bool| {
            self.gap_counts
                .iter()
                .filter(|&(length, _)| length_predicate(*length))
                .map(|(length, count)| (length.abs() as usize, *count))
                .collect_vec()
        };
        let insertion_gap_counts = counts(|l| l > 0);
        let deletion_gap_counts = counts(|l| l < 0);
        let insertion_lambda = exponential_mle(insertion_gap_counts.into_iter());
        let insertion_prob = (-insertion_lambda).exp();
        let deletion_lambda = exponential_mle(deletion_gap_counts.into_iter());
        let deletion_prob = (-deletion_lambda).exp();
        dbg!(
            insertion_lambda,
            insertion_prob,
            deletion_lambda,
            deletion_prob
        );
        todo!()
    }

    fn hop_params(&self) -> HopParams {
        let counts = |base, length_predicate: fn(i16) -> bool| {
            self.homopolymer_counts
                .iter()
                .filter(|((char, length), _)| *char == base && length_predicate(*length))
                .map(|((_, length), count)| (length.abs() as usize, *count))
                .collect_vec()
        };
        let probs = [b'A', b'C', b'G', b'T']
            .iter()
            .map(|base| {
                let insertion_counts = counts(*base, |l| l > 0);
                let deletion_counts = counts(*base, |l| l < 0);
                let insertion_lambda = exponential_mle(insertion_counts.into_iter());
                let insertion_prob = (-insertion_lambda).exp();
                let deletion_lambda = exponential_mle(deletion_counts.into_iter());
                let deletion_prob = (-deletion_lambda).exp();
                dbg!(
                    insertion_lambda,
                    insertion_prob,
                    deletion_lambda,
                    deletion_prob
                );
            })
            .collect_vec();
        todo!()
    }
}

/// Expected insert size in terms of mean and standard deviation.
/// This should be estimated from unsorted(!) bam files to avoid positional biases.
#[derive(Copy, Clone, Debug, Default, Serialize, Deserialize)]
pub(crate) struct InsertSize {
    pub(crate) mean: f64,
    pub(crate) sd: f64,
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

        let props = AlignmentProperties::estimate(&mut bam, false, &mut reference_buffer).unwrap();
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

        let props = AlignmentProperties::estimate(&mut bam, false, &mut reference_buffer).unwrap();
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

        let props = AlignmentProperties::estimate(&mut bam, false, &mut reference_buffer).unwrap();
        println!("{:?}", props);

        assert!(props.insert_size.is_none());
        assert_eq!(props.max_del_cigar_len, None);
        assert_eq!(props.max_ins_cigar_len, None);
        assert_eq!(props.frac_max_softclip, Some(0.03));
    }
}
