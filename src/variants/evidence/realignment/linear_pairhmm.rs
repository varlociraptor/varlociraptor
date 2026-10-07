// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Forward algorithm of the pair HMM computed in linear probability space.
//!
//! This mirrors `bio::stats::pairhmm::PairHMM::prob_related` for the semiglobal case (free start
//! and end gaps on the reference side), including the banding by minimal edit distance, but keeps
//! the three state matrices as plain `f64` probabilities. This replaces the three log-space
//! additions per cell (each an `exp` and a `log1p`) by multiplications and additions.
//!
//! Gap extensions emit their base and a gap state is left with the complement of its own
//! extension probability, and cells outside the band are cleared; see rust-bio/rust-bio#700 for
//! the corresponding deviations of the log-space implementation.
//!
//! Only the cells inside the band are visited: a cell can be inside the band only if its top
//! left or left neighbour (previous column) or its top neighbour (this column) is, so each column
//! is scanned over the ranges written in the previous column, extended by one, plus the chain of
//! top neighbours, and cells outside of those ranges are known to be empty. An alignment is only
//! allowed to start at columns from which the read can still be consumed within the band.
//!
//! Linear space is safe because of the band: every computed cell is reached by a path with at
//! most `max_edit_dist` edits, so its value is at least the product of that path, i.e. no smaller
//! than `min_edit_emission^max_edit_dist * min_match_emission^len_y`. With the smallest gap
//! probability around 1e-6 and base qualities down to 2 this stays far above the smallest normal
//! `f64` for bands of up to `MAX_EDIT_DIST` edits and read windows of up to `MAX_LEN_Y` bases.
//! Wider bands, longer windows and unbanded computations have to use the log-space
//! implementation.

use std::cmp;

use bio::stats::pairhmm::GapParameters;
use bio::stats::LogProb;

use super::pairhmm::{GapParams, ReadEmission};

/// Widest band for which the linear-space forward algorithm cannot underflow (together with
/// `MAX_LEN_Y`).
pub(crate) const MAX_EDIT_DIST: usize = 40;

/// Longest read window for which the linear-space forward algorithm cannot underflow (together
/// with `MAX_EDIT_DIST`).
pub(crate) const MAX_LEN_Y: usize = 128;

#[derive(Debug, Clone)]
pub(crate) struct LinearPairHMM {
    // state M (match/mismatch), X (reference base emitted alone) and Y (read base emitted alone),
    // as two alternating columns indexed by read position + 1
    m: [Vec<f64>; 2],
    x: [Vec<f64>; 2],
    y: [Vec<f64>; 2],
    min_edit_dist: [Vec<usize>; 2],
    // inclusive ranges of the cells (read position + 1) written in each of the two columns; every
    // other cell of a column holds probability 0 and an infinite edit distance
    live: [Vec<(usize, usize)>; 2],
    // scratch space for the candidate ranges of a column
    candidates: Vec<(usize, usize)>,
    // probability of the alignments ending in each column (free end gap in x)
    cols: Vec<f64>,
    no_gap: f64,
    no_gap_x_extend: f64,
    no_gap_y_extend: f64,
    gap_x: f64,
    gap_y: f64,
    gap_x_extend: f64,
    gap_y_extend: f64,
    do_gap_x_extend: bool,
    do_gap_y_extend: bool,
}

impl LinearPairHMM {
    pub(crate) fn new(gap_params: &GapParams) -> Self {
        let lin = |p: LogProb| p.exp();
        LinearPairHMM {
            m: [Vec::new(), Vec::new()],
            x: [Vec::new(), Vec::new()],
            y: [Vec::new(), Vec::new()],
            min_edit_dist: [Vec::new(), Vec::new()],
            live: [Vec::new(), Vec::new()],
            candidates: Vec::new(),
            cols: Vec::new(),
            no_gap: lin(gap_params
                .prob_gap_x()
                .ln_add_exp(gap_params.prob_gap_y())
                .ln_one_minus_exp()),
            no_gap_x_extend: lin(gap_params.prob_gap_x_extend().ln_one_minus_exp()),
            no_gap_y_extend: lin(gap_params.prob_gap_y_extend().ln_one_minus_exp()),
            gap_x: lin(gap_params.prob_gap_x()),
            gap_y: lin(gap_params.prob_gap_y()),
            gap_x_extend: lin(gap_params.prob_gap_x_extend()),
            gap_y_extend: lin(gap_params.prob_gap_y_extend()),
            do_gap_x_extend: gap_params.prob_gap_x_extend() != LogProb::ln_zero(),
            do_gap_y_extend: gap_params.prob_gap_y_extend() != LogProb::ln_zero(),
        }
    }

    /// Probability that the read window (y) is related to the allele window `allele` (x, upper
    /// case bases) via any semiglobal alignment, banded to alignments with at most `max_edit_dist`
    /// edits (which must not exceed `MAX_EDIT_DIST`, see the module documentation).
    pub(crate) fn prob_related(
        &mut self,
        read: &ReadEmission,
        allele: &[u8],
        max_edit_dist: usize,
    ) -> LogProb {
        assert!(
            max_edit_dist <= MAX_EDIT_DIST,
            "bug: band of {} edits is too wide for the linear-space pair HMM",
            max_edit_dist
        );
        let len_x = allele.len();
        let len_y = read.len();
        assert!(
            len_y <= MAX_LEN_Y,
            "bug: read window of {} bases is too long for the linear-space pair HMM",
            len_y
        );

        for k in 0..2 {
            self.m[k].clear();
            self.x[k].clear();
            self.y[k].clear();
            self.min_edit_dist[k].clear();
            self.m[k].resize(len_y + 1, 0.0);
            self.x[k].resize(len_y + 1, 0.0);
            self.y[k].resize(len_y + 1, 0.0);
            self.min_edit_dist[k].resize(len_y + 1, usize::MAX);
            self.live[k].clear();
        }
        self.cols.clear();
        self.cols.reserve(len_x);

        let mut prev = 0;
        let mut curr = 1;
        self.m[prev][0] = 1.0;

        // METHOD: the alignment may start at any reference offset (free start gap in x), but an
        // alignment starting at column c has to consume the whole read within the remaining
        // len_x - c columns, which needs at least len_y - (len_x - c) insertions. Beyond
        // last_start it cannot stay within the band, so its cells would only be computed to be
        // pruned: skipping those starts leaves the result unchanged.
        let last_start = (len_x + max_edit_dist).saturating_sub(len_y);

        for (i, &ref_base) in allele.iter().enumerate() {
            // Clear what this buffer holds from two columns ago, so that every cell outside of
            // the ranges written below is zero (see `live`).
            for &(a, b) in &self.live[curr] {
                for v in &mut self.m[curr][a..=b] {
                    *v = 0.0;
                }
                for v in &mut self.x[curr][a..=b] {
                    *v = 0.0;
                }
                for v in &mut self.y[curr][a..=b] {
                    *v = 0.0;
                }
                for v in &mut self.min_edit_dist[curr][a..=b] {
                    *v = usize::MAX;
                }
            }
            self.live[curr].clear();
            // cell 0 of this column: no alignment starts here (it gets its start probability
            // once it becomes the previous column)
            self.m[curr][0] = 0.0;

            let fresh_start = i <= last_start;
            if fresh_start {
                self.m[prev][0] += 1.0;
                self.min_edit_dist[prev][0] = 0;
            } else {
                self.m[prev][0] = 0.0;
                self.min_edit_dist[prev][0] = usize::MAX;
            }
            self.min_edit_dist[curr][0] = usize::MAX;

            // METHOD: a cell can only be inside the band if its top left or left neighbour
            // (previous column) or its top neighbour (this column) is. So the candidates are
            // the ranges of the previous column extended by one to the right, cell 1 if the
            // alignment may start here, and whatever the chain of top neighbours reaches.
            let mut candidates = std::mem::take(&mut self.candidates);
            candidates.clear();
            if fresh_start {
                candidates.push((1, 1));
            }
            candidates.extend(
                self.live[prev]
                    .iter()
                    .map(|&(a, b)| (a, (b + 1).min(len_y))),
            );
            let mut next = 1;
            for &(a, b) in &candidates {
                let mut j_ = a.max(next);
                if j_ > b {
                    continue;
                }
                loop {
                    let j = j_ - 1;
                    let med_topleft = self.min_edit_dist[prev][j];
                    let med_top = self.min_edit_dist[curr][j];
                    let med_left = self.min_edit_dist[prev][j_];
                    let live = cmp::min(med_topleft, cmp::min(med_top, med_left)) <= max_edit_dist;

                    if live {
                        let (emit_xy, is_match) = read.prob_match_mismatch_linear(j, ref_base);

                        // match or mismatch, coming from M, X (extended with gap_y_extend) or Y
                        // (extended with gap_x_extend)
                        let m = emit_xy
                            * (self.no_gap * self.m[prev][j]
                                + self.no_gap_y_extend * self.x[prev][j]
                                + self.no_gap_x_extend * self.y[prev][j]);

                        // reference base emitted alone (gap in y); emission probability of x is 1
                        let mut x = self.gap_y * self.m[prev][j_];
                        if self.do_gap_y_extend {
                            x += self.gap_y_extend * self.x[prev][j_];
                        }

                        // read base emitted alone (gap in x)
                        let mut y = self.gap_x * self.m[curr][j];
                        if self.do_gap_x_extend {
                            y += self.gap_x_extend * self.y[curr][j];
                        }
                        y *= read.prob_insertion_linear(j);

                        self.m[curr][j_] = m;
                        self.x[curr][j_] = x;
                        self.y[curr][j_] = y;
                        self.min_edit_dist[curr][j_] = cmp::min(
                            if is_match {
                                med_topleft
                            } else {
                                med_topleft.saturating_add(1)
                            },
                            cmp::min(med_left.saturating_add(1), med_top.saturating_add(1)),
                        );
                        match self.live[curr].last_mut() {
                            Some((_, last)) if *last + 1 == j_ => *last = j_,
                            _ => self.live[curr].push((j_, j_)),
                        }
                    }
                    // outside of the band the cell already holds zero

                    j_ += 1;
                    if j_ > len_y || (j_ > b && !live) {
                        break;
                    }
                }
                next = j_;
            }
            self.candidates = candidates;

            // the alignment may end at any reference offset (free end gap in x)
            self.cols
                .push(self.m[curr][len_y] + self.x[curr][len_y] + self.y[curr][len_y]);

            std::mem::swap(&mut curr, &mut prev);
        }

        let p = LogProb(self.cols.iter().sum::<f64>().ln());
        assert!(!p.is_nan());
        // the sum over all paths can exceed 1.0, especially in repeats
        if p > LogProb::ln_one() {
            LogProb::ln_one()
        } else {
            p
        }
    }
}

#[cfg(test)]
mod tests {
    use bio::stats::{LogProb, Prob};
    use rand::prelude::*;
    use rust_htslib::bam;

    use super::*;

    fn record(seq: &[u8], qual: &[u8]) -> bam::Record {
        let mut rec = bam::Record::new();
        let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(seq.len() as u32)]);
        rec.set(b"read", Some(&cigar), seq, qual);
        rec
    }

    /// Log-sum-exp with the exact `exp` (`LogProb::ln_sum_exp` uses an approximation).
    fn lse(probs: &[LogProb]) -> LogProb {
        let pmax = probs
            .iter()
            .cloned()
            .fold(LogProb::ln_zero(), |a, b| if b > a { b } else { a });
        if pmax == LogProb::ln_zero() {
            return pmax;
        }
        LogProb(*pmax + probs.iter().map(|p| (**p - *pmax).exp()).sum::<f64>().ln())
    }

    /// Exact log-space forward (no approximations) with the same recurrences and band, as the
    /// ground truth.
    fn exact(
        gap_params: &GapParams,
        read: &ReadEmission,
        allele: &[u8],
        max_edit_dist: usize,
    ) -> LogProb {
        let ln = |p: f64| LogProb(p.ln());
        let hmm = LinearPairHMM::new(gap_params);
        let len_y = read.len();
        let zero = LogProb::ln_zero();
        let mut m = vec![vec![zero; len_y + 1]; 2];
        let mut x = m.clone();
        let mut y = m.clone();
        let mut med = vec![vec![usize::MAX; len_y + 1]; 2];
        let (mut prev, mut curr) = (0, 1);
        m[prev][0] = LogProb::ln_one();
        let mut cols = Vec::new();
        for &ref_base in allele {
            m[prev][0] = lse(&[m[prev][0], LogProb::ln_one()]);
            med[prev][0] = 0;
            med[curr][0] = 0;
            for j in 0..len_y {
                let j_ = j + 1;
                let (tl, top, left) = (med[prev][j], med[curr][j], med[prev][j_]);
                if tl.min(top).min(left) > max_edit_dist {
                    m[curr][j_] = zero;
                    x[curr][j_] = zero;
                    y[curr][j_] = zero;
                    med[curr][j_] = usize::MAX;
                    continue;
                }
                let (e, is_match) = read.prob_match_mismatch_linear(j, ref_base);
                m[curr][j_] = ln(e)
                    + lse(&[
                        ln(hmm.no_gap) + m[prev][j],
                        ln(hmm.no_gap_y_extend) + x[prev][j],
                        ln(hmm.no_gap_x_extend) + y[prev][j],
                    ]);
                x[curr][j_] = lse(&[
                    ln(hmm.gap_y) + m[prev][j_],
                    ln(hmm.gap_y_extend) + x[prev][j_],
                ]);
                y[curr][j_] = ln(read.prob_insertion_linear(j))
                    + lse(&[
                        ln(hmm.gap_x) + m[curr][j],
                        ln(hmm.gap_x_extend) + y[curr][j],
                    ]);
                med[curr][j_] = (if is_match { tl } else { tl.saturating_add(1) })
                    .min(left.saturating_add(1))
                    .min(top.saturating_add(1));
            }
            cols.push(lse(&[m[curr][len_y], x[curr][len_y], y[curr][len_y]]));
            std::mem::swap(&mut prev, &mut curr);
            for v in m[curr]
                .iter_mut()
                .chain(x[curr].iter_mut())
                .chain(y[curr].iter_mut())
            {
                *v = zero;
            }
            for v in med[curr].iter_mut() {
                *v = usize::MAX;
            }
        }
        let p = lse(&cols);
        if p > LogProb::ln_one() {
            LogProb::ln_one()
        } else {
            p
        }
    }

    fn mutate(rng: &mut StdRng, seq: &[u8], edits: usize) -> Vec<u8> {
        let mut seq = seq.to_vec();
        for _ in 0..edits {
            let pos = rng.gen_range(0..seq.len());
            match rng.gen_range(0..3) {
                0 => seq[pos] = *b"ACGT".choose(rng).unwrap(),
                1 => seq.insert(pos, *b"ACGT".choose(rng).unwrap()),
                _ => {
                    seq.remove(pos);
                }
            }
        }
        seq
    }

    fn assert_close(linear: LogProb, reference: LogProb, tol: f64, what: &str) {
        assert!(
            (*linear - *reference).abs() <= tol,
            "{}: linear {} vs exact {} (diff {})",
            what,
            *linear,
            *reference,
            (*linear - *reference).abs()
        );
    }

    /// Random windows with up to `max_edits` edits, banded to `edits + band`.
    fn check_random(seed: u64, gap_params: &GapParams, max_edits: usize, band: usize) {
        let mut rng = StdRng::seed_from_u64(seed);
        let mut hmm = LinearPairHMM::new(gap_params);
        for case in 0..200 {
            let len = rng.gen_range(20..160);
            let allele: Vec<u8> = (0..len)
                .map(|_| *b"ACGT".choose(&mut rng).unwrap())
                .collect();
            let start = rng.gen_range(0..len / 2);
            let end = rng.gen_range(start + 10..=len);
            let edits = rng.gen_range(0..=max_edits);
            let window = allele[start..end].to_vec();
            let read = mutate(&mut rng, &window, edits);
            let qual: Vec<u8> = (0..read.len()).map(|_| rng.gen_range(2..=40)).collect();
            let rec = record(&read, &qual);
            let emission = ReadEmission::new(rec.seq(), rec.qual(), None, None);
            let max_edit_dist = edits + band;

            let truth = exact(gap_params, &emission, &allele, max_edit_dist);
            let linear = hmm.prob_related(&emission, &allele, max_edit_dist);
            assert_close(
                linear,
                truth,
                1e-9 * truth.abs().max(1.0),
                &format!("case {}", case),
            );
        }
    }

    #[test]
    fn matches_exact_forward_on_random_windows() {
        check_random(7, &GapParams::default(), 6, 4);
    }

    #[test]
    fn matches_exact_forward_with_gap_extension() {
        let gap_params = GapParams {
            prob_insertion_artifact: LogProb::from(Prob(1e-3)),
            prob_deletion_artifact: LogProb::from(Prob(2e-3)),
            prob_insertion_extend_artifact: LogProb::from(Prob(0.3)),
            prob_deletion_extend_artifact: LogProb::from(Prob(0.4)),
        };
        check_random(11, &gap_params, 8, 4);
    }

    #[test]
    fn widest_band_stays_finite() {
        // The worst case the band bound allows: base quality 2 everywhere, as many mismatches
        // as the band permits, on the longest read window.
        let gap_params = GapParams::default();
        let allele: Vec<u8> = (0..200).map(|i| b"ACGT"[i % 4]).collect();
        let mut read = allele[30..158].to_vec();
        for k in 0..MAX_EDIT_DIST {
            read[k * 3] = b'A' + ((read[k * 3] - b'A' + 1) % 4);
        }
        let qual = vec![2u8; read.len()];
        let rec = record(&read, &qual);
        let emission = ReadEmission::new(rec.seq(), rec.qual(), None, None);
        let mut hmm = LinearPairHMM::new(&gap_params);
        let linear = hmm.prob_related(&emission, &allele, MAX_EDIT_DIST);
        assert!(linear.is_finite(), "{}", *linear);
        let truth = exact(&gap_params, &emission, &allele, MAX_EDIT_DIST);
        assert_close(linear, truth, 1e-9 * truth.abs(), "widest band");
    }

    #[test]
    fn perfect_match_is_close_to_one() {
        let gap_params = GapParams::default();
        let allele: Vec<u8> = b"ACGTTGCAAGGCTTACGGATCCATGCATGCAAGTC".to_vec();
        let qual = vec![40u8; allele.len()];
        let rec = record(&allele, &qual);
        let emission = ReadEmission::new(rec.seq(), rec.qual(), None, None);
        let mut hmm = LinearPairHMM::new(&gap_params);
        let p = hmm.prob_related(&emission, &allele, 4);
        // product of the no-miscall probabilities, up to the alternative paths
        assert!(*p > (allele.len() as f64) * (1.0f64 - 1e-4).ln() - 1e-3 && *p <= 0.0);
    }
}
