//! This module provides many Levenshtein distance routines.
//!
//! These distance functions share the same efficient underlying SIMD-accelerated implementation:
//! * `levenshtein_exp` for low number of edits, otherwise `levenshtein`
//! * `rdamerau_exp` for low number of edits, otherwise `rdamerau`
//! * `levenshtein_simd_k`
//! * `levenshtein_simd_k_with_opts`
//!
//! These search functions share the same efficient underlying SIMD-accelerated implementation:
//! * `levenshtein_search`
//! * `levenshtein_search_simd`
//! * `levenshtein_search_simd_with_opts`

use std::*;
use super::*;
use super::jewel::*;

/// A struct holding the edit costs for mismatches, gaps, and possibly transpositions.
///
/// This should be used as a parameter for Levenshtein distance or search routines.
#[derive(Copy, Clone, Debug)]
pub struct EditCosts {
    mismatch_cost: u8,
    gap_cost: u8,
    start_gap_cost: u8,
    transpose_cost: Option<u8>
}

impl EditCosts {
    /// Create a new `EditCosts` struct, checking for whether the specified costs are valid.
    ///
    /// # Arguments
    /// * `mismatch_cost` - cost of a mismatch edit, which must be positive
    /// * `gap_cost` - cost of a gap, which must be positive
    /// * `start_gap_cost` - additional cost of starting a gap, for affine gap costs; this can
    /// be zero for linear gap costs
    /// * `transpose_cost` - cost of a transpose, which must be cheaper than doing the equivalent
    /// operation with mismatches and gaps
    pub fn new(mismatch_cost: u8, gap_cost: u8, start_gap_cost: u8, transpose_cost: Option<u8>) -> Self {
        assert!(mismatch_cost > 0);
        assert!(gap_cost > 0);

        if let Some(cost) = transpose_cost {
            assert!(cost > 0);
            // transpose cost must be cheaper than doing the equivalent with other edits
            assert!((cost >> 1) < mismatch_cost);
            assert!((cost >> 1) < gap_cost);
        }

        Self{
            mismatch_cost: mismatch_cost,
            gap_cost: gap_cost,
            start_gap_cost: start_gap_cost,
            transpose_cost: transpose_cost
        }
    }

    /// For Levenshtein searches, the cost of transpositions must be less than or equal to cost of
    /// gaps.
    ///
    /// This is important for free gaps at the beginning of the needle to be unable to take priority
    /// over transpositions, as it is possible to emulate a transposition with two gaps.
    fn check_search(&self) {
        if let Some(cost) = self.transpose_cost {
            assert!(cost <= self.start_gap_cost + self.gap_cost);
        }
    }
}

/// Costs for Levenshtein distance, where mismatches and gaps both have a cost of 1, and
/// transpositions are not allowed.
pub const LEVENSHTEIN_COSTS: EditCosts = EditCosts{mismatch_cost: 1, gap_cost: 1, start_gap_cost: 0, transpose_cost: None};
/// Costs for restricted Damerau-Levenshtein distance, where mismatches, gaps, and transpositions
/// all have a cost of 1.
pub const RDAMERAU_COSTS: EditCosts = EditCosts{mismatch_cost: 1, gap_cost: 1, start_gap_cost: 0, transpose_cost: Some(1)};

/// Returns the Levenshtein distance between two strings using the naive scalar algorithm.
///
/// # Arguments
/// * `a` - first string (slice)
/// * `b` - second string (slice)
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::levenshtein::*;
/// let dist = levenshtein_naive(b"abc", b"ab");
///
/// assert!(dist == 1);
/// ```
pub fn levenshtein_naive<T: PartialEq>(a: &[T], b: &[T]) -> u32 {
    levenshtein_naive_with_opts(a, b, false, LEVENSHTEIN_COSTS).0
}

/// Returns the Levenshtein distance between two strings using the naive scalar algorithm.
///
/// # Arguments
/// * `a` - first string (&str)
/// * `b` - second string (&str)
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::levenshtein::*;
/// let dist = levenstein_naive_str("abc", "ab");
///
/// assert!(dist == 1);
/// ```
pub fn levenstein_naive_str(a: &str, b: &str) -> u32 {
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    levenshtein_naive(&a, &b)
}

/// Returns the Levenshtein distance between two strings and optionally, the edit traceback,
/// using the naive scalar algorithm, with extra options.
///
/// # Arguments
/// * `a` - first string (slice)
/// * `b` - second string (slice)
/// * `trace_on` - whether to return the traceback, the sequence of edits between `a` and `b`
/// * `costs` - `EditCosts` struct for the cost of each edit operation
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::levenshtein::*;
/// let dist = levenshtein_naive_with_opts(b"abc", b"ab", true, LEVENSHTEIN_COSTS);
///
/// assert!(dist == (1, Some(vec![Edit{edit: EditType::Match, count: 2},
///                               Edit{edit: EditType::BGap, count: 1}])));
/// ```
#[inline]
pub fn levenshtein_naive_with_opts<T>(a: &[T], b: &[T], trace_on: bool, costs: EditCosts) -> (u32, Option<Vec<Edit>>)
    where T: PartialEq
{
    let swap = a.len() > b.len(); // swap so that a len <= b len
    let a_new = if swap {b} else {a};
    let a_new_len = a_new.len();
    let b_new = if swap {a} else {b};
    let b_new_len = b_new.len();
    let mismatch_cost = costs.mismatch_cost as u32;
    let gap_cost = costs.gap_cost as u32;
    let start_gap_cost = costs.start_gap_cost as u32;
    let transpose_cost = match costs.transpose_cost {
        Some(cost) => cost as u32,
        None => 0
    };
    let allow_transpose = costs.transpose_cost.is_some();

    let len = a_new_len + 1;
    let mut dp0 = vec![0u32; len];
    let mut dp1 = vec![0u32; len]; // in each iteration, dp0 and dp1 are already calculated
    let mut dp2 = vec![0u32; len]; // dp2 the currently calculated column
    let mut a_gap_dp = vec![u32::MAX; len];
    let mut b_gap_dp = vec![u32::MAX; len];
    let mut traceback = if trace_on {vec![0u8; (b_new_len + 1) * len]} else {vec![]};

    for i in 0..len {
        dp1[i] = (i as u32) * gap_cost + if i == 0 {0} else {start_gap_cost};

        if trace_on {
            traceback[0 * len + i] = 2u8;
        }
    }

    for i in 1..(b_new_len + 1) {
        a_gap_dp[0] = (i as u32) * gap_cost + start_gap_cost;
        dp2[0] = (i as u32) * gap_cost + start_gap_cost;

        if trace_on {
            traceback[i * len + 0] = 1u8;
        }

        for j in 1..len {
            let sub = dp1[j - 1] + ((a_new[j - 1] != b_new[i - 1]) as u32) * mismatch_cost;
            a_gap_dp[j] = cmp::min(dp1[j] + start_gap_cost + gap_cost, a_gap_dp[j].saturating_add(gap_cost));
            b_gap_dp[j] = cmp::min(dp2[j - 1] + start_gap_cost + gap_cost, b_gap_dp[j - 1].saturating_add(gap_cost));
            let traceback_idx = i * len + j;

            dp2[j] = a_gap_dp[j];

            if trace_on {
                traceback[traceback_idx] = 1u8;
            }

            if b_gap_dp[j] < dp2[j] {
                dp2[j] = b_gap_dp[j];

                if trace_on {
                    traceback[traceback_idx] = 2u8;
                }
            }

            if sub <= dp2[j] {
                dp2[j] = sub;

                if trace_on {
                    traceback[traceback_idx] = 0u8;
                }
            }

            if allow_transpose && i > 1 && j > 1
                && a_new[j - 1] == b_new[i - 2] && a_new[j - 2] == b_new[i - 1] {
                let transpose = dp0[j - 2] + transpose_cost;

                if transpose <= dp2[j] {
                    dp2[j] = transpose;

                    if trace_on {
                        traceback[traceback_idx] = 3u8;
                    }
                }
            }
        }

        mem::swap(&mut dp0, &mut dp1);
        mem::swap(&mut dp1, &mut dp2);
    }

    if trace_on {
        // estimate an upper bound for the number of Edits
        let mut upper_bound_edits = dp1[a_new_len] / cmp::min(mismatch_cost, gap_cost);

        if allow_transpose {
            upper_bound_edits = cmp::max(upper_bound_edits, (dp1[a_new_len] >> 1) / transpose_cost + 1);
        }

        let mut res: Vec<Edit> = Vec::with_capacity(((upper_bound_edits << 1) + 1) as usize);
        let mut i = b_new_len;
        let mut j = a_new_len;

        while i > 0 || j > 0 {
            let edit = traceback[i * len + j];

            let e = match edit {
                0 => {
                    i -= 1;
                    j -= 1;
                    if a_new[j] == b_new[i] {EditType::Match} else {EditType::Mismatch}
                },
                1 => {
                    i -= 1;
                    if swap {EditType::BGap} else {EditType::AGap}
                },
                2 => {
                    j -= 1;
                    if swap {EditType::AGap} else {EditType::BGap}
                },
                3 => {
                    i -= 2;
                    j -= 2;
                    EditType::Transpose
                },
                _ => unreachable!()
            };

            if res.len() > 0 && res.last().unwrap().edit == e {
                res.last_mut().unwrap().count += 1;
            }else{
                res.push(Edit{edit: e, count: 1});
            }
        }

        res.reverse();
        (dp1[a_new_len], Some(res))
    }else{
        (dp1[a_new_len], None)
    }
}

/// Returns the Levenshtein distance, bounded by a cost threshold `k`, between two strings, using the
/// naive scalar algorithm.
///
/// This will return `None` if the Levenshtein distance between `a` and `b` is greater than the
/// threshold `k`.
/// This should be much faster than `levenshtein_naive` if `k` is small compared to the lengths of
/// `a` and `b`.
///
/// # Arguments
/// * `a` - first string (slice)
/// * `b` - second string (slice)
/// * `k` - maximum number of edits allowed between `a` and `b`
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::levenshtein::*;
/// let dist = levenshtein_naive_k(b"abc", b"ab", 1);
///
/// assert!(dist.unwrap() == 1);
/// ```
pub fn levenshtein_naive_k(a: &[u8], b: &[u8], k: u32) -> Option<u32> {
    let res = levenshtein_naive_k_with_opts(a, b, k, false, LEVENSHTEIN_COSTS);

    match res {
        Some((edits, _)) => Some(edits),
        None => None
    }
}

/// Returns the Levenshtein distance, bounded by a cost threshold `k`, between two strings and optionally,
/// the edit traceback, using the naive scalar algorithm, with extra options.
///
/// This will return `None` if the Levenshtein distance between `a` and `b` is greater than the
/// threshold `k`.
/// This should be much faster than `levenshtein_naive_with_opts` if `k` is small compared to the lengths of
/// `a` and `b`.
///
/// # Arguments
/// * `a` - first string (slice)
/// * `b` - second string (slice)
/// * `k` - maximum number of cost allowed between `a` and `b`
/// * `trace_on` - whether to return the traceback, the sequence of edits between `a` and `b`
/// * `costs` - `EditCosts` struct for the cost of each edit operation
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::levenshtein::*;
/// let dist = levenshtein_naive_k_with_opts(b"abc", b"ab", 1, true, LEVENSHTEIN_COSTS);
///
/// assert!(dist.unwrap() == (1, Some(vec![Edit{edit: EditType::Match, count: 2},
///                                        Edit{edit: EditType::BGap, count: 1}])));
/// ```
#[inline]
pub fn levenshtein_naive_k_with_opts<T>(a: &[T], b: &[T], k: u32, trace_on: bool, costs: EditCosts) -> Option<(u32, Option<Vec<Edit>>)>
    where T: PartialEq
{
    let swap = a.len() > b.len(); // swap so that a len <= b len
    let a_new = if swap {b} else {a};
    let a_new_len = a_new.len();
    let b_new = if swap {a} else {b};
    let b_new_len = b_new.len();
    let mismatch_cost = costs.mismatch_cost as u32;
    let gap_cost = costs.gap_cost as u32;
    let start_gap_cost = costs.start_gap_cost as u32;
    let transpose_cost = match costs.transpose_cost {
        Some(cost) => cost as u32,
        None => 0
    };
    let allow_transpose = costs.transpose_cost.is_some();
    // upper bound on the number of edits, in case k is too large
    let max_k = cmp::min((a_new_len as u32) * mismatch_cost,
                         ((a_new_len as u32) << 1) * gap_cost
                         + if a_new_len == 0 {0} else {start_gap_cost
                             + if b_new_len == a_new_len {start_gap_cost} else {0}});
    let max_k = cmp::min(k,
                         max_k + ((b_new_len - a_new_len) as u32) * gap_cost
                         + if b_new_len == a_new_len {0} else {start_gap_cost});
    // farthest we can stray from the main diagonal
    // have to start at least one gap
    let unit_k = (max_k.saturating_sub(start_gap_cost) / gap_cost) as usize;

    if b_new_len - a_new_len > unit_k {
        return None;
    }

    let len = a_new_len + 1;
    let mut lo = 0usize;
    let mut hi = cmp::min(unit_k + 1, b_new_len + 1);
    let mut prev_lo0;
    let mut prev_lo1 = 0; // unused value
    let mut prev_hi;
    let k_len = cmp::min((unit_k << 1) + 1, b_new_len + 1);
    let mut dp0 = vec![0u32; k_len];
    let mut dp1 = vec![0u32; k_len]; // in each iteration, dp0 and dp1 are already calculated
    let mut dp2 = vec![0u32; k_len]; // dp2 the currently calculated row
    let mut a_gap_dp = vec![u32::MAX; k_len];
    let mut b_gap_dp = vec![u32::MAX; k_len];
    let mut traceback = if trace_on {vec![0u8; len * k_len]} else {vec![]};

    for i in 0..(hi - lo) {
        dp1[i] = (i as u32) * gap_cost + if i == 0 {0} else {start_gap_cost};

        if trace_on {
            traceback[0 * k_len + i] = 1u8;
        }
    }

    for i in 1..len {
        // keep track of prev_lo for the offset of previous rows
        prev_lo0 = prev_lo1;
        prev_lo1 = lo;
        prev_hi = hi;
        hi = cmp::min(hi + 1, b_new_len + 1);

        if i > unit_k {
            lo += 1;
        }

        for j in 0..(hi - lo) {
            let idx = lo + j;
            let sub = if idx == 0 {
                u32::MAX
            }else{
                dp1[idx - 1 - prev_lo1] + ((a_new[i - 1] != b_new[idx - 1]) as u32) * mismatch_cost
            };
            a_gap_dp[j] = if j == 0 {
                u32::MAX
            }else{
                cmp::min(dp2[j - 1] + start_gap_cost + gap_cost, a_gap_dp[j - 1].saturating_add(gap_cost))
            };
            b_gap_dp[j] = if idx >= prev_hi {
                u32::MAX
            }else{
                cmp::min(dp1[idx - prev_lo1] + start_gap_cost + gap_cost, b_gap_dp[idx - prev_lo1].saturating_add(gap_cost))
            };

            dp2[j] = sub;

            let traceback_idx = i * k_len + j;

            if trace_on {
                traceback[traceback_idx] = 0u8;
            }

            if a_gap_dp[j] < dp2[j] {
                dp2[j] = a_gap_dp[j];

                if trace_on {
                    traceback[traceback_idx] = 1u8;
                }
            }

            if b_gap_dp[j] < dp2[j] {
                dp2[j] = b_gap_dp[j];

                if trace_on {
                    traceback[traceback_idx] = 2u8;
                }
            }

            if allow_transpose && i > 1 && idx > 1
                && a_new[i - 1] == b_new[idx - 2] && a_new[i - 2] == b_new[idx - 1] {
                let transpose = dp0[idx - prev_lo0 - 2] + transpose_cost;

                if transpose <= dp2[j] {
                    dp2[j] = transpose;

                    if trace_on {
                        traceback[traceback_idx] = 3u8;
                    }
                }
            }
        }

        mem::swap(&mut dp0, &mut dp1);
        mem::swap(&mut dp1, &mut dp2);
    }

    if dp1[hi - lo - 1] > max_k {
        return None;
    }

    if !trace_on {
        return Some((dp1[hi - lo - 1], None));
    }

    // estimate an upper bound for the number of Edits
    let mut upper_bound_edits = dp1[hi - lo - 1] / cmp::min(mismatch_cost, gap_cost);

    if allow_transpose {
        upper_bound_edits = cmp::max(upper_bound_edits, (dp1[hi - lo - 1] >> 1) / transpose_cost + 1);
    }

    let mut res: Vec<Edit> = Vec::with_capacity(((upper_bound_edits << 1) + 1) as usize);
    let mut i = a_new_len;
    let mut j = b_new_len;

    while i > 0 || j > 0 {
        let edit = traceback[i * k_len + (j - (if i > unit_k {i - unit_k} else {0}))];

        let e = match edit {
            0 => {
                i -= 1;
                j -= 1;
                if a_new[i] == b_new[j] {EditType::Match} else {EditType::Mismatch}
            },
            1 => {
                j -= 1;
                if swap {EditType::BGap} else {EditType::AGap}
            },
            2 => {
                i -= 1;
                if swap {EditType::AGap} else {EditType::BGap}
            },
            3 => {
                i -= 2;
                j -= 2;
                EditType::Transpose
            },
            _ => unreachable!()
        };

        if res.len() > 0 && res.last().unwrap().edit == e {
            res.last_mut().unwrap().count += 1;
        }else{
            res.push(Edit{edit: e, count: 1});
        }
    }

    res.reverse();
    Some((dp1[hi - lo - 1], Some(res)))
}

fn translate_str(chars: &mut Vec<char>, s: &str) -> Option<Vec<u8>> {
    s.chars().map(|c| match chars.iter().position(|&d| c == d) {
        Some(i) => Some(i as u8),
        None => {
            let idx = chars.len();
            if idx < 256 {
                chars.push(c);
                Some(idx as u8)
            } else {
                None
            }
        }
    }).collect()
}

/// Returns the Levenshtein distance, bounded by a cost threshold `k`, between two utf8 encoded strings, using
/// SIMD acceleration.
/// # Arguments
/// * `a` - first string (&str)
/// * `b` - second string (&str)
/// * `k` - maximum number of edits allowed between `a` and `b`
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::levenshtein::*;
/// let dist = levenshtein_simd_k_str("abc", "ab", 1);
///
/// assert!(dist.unwrap() == 1);
/// ```
pub fn levenshtein_simd_k_str(a: &str, b: &str, k: u32) -> Option<u32> {
    if a.is_ascii() && b.is_ascii() {
        levenshtein_simd_k(a.as_bytes(), b.as_bytes(), k)
    } else {
        let mut chars = Vec::with_capacity(256);

        let a = translate_str(&mut chars, a)?;
        let b = translate_str(&mut chars, b)?;
        levenshtein_simd_k(&a, &b, k)
    }
}

/// Returns the Levenshtein distance, bounded by a cost threshold `k`, between two strings, using
/// SIMD acceleration.
///
/// This will return `None` if the Levenshtein distance between `a` and `b` is greater than the
/// threshold `k`.
/// This should be much faster than `levenshtein_naive` and `levenshtein_naive_k`.
/// Internally, this will automatically use AVX or SSE vectors with 8-bit, 16-bit, or 32-bit elements
/// to represent anti-diagonals in the dynamic programming matrix for calculating Levenshtein distance.
/// If AVX2 or SSE4.1 is not supported, then this will automatically fall back to
/// `levenshtein_naive_k_with_opts`.
///
/// # Arguments
/// * `a` - first string (slice)
/// * `b` - second string (slice)
/// * `k` - maximum number of edits allowed between `a` and `b`
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::levenshtein::*;
/// let dist = levenshtein_simd_k(b"abc", b"ab", 1);
///
/// assert!(dist.unwrap() == 1);
/// ```
pub fn levenshtein_simd_k(a: &[u8], b: &[u8], k: u32) -> Option<u32> {
    let res = levenshtein_simd_k_with_opts(a, b, k, false, LEVENSHTEIN_COSTS);

    match res {
        Some((edits, _)) => Some(edits),
        None => None
    }
}

/// Returns the Levenshtein distance, bounded by a cost threshold `k`, between two strings and optionally,
/// the edit traceback, using SIMD acceleration, with extra options.
///
/// This will return `None` if the Levenshtein distance between `a` and `b` is greater than the
/// threshold `k`.
/// This should be much faster than `levenshtein_naive_with_opts` and
/// `levenshtein_naive_k_with_opts`.
/// Internally, this will automatically use AVX or SSE vectors with 8-bit, 16-bit, or 32-bit elements
/// to represent anti-diagonals in the dynamic programming matrix for calculating Levenshtein distance.
/// If AVX2 or SSE4.1 is not supported, then this will automatically fall back to
/// `levenshtein_naive_k_with_opts`.
///
/// # Arguments
/// * `a` - first string (slice)
/// * `b` - second string (slice)
/// * `k` - maximum number of cost allowed between `a` and `b`
/// * `trace_on` - whether to return the traceback, the sequence of edits between `a` and `b`
/// * `costs` - `EditCosts` struct for the cost of each edit operation
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::levenshtein::*;
/// let dist = levenshtein_simd_k_with_opts(b"abc", b"ab", 1, true, LEVENSHTEIN_COSTS);
///
/// assert!(dist.unwrap() == (1, Some(vec![Edit{edit: EditType::Match, count: 2},
///                                        Edit{edit: EditType::BGap, count: 1}])));
/// ```
pub fn levenshtein_simd_k_with_opts(a: &[u8], b: &[u8], k: u32, trace_on: bool, costs: EditCosts) -> Option<(u32, Option<Vec<Edit>>)> {
    if a.len() == 0 && b.len() == 0 {
        return if trace_on {Some((0u32, Some(vec![])))} else {Some((0u32, None))};
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        let min_len = cmp::min(a.len(), b.len()) as u32;
        let max_len = cmp::max(a.len(), b.len()) as u32;
        // upper bound on the number of edits, in case k is too large
        let max_k = cmp::min(min_len * (costs.mismatch_cost as u32),
                             (min_len << 1) * (costs.gap_cost as u32)
                             + if min_len == 0 {0} else {costs.start_gap_cost as u32
                                 + if max_len == min_len {costs.start_gap_cost as u32} else {0}});
        let max_k = cmp::min(k,
                             max_k + (max_len - min_len) * (costs.gap_cost as u32)
                             + if max_len == min_len {0} else {costs.start_gap_cost as u32});
        // farthest we can stray from the main diagonal
        // have to start at least one gap
        let unit_k = cmp::min(max_k.saturating_sub(costs.start_gap_cost as u32) / (costs.gap_cost as u32), max_len);

        // note: do not use the MAX value, because it indicates overflow/inaccuracy
        if cfg!(feature = "jewel-avx") && is_x86_feature_detected!("avx2") {
            if cfg!(feature = "jewel-8bit") && unit_k <= (Avx1x32x8::static_upper_bound() as u32 - 2) && max_k <= ((u8::MAX - 1) as u32) {
                return unsafe {levenshtein_simd_core_avx_1x32x8(a, b, max_k, trace_on, costs)};
            }else if cfg!(feature = "jewel-8bit") && unit_k <= (Avx2x32x8::static_upper_bound() as u32 - 2) && max_k <= ((u8::MAX - 1) as u32) {
                return unsafe {levenshtein_simd_core_avx_2x32x8(a, b, max_k, trace_on, costs)};
            }else if cfg!(feature = "jewel-8bit") && unit_k <= (Avx4x32x8::static_upper_bound() as u32 - 2) && max_k <= ((u8::MAX - 1) as u32) {
                return unsafe {levenshtein_simd_core_avx_4x32x8(a, b, max_k, trace_on, costs)};
            }else if cfg!(feature = "jewel-8bit") && unit_k <= (Avx8x32x8::static_upper_bound() as u32 - 2) && max_k <= ((u8::MAX - 1) as u32) {
                return unsafe {levenshtein_simd_core_avx_8x32x8(a, b, max_k, trace_on, costs)};
            }else if cfg!(feature = "jewel-16bit") && max_k <= ((u16::MAX - 1) as u32) {
                return unsafe {levenshtein_simd_core_avx_nx16x16(a, b, max_k, trace_on, costs)};
            }else if cfg!(feature = "jewel-32bit") {
                return unsafe {levenshtein_simd_core_avx_nx8x32(a, b, max_k, trace_on, costs)};
            }
        }else if cfg!(feature = "jewel-sse") && is_x86_feature_detected!("sse4.1") {
            if cfg!(feature = "jewel-8bit") && unit_k <= (Sse1x16x8::static_upper_bound() as u32 - 2) && max_k <= ((u8::MAX - 1) as u32) {
                return unsafe {levenshtein_simd_core_sse_1x16x8(a, b, max_k, trace_on, costs)};
            }else if cfg!(feature = "jewel-8bit") && unit_k <= (Sse2x16x8::static_upper_bound() as u32 - 2) && max_k <= ((u8::MAX - 1) as u32) {
                return unsafe {levenshtein_simd_core_sse_2x16x8(a, b, max_k, trace_on, costs)};
            }else if cfg!(feature = "jewel-8bit") && unit_k <= (Sse4x16x8::static_upper_bound() as u32 - 2) && max_k <= ((u8::MAX - 1) as u32) {
                return unsafe {levenshtein_simd_core_sse_4x16x8(a, b, max_k, trace_on, costs)};
            }else if cfg!(feature = "jewel-8bit") && unit_k <= (Sse8x16x8::static_upper_bound() as u32 - 2) && max_k <= ((u8::MAX - 1) as u32) {
                return unsafe {levenshtein_simd_core_sse_8x16x8(a, b, max_k, trace_on, costs)};
            }else if cfg!(feature = "jewel-8bit") && unit_k <= (Sse16x16x8::static_upper_bound() as u32 - 2) && max_k <= ((u8::MAX - 1) as u32) {
                return unsafe {levenshtein_simd_core_sse_16x16x8(a, b, max_k, trace_on, costs)};
            }else if cfg!(feature = "jewel-16bit") && max_k <= ((u16::MAX - 1) as u32) {
                return unsafe {levenshtein_simd_core_sse_nx8x16(a, b, max_k, trace_on, costs)};
            }else if cfg!(feature = "jewel-32bit") {
                return unsafe {levenshtein_simd_core_sse_nx4x32(a, b, max_k, trace_on, costs)};
            }
        }
    }

    levenshtein_naive_k_with_opts(a, b, k, trace_on, costs)
}

macro_rules! create_levenshtein_simd_core {
    ($name:ident, $traceback_name:ident, $jewel:ty, $target:literal) => {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        #[target_feature(enable = $target)]
        unsafe fn $name(a_old: &[u8], b_old: &[u8], k: u32, trace_on: bool, costs: EditCosts) -> Option<(u32, Option<Vec<Edit>>)> {
            #[cfg(feature = "debug")]
            {
                println!("Debug: Levenshtein Jewel vector type {} for target {}.", stringify!($jewel), stringify!($target));
            }

            // swap a and b so that a is shorter than b, if applicable
            // makes operations later on slightly easier, since length of a <= length of b
            let swap = a_old.len() > b_old.len();
            let a = if swap {b_old} else {a_old};
            let a_len = a.len();
            let b = if swap {a_old} else {b_old};
            let b_len = b.len();
            let unit_k = cmp::min((k.saturating_sub(costs.start_gap_cost as u32) / (costs.gap_cost as u32)) as usize, b_len);

            if b_len - a_len > unit_k {
                return None;
            }

            // initialized with max values
            // must use saturated additions afterwards to not overflow
            let mut dp1 = <$jewel>::repeating_max((unit_k + 2) as usize);
            let max_len = dp1.upper_bound();
            let mut dp2 = <$jewel>::repeating_max(max_len);
            let mut dp0 = <$jewel>::repeating_max(max_len);
            let mut dp_temp = <$jewel>::repeating_max(max_len);
            // dp0 -> dp_temp -> dp1 -> dp2 -> current diagonal

            // dp for whether to extend gap or start new gap
            let mut a_gap_dp = <$jewel>::repeating_max(max_len);
            let mut b_gap_dp = <$jewel>::repeating_max(max_len);

            // lengths of the (anti) diagonals
            // assumes max_len is even
            let k1 = max_len - 1;
            let k1_div2 = k1 >> 1;
            let k2 = max_len - 2;
            let k2_div2 = k2 >> 1;

            // set dp[0][0] = 0
            dp1.slow_insert(k1_div2, 0);
            // set dp[0][1] = start_gap_cost + gap_cost and dp[1][0] = start_gap_cost + gap_cost
            dp2.slow_insert(k2_div2 - 1, costs.start_gap_cost as u32 + costs.gap_cost as u32);
            dp2.slow_insert(k2_div2, costs.start_gap_cost as u32 + costs.gap_cost as u32);
            b_gap_dp.slow_insert(k2_div2 - 1, costs.start_gap_cost as u32 + costs.gap_cost as u32);
            a_gap_dp.slow_insert(k2_div2, costs.start_gap_cost as u32 + costs.gap_cost as u32);

            // a_k1_window and a_k2_window represent reversed portions of the string a
            // copy in half of k1/k2 number of characters
            // these characters are placed in the second half of b windows
            // since a windows are reversed, the characters are placed in reverse in the first half of b windows
            let mut a_k1_window = <$jewel>::repeating(0, max_len);
            a_k1_window.slow_loadu(k1_div2 - 1, a.as_ptr(), cmp::min(k1_div2, a_len), true);

            let mut b_k1_window = <$jewel>::repeating(0, max_len);
            b_k1_window.slow_loadu(k1_div2 + 1, b.as_ptr(), cmp::min(k1_div2, b_len), false);

            let mut a_k2_window = <$jewel>::repeating(0, max_len);
            a_k2_window.slow_loadu(k2_div2 - 1, a.as_ptr(), cmp::min(k2_div2, a_len), true);

            let mut b_k2_window = <$jewel>::repeating(0, max_len);
            b_k2_window.slow_loadu(k2_div2, b.as_ptr(), cmp::min(k2_div2, b_len), false);

            // used to keep track of the next characters to place in the windows
            let mut k1_idx = k1_div2 - 1;
            let mut k2_idx = k2_div2 - 1;

            let len_diff = b_len - a_len;
            let len = a_len + b_len + 1;
            let len_div2 = (len >> 1) + (len & 1);
            let ends_with_k2 = len & 1 == 0;
            // every diff between the length of a and b results in a shift from the main diagonal
            let final_idx = {
                if ends_with_k2 { // divisible by 2, ends with k2
                    k2_div2 + ((len_diff - 1) >> 1)
                }else{ // not divisible by 2, ends with k1
                    k1_div2 + (len_diff >> 1)
                }
            };

            // 0 = match/mismatch, 1 = a gap, 2 = b gap, 3 = transpose
            let mut traceback_arr = if trace_on {Vec::with_capacity(len + (len & 1))} else {vec![]};

            if trace_on {
                traceback_arr.push(<$jewel>::repeating(0, max_len));
                traceback_arr.push(<$jewel>::repeating(0, max_len));
                traceback_arr.get_unchecked_mut(1).slow_insert(k2_div2 - 1, 2);
                traceback_arr.get_unchecked_mut(1).slow_insert(k2_div2, 1);
            }

            // reusable constant
            let threes = <$jewel>::repeating(3, max_len);

            // used in calculations
            let mut sub = <$jewel>::repeating(0, max_len);
            let mut match_mask0 = <$jewel>::repeating(0, max_len);
            let mut match_mask1 = <$jewel>::repeating(0, max_len);
            let mut a_gap = <$jewel>::repeating(0, max_len);
            let mut b_gap = <$jewel>::repeating(0, max_len);
            let mut transpose = <$jewel>::repeating(0, max_len);

            let mismatch_cost = <$jewel>::repeating(costs.mismatch_cost as u32, max_len);
            let gap_cost = <$jewel>::repeating(costs.gap_cost as u32, max_len);
            let start_gap_cost = <$jewel>::repeating(costs.start_gap_cost as u32 + costs.gap_cost as u32, max_len);
            let transpose_cost = match costs.transpose_cost {
                Some(cost) => <$jewel>::repeating(cost as u32, max_len),
                None => <$jewel>::repeating(0, max_len) // value does not matter
            };
            let allow_transpose = costs.transpose_cost.is_some();

            // example: allow k = 2 edits for two strings of length 3
            //      -b--   
            // | xx */*
            // a  x /*/*
            // |    */*/ x
            // |     */* xx
            //
            // each (anti) diagonal is represented with '*' or '/'
            // '/', use k2 = 2
            // '*', use k1 = 3
            // 'x' represents cells not in the "traditional" dp array
            // these out of bounds dp cells are shown because they represent
            // a horizontal sliding window of length 5 (2 * k + 1)
            //
            // dp2 is one diagonal before current
            // dp1 is two diagonals before current
            // dp0 is four diagonals before current
            // dp0 is useful for transpositions
            // we are trying to calculate the "current" diagonal
            // note that a k1 '*' dp diagonal has its center cell on the main diagonal
            // in general, the diagonals are centered on the main diagonal
            // each diagonal is represented using a Jewel vector
            // each vector goes from bottom-left to top-right
            //
            // the a windows and b windows are queues of a fixed length
            // a is reversed, so that elementwise comparison can be done between a and b
            // this operation obtains the comparison of characters along the (anti) diagonal
            // if transpositions are allowed, then previous match_masks must be saved to calculate
            // a[i - 1] == b[j] and a[i] == b[j - 1]
            // for speed, transpositions are done by directly blending using the mask, without calculating
            // the minimum cost compared to the other edit operations
            //
            // example of moving the windows:
            // a windows: [5 4 3 2 1] -> [6 5 4 3 2] (right shift + insert)
            // b windows: [1 2 3 4 5] -> [2 3 4 5 6] (left shift + insert)
            //
            // initially:
            // a windows: [2 1 0 0 0]
            // b windows: [0 0 0 1 2]
            //
            // note that there will be left over cells not filled in the Jewel vector
            // this is because k1 and k2 are not long enough
            // all of these empty cells should be at the end of the SIMD vectors
            //
            // each iteration of the loop below results in processing both a k1 diagonal and a k2 diagonal
            // this could be done with an alternating state flag but it is unrolled for less branching
            //
            // note: in traditional dp array
            // dp[i][j] -> dp[i + 1][j] is a gap in string b
            // dp[i][j] -> dp[i][j + 1] is a gap in string a

            for _ in 1..len_div2 {
                // move indexes in strings forward
                k1_idx += 1;
                k2_idx += 1;

                // move windows for the strings a and b
                a_k1_window.shift_right_1_mut();

                if k1_idx < a_len {
                    a_k1_window.insert_first(*a.get_unchecked(k1_idx) as u32);
                }

                b_k1_window.shift_left_1_mut();

                if k1_idx < b_len {
                    b_k1_window.insert_last_1(*b.get_unchecked(k1_idx) as u32); // k1 - 1
                }

                a_k2_window.shift_right_1_mut();

                if k2_idx < a_len {
                    a_k2_window.insert_first(*a.get_unchecked(k2_idx) as u32);
                }

                b_k2_window.shift_left_1_mut();

                if k2_idx < b_len {
                    b_k2_window.insert_last_2(*b.get_unchecked(k2_idx) as u32); // k2 - 1
                }

                // (anti) diagonal that matches in the a and b windows
                <$jewel>::cmpeq(&a_k1_window, &b_k1_window, &mut match_mask1);
                <$jewel>::andnot(&match_mask1, &mismatch_cost, &mut sub);
                sub.adds_mut(&dp1);
                // cost of gaps in a
                // start new gap
                <$jewel>::adds(&dp2, &start_gap_cost, &mut a_gap);
                // continue gap
                a_gap_dp.adds_mut(&gap_cost);
                a_gap_dp.min_mut(&a_gap);
                a_gap_dp.shift_right_1_mut();
                a_gap_dp.insert_first_max();
                // cost of gaps in b
                // start new gap
                <$jewel>::adds(&dp2, &start_gap_cost, &mut b_gap);
                // continue gap
                b_gap_dp.adds_mut(&gap_cost);
                b_gap_dp.min_mut(&b_gap);

                if allow_transpose {
                    <$jewel>::shift_right_1(&match_mask0, &mut transpose); // reuse transpose, zeros shifted in
                    transpose.and_mut(&match_mask0);
                    // make sure that current matching locations are excluded
                    <$jewel>::andnot(&match_mask1, &transpose, &mut match_mask0); // reuse match_mask0 to represent transpose mask
                    <$jewel>::adds(&dp0, &transpose_cost, &mut transpose);
                }

                // min of the cost of all three edit operations
                if trace_on {
                    let mut args = <$jewel>::triple_argmin(&sub, &a_gap_dp, &b_gap_dp, &mut dp0);

                    if allow_transpose {
                        // blend using transpose mask
                        dp0.blendv_mut(&transpose, &match_mask0);
                        args.blendv_mut(&threes, &match_mask0);
                        mem::swap(&mut match_mask0, &mut match_mask1);
                    }

                    traceback_arr.push(args);
                }else{
                    <$jewel>::min(&a_gap_dp, &b_gap_dp, &mut dp0);
                    dp0.min_mut(&sub);

                    if allow_transpose {
                        // blend using transpose mask
                        dp0.blendv_mut(&transpose, &match_mask0);
                        mem::swap(&mut match_mask0, &mut match_mask1);
                    }
                }

                mem::swap(&mut dp0, &mut dp_temp);
                mem::swap(&mut dp_temp, &mut dp1);
                mem::swap(&mut dp1, &mut dp2);

                // (anti) diagonal that matches in the a and b windows
                <$jewel>::cmpeq(&a_k2_window, &b_k2_window, &mut match_mask1);
                <$jewel>::andnot(&match_mask1, &mismatch_cost, &mut sub);
                sub.adds_mut(&dp1);
                // cost of gaps in b
                // start new gap
                <$jewel>::adds(&dp2, &start_gap_cost, &mut b_gap);
                // continue gap
                b_gap_dp.adds_mut(&gap_cost);
                b_gap_dp.min_mut(&b_gap);
                b_gap_dp.shift_left_1_mut();
                b_gap_dp.insert_last_max(); // k1, shift in max value
                // cost of gaps in a
                // start new gap
                <$jewel>::adds(&dp2, &start_gap_cost, &mut a_gap);
                a_gap_dp.adds_mut(&gap_cost);
                // continue gap
                a_gap_dp.min_mut(&a_gap);

                if allow_transpose {
                    <$jewel>::shift_left_1(&match_mask0, &mut transpose); // reuse transpose, zeros shifted in
                    transpose.and_mut(&match_mask0);
                    // make sure that current matching locations are excluded
                    <$jewel>::andnot(&match_mask1, &transpose, &mut match_mask0); // reuse match_mask0 to represent transpose mask
                    <$jewel>::adds(&dp0, &transpose_cost, &mut transpose);
                }

                // min of the cost of all three edit operations
                if trace_on {
                    let mut args = <$jewel>::triple_argmin(&sub, &a_gap_dp, &b_gap_dp, &mut dp0);

                    if allow_transpose {
                        // blend using transpose mask
                        dp0.blendv_mut(&transpose, &match_mask0);
                        args.blendv_mut(&threes, &match_mask0);
                        mem::swap(&mut match_mask0, &mut match_mask1);
                    }

                    traceback_arr.push(args);
                }else{
                    <$jewel>::min(&a_gap_dp, &b_gap_dp, &mut dp0);
                    dp0.min_mut(&sub);

                    if allow_transpose {
                        // blend using transpose mask
                        dp0.blendv_mut(&transpose, &match_mask0);
                        mem::swap(&mut match_mask0, &mut match_mask1);
                    }
                }

                mem::swap(&mut dp0, &mut dp_temp);
                mem::swap(&mut dp_temp, &mut dp1);
                mem::swap(&mut dp1, &mut dp2);
            }

            let final_res = if ends_with_k2 {
                dp2.slow_extract(final_idx)
            }else{
                dp1.slow_extract(final_idx)
            };

            if final_res > k {
                return None;
            }

            if !trace_on {
                return Some((final_res, None));
            }

            // upper bound the number of edit operations, to reduce memory allocations for saving the traceback
            let mut upper_bound_edits = final_res / (cmp::min(costs.mismatch_cost, costs.gap_cost) as u32);

            if let Some(cost) = costs.transpose_cost {
                upper_bound_edits = cmp::max(upper_bound_edits, (final_res >> 1) / (cost as u32) + 1);
            }

            Some((final_res, Some($traceback_name(&traceback_arr, upper_bound_edits as usize, final_idx, a, b, swap, ends_with_k2))))
        }

        unsafe fn $traceback_name(arr: &[$jewel], k: usize, mut idx: usize, a: &[u8], b: &[u8], swap: bool, mut is_k2: bool) -> Vec<Edit> {
            // keep track of position in traditional dp array and strings
            let mut i = a.len(); // index in a
            let mut j = b.len(); // index in b

            // last diagonal may overshoot, so ignore it
            let mut arr_idx = arr.len() - 1 - (if is_k2 {0} else {1});
            let mut res: Vec<Edit> = Vec::with_capacity((k << 1) + 1);

            while arr_idx > 0 {
                // each Jewel vector in arr is only visited once, so extract (which is costly) is fine
                let edit = arr.get_unchecked(arr_idx).slow_extract(idx);

                let e = match edit {
                    0 => { // match/mismatch
                        arr_idx -= 2;
                        i -= 1;
                        j -= 1;
                        if *a.get_unchecked(i) == *b.get_unchecked(j) {EditType::Match} else {EditType::Mismatch}
                    },
                    1 => { // a gap
                        arr_idx -= 1;

                        if !is_k2 {
                            idx -= 1;
                        }

                        j -= 1;
                        is_k2 = !is_k2; // must account for alternating k1/k2 diagonals
                        if swap {EditType::BGap} else {EditType::AGap} // account for the swap in the beginning
                    },
                    2 => { // b gap
                        arr_idx -= 1;

                        if is_k2 {
                            idx += 1;
                        }

                        i -= 1;
                        is_k2 = !is_k2;
                        if swap {EditType::AGap} else {EditType::BGap}
                    },
                    3 => { // transpose
                        arr_idx -= 4;
                        i -= 2;
                        j -= 2;
                        EditType::Transpose
                    },
                    _ => unreachable!()
                };

                if res.len() > 0 && res.last().unwrap().edit == e {
                    res.last_mut().unwrap().count += 1;
                }else{
                    res.push(Edit{edit: e, count: 1});
                }
            }

            res.reverse();
            res
        }
    };
}

// create a version of the functions for each Jewel vector
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_simd_core!(levenshtein_simd_core_avx_1x32x8, traceback_avx_1x32x8, Avx1x32x8, "avx2");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_simd_core!(levenshtein_simd_core_avx_2x32x8, traceback_avx_2x32x8, Avx2x32x8, "avx2");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_simd_core!(levenshtein_simd_core_avx_4x32x8, traceback_avx_4x32x8, Avx4x32x8, "avx2");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_simd_core!(levenshtein_simd_core_avx_8x32x8, traceback_avx_8x32x8, Avx8x32x8, "avx2");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_simd_core!(levenshtein_simd_core_avx_nx16x16, traceback_avx_nx16x16, AvxNx16x16, "avx2");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_simd_core!(levenshtein_simd_core_avx_nx8x32, traceback_avx_nx8x32, AvxNx8x32, "avx2");

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_simd_core!(levenshtein_simd_core_sse_1x16x8, traceback_sse_1x16x8, Sse1x16x8, "sse4.1");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_simd_core!(levenshtein_simd_core_sse_2x16x8, traceback_sse_2x16x8, Sse2x16x8, "sse4.1");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_simd_core!(levenshtein_simd_core_sse_4x16x8, traceback_sse_4x16x8, Sse4x16x8, "sse4.1");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_simd_core!(levenshtein_simd_core_sse_8x16x8, traceback_sse_8x16x8, Sse8x16x8, "sse4.1");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_simd_core!(levenshtein_simd_core_sse_16x16x8, traceback_sse_16x16x8, Sse16x16x8, "sse4.1");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_simd_core!(levenshtein_simd_core_sse_nx8x16, traceback_sse_nx8x16, SseNx8x16, "sse4.1");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_simd_core!(levenshtein_simd_core_sse_nx4x32, traceback_sse_nx4x32, SseNx4x32, "sse4.1");

/// Returns the Levenshtein distance between two strings using SIMD acceleration.
///
/// Note that `levenshtein_exp` may be much faster if the number of edits between the two strings
/// is expected to be small.
/// Internally, this will call `levenshtein_simd_k`.
/// If AVX2 or SSE4.1 is not supported, then this will automatically fall back to a scalar alternative.
///
/// # Arguments
/// * `a` - first string (slice)
/// * `b` - second string (slice)
///
/// # Example
/// ```
/// # use triple_accel::*;
/// let dist = levenshtein(b"abc", b"ab");
///
/// assert!(dist == 1);
/// ```
pub fn levenshtein(a: &[u8], b: &[u8]) -> u32 {
    levenshtein_simd_k(a, b, u32::MAX).unwrap()
}

/// Returns the restricted Damerau-Levenshtein distance between two strings using SIMD acceleration.
///
/// Note that `rdamerau_exp` may be much faster if the number of edits between the two strings
/// is expected to be small.
/// Internally, this will call `levenshtein_simd_k_with_opts`.
/// If AVX2 or SSE4.1 is not supported, then this will automatically fall back to a scalar alternative.
///
/// # Arguments
/// * `a` - first string (slice)
/// * `b` - second string (slice)
///
/// # Example
/// ```
/// # use triple_accel::*;
/// let dist = rdamerau(b"abc", b"acb");
///
/// assert!(dist == 1);
/// ```
pub fn rdamerau(a: &[u8], b: &[u8]) -> u32 {
    levenshtein_simd_k_with_opts(a, b, u32::MAX, false, RDAMERAU_COSTS).unwrap().0
}

/// Returns the Levenshtein distance between two strings using exponential search and SIMD
/// acceleration.
///
/// This may be much more efficient than `levenshtein` if the number of edits between `a` and `b`
/// is expected to be small.
/// Internally, this will call `levenshtein_simd_k` with values of `k` determined through
/// exponential search.
/// If AVX2 or SSE4.1 is not supported, then this will automatically fall back to a scalar alternative.
///
/// # Arguments
/// * `a` - first string (slice)
/// * `b` - second string (slice)
///
/// # Example
/// ```
/// # use triple_accel::*;
/// let dist = levenshtein_exp(b"abc", b"ab");
///
/// assert!(dist == 1);
/// ```
pub fn levenshtein_exp(a: &[u8], b: &[u8]) -> u32 {
    let mut k = 30;
    let mut res = levenshtein_simd_k(a, b, k);

    // exponential search
    while res.is_none() {
        k <<= 1;
        res = levenshtein_simd_k(a, b, k);
    }

    // should not panic
    res.unwrap()
}

/// Returns the restricted Damerau-Levenshtein distance between two strings using exponential
/// search and SIMD acceleration.
///
/// This may be much more efficient than `rdamerau` if the number of edits between `a` and `b`
/// is expected to be small.
/// Internally, this will call `levenshtein_simd_k_with_opts` with values of `k` determined through
/// exponential search.
/// If AVX2 or SSE4.1 is not supported, then this will automatically fall back to a scalar alternative.
///
/// # Arguments
/// * `a` - first string (slice)
/// * `b` - second string (slice)
///
/// # Example
/// ```
/// # use triple_accel::*;
/// let dist = rdamerau_exp(b"abc", b"acb");
///
/// assert!(dist == 1);
/// ```
pub fn rdamerau_exp(a: &[u8], b: &[u8]) -> u32 {
    let mut k = 30;
    let mut res = levenshtein_simd_k_with_opts(a, b, k, false, RDAMERAU_COSTS);

    // exponential search
    while res.is_none() {
        k <<= 1;
        res = levenshtein_simd_k_with_opts(a, b, k, false, RDAMERAU_COSTS);
    }

    // should not panic
    res.unwrap().0
}

/// Returns an iterator over the best `Match`s by searching through the text `haystack` for the
/// pattern `needle` using the naive algorithm.
///
/// The best matches are the matches with the lowest Levenshtein distance.
/// If multiple best matches end at the same position or fully overlap, then the longest match is chosen.
/// If `needle` is empty, then no `Match`es are returned.
/// Each returned `Match` requires at least half or more bytes of the `needle` to match
/// somewhere in the `haystack`.
///
/// # Arguments
/// * `needle` - pattern string (slice)
/// * `haystack` - text string (slice)
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::levenshtein::*;
/// let matches: Vec<Match> = levenshtein_search_naive(b"abc", b"  abd").collect();
///
/// assert!(matches == vec![Match{start: 2, end: 5, k: 1}]);
/// ```
pub fn levenshtein_search_naive<'a>(needle: &'a [u8], haystack: &'a [u8]) -> Box<dyn Iterator<Item = Match> + 'a> {
    levenshtein_search_naive_with_opts(needle, haystack, ((needle.len() as u32) >> 1) + ((needle.len() as u32) & 1), SearchType::Best, LEVENSHTEIN_COSTS, false)
}

/// Returns an iterator over `Match`s by searching through the text `haystack` for the
/// pattern `needle` using the naive algorithm, with extra options.
///
/// Note that overlapping matches may be returned.
/// If multiple matches end at the same position, then the longest match is chosen.
/// If `needle` is empty and `anchored` is false, then no `Match`es are returned.
///
/// # Arguments
/// * `needle` - pattern string (slice)
/// * `haystack` - text string (slice)
/// * `k` - maximum cost threshold for a match to be returned
/// * `search_type` - indicates whether to return all matches (within a cost of `k`), or the best matches with
/// the lowest cost (additionally, only the longest matches are retained for matches that fully overlap)
/// * `costs` - `EditCosts` struct for the cost of each edit operation
/// * `anchored` - whether the `needle` should be anchored to the start of the `haystack` string,
/// causing any shifts to cost gap edits
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::levenshtein::*;
/// let matches: Vec<Match> = levenshtein_search_naive_with_opts(b"abc", b"  acb", 1, SearchType::All, RDAMERAU_COSTS, false).collect();
///
/// // note: it is possible to end the match at two different positions
/// assert!(matches == vec![Match{start: 2, end: 4, k: 1}, Match{start: 2, end: 5, k: 1}]);
/// ```
pub fn levenshtein_search_naive_with_opts<'a>(needle: &'a [u8], haystack: &'a [u8], k: u32, search_type: SearchType, costs: EditCosts, anchored: bool) -> Box<dyn Iterator<Item = Match> + 'a> {
    let needle_len = needle.len();
    let haystack_len = haystack.len();

    if needle_len == 0 {
        // special case when anchored is true: return possible matches
        if anchored {
            return match search_type {
                SearchType::All => {
                    let mut i = 0;
                    let mut cost = costs.start_gap_cost as u32;
                    let mut start = true;

                    Box::new(iter::from_fn(move || {
                        if start {
                            start = false;
                            return Some(Match{start: 0, end: 0, k: 0});
                        }

                        if i < haystack_len {
                            i += 1;
                            cost += costs.gap_cost as u32;

                            if cost <= k {
                                return Some(Match{start: 0, end: i, k: cost});
                            }
                        }

                        None
                    }))
                },
                SearchType::Best => Box::new(iter::once(Match{start: 0, end: 0, k: 0}))
            };
        }else{
            return Box::new(iter::empty());
        }
    }

    // enforce another constraint on the costs
    costs.check_search();

    let len = needle_len + 1;
    let iter_len = if anchored {
        // start_gap_cost must be incurred at least once
        cmp::min(haystack_len, needle_len.saturating_add(
                (k.saturating_sub(costs.start_gap_cost as u32) as usize) / (costs.gap_cost as usize)))
    }else{
        haystack_len
    };

    let mut dp0 = vec![0u32; len];
    let mut dp1 = vec![0u32; len];
    let mut dp2 = vec![0u32; len];
    let mut needle_gap_dp = vec![u32::MAX; len];
    let mut haystack_gap_dp = vec![u32::MAX; len];
    let mut length0 = vec![0usize; len];
    let mut length1 = vec![0usize; len];
    let mut length2 = vec![0usize; len];
    let mut needle_gap_length = vec![0usize; len];
    let mut haystack_gap_length = vec![0usize; len];
    let mut curr_k = k;
    let mismatch_cost = costs.mismatch_cost as u32;
    let gap_cost = costs.gap_cost as u32;
    let start_gap_cost = costs.start_gap_cost as u32;
    let transpose_cost = match costs.transpose_cost {
        Some(cost) => cost as u32,
        None => 0
    };
    let allow_transpose = costs.transpose_cost.is_some();
    let mut first = true;
    let mut i = 0;

    let res = iter::from_fn(move || {
        if first {
            first = false;

            for j in 0..len {
                dp1[j] = (j as u32) * gap_cost + if j == 0 {0} else {start_gap_cost};
            }

            if dp1[len - 1] <= curr_k {
                if search_type == SearchType::Best {
                    curr_k = dp1[len - 1];
                }

                return Some((Match{start: 0, end: 0, k: dp1[len - 1]}, curr_k));
            }
        }

        while i < iter_len {
            needle_gap_dp[0] = if anchored {(i as u32 + 1) * gap_cost + start_gap_cost} else {0};
            dp2[0] = if anchored {(i as u32 + 1) * gap_cost + start_gap_cost} else {0};
            needle_gap_length[0] = 0;
            length2[0] = 0;

            for j in 1..len {
                let sub = dp1[j - 1] + ((needle[j - 1] != haystack[i]) as u32) * mismatch_cost;

                let new_gap = dp1[j] + start_gap_cost + gap_cost;
                let cont_gap = needle_gap_dp[j].saturating_add(gap_cost);
                if new_gap < cont_gap {
                    needle_gap_dp[j] = new_gap;
                    needle_gap_length[j] = length1[j] + 1;
                }else if new_gap > cont_gap {
                    needle_gap_dp[j] = cont_gap;
                    needle_gap_length[j] += 1;
                }else{
                    needle_gap_dp[j] = cont_gap;
                    needle_gap_length[j] = cmp::max(length1[j], needle_gap_length[j]) + 1;
                }

                let new_gap = dp2[j - 1] + start_gap_cost + gap_cost;
                let cont_gap = haystack_gap_dp[j - 1].saturating_add(gap_cost);
                if new_gap < cont_gap {
                    haystack_gap_dp[j] = new_gap;
                    haystack_gap_length[j] = length2[j - 1];
                }else if new_gap > cont_gap {
                    haystack_gap_dp[j] = cont_gap;
                    haystack_gap_length[j] = haystack_gap_length[j - 1];
                }else{
                    haystack_gap_dp[j] = cont_gap;
                    haystack_gap_length[j] = cmp::max(length2[j - 1], haystack_gap_length[j - 1]);
                }

                dp2[j] = needle_gap_dp[j];
                length2[j] = needle_gap_length[j];

                if (haystack_gap_dp[j] < dp2[j]) || (haystack_gap_dp[j] == dp2[j] && length2[j - 1] > length2[j]) {
                    dp2[j] = haystack_gap_dp[j];
                    length2[j] = haystack_gap_length[j];
                }

                if (sub < dp2[j]) || (sub == dp2[j] && (length1[j - 1] + 1) > length2[j]) {
                    dp2[j] = sub;
                    length2[j] = length1[j - 1] + 1;
                }

                if allow_transpose && i > 0 && j > 1
                    && needle[j - 1] == haystack[i - 1] && needle[j - 2] == haystack[i] {
                        let transpose = dp0[j - 2] + transpose_cost;

                        if transpose <= dp2[j] {
                            dp2[j] = transpose;
                            length2[j] = length0[j - 2] + 2;
                        }
                    }
            }

            let final_res = dp2[len - 1];
            let final_length = length2[len - 1];

            mem::swap(&mut dp0, &mut dp1);
            mem::swap(&mut dp1, &mut dp2);
            mem::swap(&mut length0, &mut length1);
            mem::swap(&mut length1, &mut length2);

            i += 1;

            if final_res <= curr_k {
                match search_type {
                    SearchType::Best => curr_k = final_res,
                    _ => ()
                }

                return Some((Match{start: i - final_length, end: i, k: final_res}, curr_k));
            }
        }

        None
    });

    if search_type == SearchType::Best {
        // estimate the number of Matches
        let mut res_vec = Vec::with_capacity(iter_len / needle_len);

        for m in res {
            match res_vec.len() {
                0 => res_vec.push(m.0),
                _ => {
                    let last = res_vec.last_mut().unwrap();

                    // replace previous if fully overlapping
                    if m.0.start <= last.start {
                        *last = m.0;
                    }else{
                        res_vec.push(m.0);
                    }
                }
            }

            curr_k = m.1;
        }

        return Box::new(res_vec.into_iter().filter(move |m| m.k == curr_k));
    }

    Box::new(res.map(|m| m.0))
}

/// Returns an iterator over the best `Match`s by searching through the text `haystack` for the
/// pattern `needle` using SIMD acceleration.
///
/// The best matches are the matches with the lowest Levenshtein distance.
/// If multiple best matches end at the same position or fully overlap, then the longest match is chosen.
/// If `needle` is empty, then no `Match`es are returned.
/// Each returned `Match` requires at least half or more bytes of the `needle` to match
/// somewhere in the `haystack`.
/// This should be much faster than `levenshtein_search_naive`.
/// Internally, this will automatically use AVX or SSE vectors with 8-bit, 16-bit, or 32-bit elements
/// to represent anti-diagonals in the dynamic programming matrix for calculating Levenshtein distance.
/// If AVX2 or SSE4.1 is not supported, then this will automatically fall back to
/// `levenshtein_search_naive_with_opts`.
///
/// # Arguments
/// * `needle` - pattern string (slice)
/// * `haystack` - text string (slice)
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::levenshtein::*;
/// let matches: Vec<Match> = levenshtein_search_simd(b"abc", b"  abd").collect();
///
/// assert!(matches == vec![Match{start: 2, end: 5, k: 1}]);
/// ```
pub fn levenshtein_search_simd<'a>(needle: &'a [u8], haystack: &'a [u8]) -> Box<dyn Iterator<Item = Match> + 'a> {
    levenshtein_search_simd_with_opts(needle, haystack, ((needle.len() >> 1) as u32) + ((needle.len() as u32) & 1), SearchType::Best, LEVENSHTEIN_COSTS, false)
}

/// Returns an iterator over `Match`s by searching through the text `haystack` for the
/// pattern `needle` using SIMD acceleration, with extra options.
///
/// Note that overlapping matches may be returned.
/// If multiple matches end at the same position, then the longest match is chosen.
/// If `needle` is empty and `anchored` is false, then no `Match`es are returned.
/// This should be much faster than `levenshtein_search_naive_with_opts`.
/// Internally, this will automatically use AVX or SSE vectors with 8-bit, 16-bit, or 32-bit elements
/// to represent anti-diagonals in the dynamic programming matrix for calculating Levenshtein distance.
/// If AVX2 or SSE4.1 is not supported, then this will automatically fall back to
/// `levenshtein_search_naive_with_opts`.
///
/// # Arguments
/// * `needle` - pattern string (slice)
/// * `haystack` - text string (slice)
/// * `k` - maximum cost threshold for a match to be returned
/// * `search_type` - indicates whether to return all matches (within a cost of `k`), or the best matches with
/// the lowest cost (additionally, only the longest matches are retained for matches that fully overlap)
/// * `costs` - `EditCosts` struct for the cost of each edit operation
/// * `anchored` - whether the `needle` should be anchored to the start of the `haystack` string,
/// causing any shifts to cost gap edits
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::levenshtein::*;
/// let matches: Vec<Match> = levenshtein_search_simd_with_opts(b"abc", b"  acb", 1, SearchType::All, RDAMERAU_COSTS, false).collect();
///
/// // note: it is possible to end the match at two different positions
/// assert!(matches == vec![Match{start: 2, end: 4, k: 1}, Match{start: 2, end: 5, k: 1}]);
/// ```
pub fn levenshtein_search_simd_with_opts<'a>(needle: &'a [u8], haystack: &'a [u8], k: u32, search_type: SearchType, costs: EditCosts, anchored: bool) -> Box<dyn Iterator<Item = Match> + 'a> {
    if needle.len() == 0 {
        // special case when anchored is true: return possible matches
        if anchored {
            return match search_type {
                SearchType::All => {
                    let mut i = 0;
                    let mut cost = costs.start_gap_cost as u32;
                    let mut start = true;

                    Box::new(iter::from_fn(move || {
                        if start {
                            start = false;
                            return Some(Match{start: 0, end: 0, k: 0});
                        }

                        if i < haystack.len() {
                            i += 1;
                            cost += costs.gap_cost as u32;

                            if cost <= k {
                                return Some(Match{start: 0, end: i, k: cost});
                            }
                        }

                        None
                    }))
                },
                SearchType::Best => Box::new(iter::once(Match{start: 0, end: 0, k: 0}))
            };
        }else{
            return Box::new(iter::empty());
        }
    }

    costs.check_search();

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        let unit_k = k.saturating_sub(costs.start_gap_cost as u32) / (costs.gap_cost as u32);
        // either the length of the match or the number of edits may exceed the maximum
        // available int size; additionally, MAX value is used to indicate overflow
        let upper_bound = cmp::max((needle.len() as u32).saturating_add(unit_k), k.saturating_add(1));

        if cfg!(feature = "jewel-avx") && is_x86_feature_detected!("avx2") {
            if cfg!(feature = "jewel-8bit") && needle.len() <= Avx1x32x8::static_upper_bound() && upper_bound <= u8::MAX as u32 {
                return unsafe {levenshtein_search_simd_core_avx_1x32x8(needle, haystack, k, search_type, costs, anchored)};
            }else if cfg!(feature = "jewel-8bit") && needle.len() <= Avx2x32x8::static_upper_bound() && upper_bound <= u8::MAX as u32 {
                return unsafe {levenshtein_search_simd_core_avx_2x32x8(needle, haystack, k, search_type, costs, anchored)};
            }else if cfg!(feature = "jewel-8bit") && needle.len() <= Avx4x32x8::static_upper_bound() && upper_bound <= u8::MAX as u32 {
                return unsafe {levenshtein_search_simd_core_avx_4x32x8(needle, haystack, k, search_type, costs, anchored)};
            }else if cfg!(feature = "jewel-8bit") && needle.len() <= Avx8x32x8::static_upper_bound() && upper_bound <= u8::MAX as u32 {
                return unsafe {levenshtein_search_simd_core_avx_8x32x8(needle, haystack, k, search_type, costs, anchored)};
            }else if cfg!(feature = "jewel-16bit") && upper_bound <= u16::MAX as u32 {
                return unsafe {levenshtein_search_simd_core_avx_nx16x16(needle, haystack, k, search_type, costs, anchored)};
            }else if cfg!(feature = "jewel-32bit") {
                return unsafe {levenshtein_search_simd_core_avx_nx8x32(needle, haystack, k, search_type, costs, anchored)};
            }
        }else if cfg!(feature = "jewel-sse") && is_x86_feature_detected!("sse4.1") {
            if cfg!(feature = "jewel-8bit") && needle.len() <= Sse1x16x8::static_upper_bound() && upper_bound <= u8::MAX as u32 {
                return unsafe {levenshtein_search_simd_core_sse_1x16x8(needle, haystack, k, search_type, costs, anchored)};
            }else if cfg!(feature = "jewel-8bit") && needle.len() <= Sse2x16x8::static_upper_bound() && upper_bound <= u8::MAX as u32 {
                return unsafe {levenshtein_search_simd_core_sse_2x16x8(needle, haystack, k, search_type, costs, anchored)};
            }else if cfg!(feature = "jewel-8bit") && needle.len() <= Sse4x16x8::static_upper_bound() && upper_bound <= u8::MAX as u32 {
                return unsafe {levenshtein_search_simd_core_sse_4x16x8(needle, haystack, k, search_type, costs, anchored)};
            }else if cfg!(feature = "jewel-8bit") && needle.len() <= Sse8x16x8::static_upper_bound() && upper_bound <= u8::MAX as u32 {
                return unsafe {levenshtein_search_simd_core_sse_8x16x8(needle, haystack, k, search_type, costs, anchored)};
            }else if cfg!(feature = "jewel-8bit") && needle.len() <= Sse16x16x8::static_upper_bound() && upper_bound <= u8::MAX as u32 {
                return unsafe {levenshtein_search_simd_core_sse_16x16x8(needle, haystack, k, search_type, costs, anchored)};
            }else if cfg!(feature = "jewel-16bit") && upper_bound <= u16::MAX as u32 {
                return unsafe {levenshtein_search_simd_core_sse_nx8x16(needle, haystack, k, search_type, costs, anchored)};
            }else if cfg!(feature = "jewel-32bit") {
                return unsafe {levenshtein_search_simd_core_sse_nx4x32(needle, haystack, k, search_type, costs, anchored)};
            }
        }
    }

    levenshtein_search_naive_with_opts(needle, haystack, k, search_type, costs, anchored)
}

macro_rules! create_levenshtein_search_simd_core {
    ($name:ident, $jewel:ty, $target:literal) => {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        #[target_feature(enable = $target)]
        unsafe fn $name<'a>(needle: &'a [u8], haystack: &'a [u8], k: u32, search_type: SearchType, costs: EditCosts, anchored: bool) -> Box<dyn Iterator<Item = Match> + 'a> {
            #[cfg(feature = "debug")]
            {
                println!("Debug: Levenshtein search Jewel vector type {} for target {}.", stringify!($jewel), stringify!($target));
            }

            let needle_len = needle.len();
            let haystack_len = haystack.len();
            let mut dp0 = <$jewel>::repeating_max(needle_len);
            let mut dp_temp = <$jewel>::repeating_max(needle_len);
            let mut dp1 = <$jewel>::repeating_max(needle_len);
            let mut dp2 = <$jewel>::repeating_max(needle_len);
            let mut needle_gap_dp = <$jewel>::repeating_max(needle_len);
            let mut haystack_gap_dp = <$jewel>::repeating_max(needle_len);
            dp2.slow_insert(dp2.upper_bound() - 1, costs.start_gap_cost as u32 + costs.gap_cost as u32); // last cell
            haystack_gap_dp.slow_insert(dp2.upper_bound() - 1, costs.start_gap_cost as u32 + costs.gap_cost as u32);

            // save length instead of start idx due to int size constraints
            let mut length0 = <$jewel>::repeating(0, needle_len);
            let mut length_temp = <$jewel>::repeating(0, needle_len);
            let mut length1 = <$jewel>::repeating(0, needle_len);
            let mut length2 = <$jewel>::repeating(0, needle_len);
            let mut needle_gap_length = <$jewel>::repeating(0, needle_len);
            let mut haystack_gap_length = <$jewel>::repeating(0, needle_len);

            let ones = <$jewel>::repeating(1, needle_len);
            let twos = <$jewel>::repeating(2, needle_len);

            // the suffix of haystack can be ignored if needle must be anchored
            let len = if anchored {
                needle_len + cmp::min(haystack_len, needle_len.saturating_add(
                        (k.saturating_sub(costs.start_gap_cost as u32) as usize) / (costs.gap_cost as usize)))
            }else{
                needle_len + haystack_len
            };

            let final_idx = dp1.upper_bound() - needle_len;

            // load needle characters into needle_window in reversed order
            let mut needle_window = <$jewel>::repeating(0, needle_len);
            needle_window.slow_loadu(needle_window.upper_bound() - 1, needle.as_ptr(), needle_len, true);

            let mut haystack_window = <$jewel>::repeating(0, needle_len);
            let mut haystack_idx = 0usize;
            let mut curr_k = k;

            // used in calculations
            let mut match_mask0 = <$jewel>::repeating(0, needle_len);
            let mut match_mask1 = <$jewel>::repeating(0, needle_len);
            let mut match_mask_cost = <$jewel>::repeating(0, needle_len);
            let mut sub = <$jewel>::repeating(0, needle_len);
            let mut sub_length = <$jewel>::repeating(0, needle_len);
            let mut needle_gap = <$jewel>::repeating(0, needle_len);
            let mut haystack_gap = <$jewel>::repeating(0, needle_len);
            let mut transpose = <$jewel>::repeating(0, needle_len);
            let mut transpose_length = <$jewel>::repeating(0, needle_len);

            let mismatch_cost = <$jewel>::repeating(costs.mismatch_cost as u32, needle_len);
            let gap_cost = <$jewel>::repeating(costs.gap_cost as u32, needle_len);
            let start_gap_cost = <$jewel>::repeating(costs.start_gap_cost as u32, needle_len);
            let transpose_cost = match costs.transpose_cost {
                Some(cost) => <$jewel>::repeating(cost as u32, needle_len),
                None => <$jewel>::repeating(0, needle_len)
            };
            let allow_transpose = costs.transpose_cost.is_some();

            //       ..i...
            //       --h---
            // |     //////xxxx
            // |    x//////xxx
            // n   xx//////xx
            // |  xxx//////x
            // | xxxx//////
            //
            // 'n' = needle, 'h' = haystack
            // each (anti) diagonal is denoted using '/' and 'x'
            // 'x' marks cells that are not in the traditional dp array
            // every (anti) diagonal is calculated simultaneously using Jewel vectors
            // note: each vector goes from bottom-left to top-right, ending at the first row in the
            // DP matrix
            // similar to levenshtein_simd_k, but without alternating anti-diagonals
            // note that the first row, which should be all zeros for searching, is not saved in the
            // Jewel vectors, for space concerns
            // therefore, starting at i = 1 does not include the first diagonal that only contains
            // a zero from the initial row of zeros
            // when left shift are required, then zeros must be shifted in
            // if anchored = true, then the number of gaps times the gap cost plus the starting gap
            // cost must be shifted in
            // for speed, transpositions are done by directly blending using the mask, without calculating
            // the minimum cost compared to the other edit operations

            let mut i = 1;

            let res = iter::from_fn(move || {
                while i < len {
                    // shift the haystack window
                    haystack_window.shift_left_1_mut();

                    if haystack_idx < haystack_len {
                        haystack_window.insert_last_0(*haystack.get_unchecked(haystack_idx) as u32);
                        haystack_idx += 1;
                    }

                    <$jewel>::cmpeq(&needle_window, &haystack_window, &mut match_mask1);
                    <$jewel>::andnot(&match_mask1, &mismatch_cost, &mut match_mask_cost);

                    // match/mismatch
                    <$jewel>::shift_left_1(&dp1, &mut sub);

                    if anchored && i > 1 {
                        // dp1 is 2 diagonals behind the current i
                        // must be capped at k to prevent overflow when inserting
                        sub.insert_last_0(cmp::min((i as u32 - 1) * (costs.gap_cost as u32) + costs.start_gap_cost as u32, k + 1));
                    }

                    sub.adds_mut(&match_mask_cost);

                    <$jewel>::shift_left_1(&length1, &mut sub_length); // zeros are shifted in
                    sub_length.add_mut(&ones);

                    // gap in needle
                    <$jewel>::adds(&dp2, &start_gap_cost, &mut needle_gap);
                    <$jewel>::double_min_length(&needle_gap, &mut needle_gap_dp, &length2, &mut needle_gap_length);
                    needle_gap_dp.adds_mut(&gap_cost);
                    needle_gap_length.add_mut(&ones);

                    // gap in haystack
                    <$jewel>::adds(&dp2, &start_gap_cost, &mut haystack_gap);
                    <$jewel>::double_min_length(&haystack_gap, &mut haystack_gap_dp, &length2, &mut haystack_gap_length);
                    haystack_gap_dp.shift_left_1_mut(); // zeros are shifted in

                    if anchored {
                        // dp2 is one diagonal behind the current i
                        haystack_gap_dp.insert_last_0(cmp::min((i as u32) * (costs.gap_cost as u32) + costs.start_gap_cost as u32, k + 1));
                    }else{
                        haystack_gap_dp.insert_last_0(costs.start_gap_cost as u32);
                    }

                    haystack_gap_dp.adds_mut(&gap_cost);
                    haystack_gap_length.shift_left_1_mut(); // zeros are shifted in

                    if allow_transpose {
                        <$jewel>::shift_left_1(&match_mask0, &mut transpose); // reuse transpose, shift in zeros
                        transpose.and_mut(&match_mask0);
                        // ensure that current matches are excluded
                        <$jewel>::andnot(&match_mask1, &transpose, &mut match_mask0); // reuse match_mask0 to represent transpose mask
                        dp0.shift_left_2_mut();

                        if anchored && i > 3 {
                            // dp0 is four diagonals behind the current i
                            dp0.insert_last_1(cmp::min((i as u32 - 3) * (costs.gap_cost as u32) + costs.start_gap_cost as u32, k + 1));
                        }

                        length0.shift_left_2_mut();
                        <$jewel>::adds(&dp0, &transpose_cost, &mut transpose);
                        <$jewel>::add(&length0, &twos, &mut transpose_length);
                    }

                    <$jewel>::triple_min_length(&sub, &needle_gap_dp, &haystack_gap_dp, &sub_length,
                                                &needle_gap_length, &haystack_gap_length, &mut dp0, &mut length0);

                    if allow_transpose {
                        // blend using transpose mask
                        dp0.blendv_mut(&transpose, &match_mask0);
                        length0.blendv_mut(&transpose_length, &match_mask0);
                        mem::swap(&mut match_mask0, &mut match_mask1);
                    }

                    mem::swap(&mut dp0, &mut dp_temp);
                    mem::swap(&mut dp_temp, &mut dp1);
                    mem::swap(&mut dp1, &mut dp2);
                    mem::swap(&mut length0, &mut length_temp);
                    mem::swap(&mut length_temp, &mut length1);
                    mem::swap(&mut length1, &mut length2);

                    i += 1;

                    if i >= needle_len {
                        let final_res = dp2.slow_extract(final_idx);
                        let final_length = length2.slow_extract(final_idx) as usize;

                        if final_res <= curr_k {
                            let end_idx = i - needle_len;

                            match search_type {
                                // if we want the best, then we can shrink the k threshold
                                SearchType::Best => curr_k = final_res,
                                _ => ()
                            }

                            return Some((Match{start: end_idx - final_length, end: end_idx, k: final_res}, curr_k));
                        }
                    }
                }

                None
            });

            if search_type == SearchType::Best {
                // estimate the number of Matches
                let mut res_vec = Vec::with_capacity((len - needle_len) / needle_len);

                for m in res {
                    match res_vec.len() {
                        0 => res_vec.push(m.0),
                        _ => {
                            let last = res_vec.last_mut().unwrap();

                            // replace previous if fully overlapping
                            if m.0.start <= last.start {
                                *last = m.0;
                            }else{
                                res_vec.push(m.0);
                            }
                        }
                    }

                    curr_k = m.1;
                }

                // only retain matches with the lowest k
                return Box::new(res_vec.into_iter().filter(move |m| m.k == curr_k));
            }

            Box::new(res.map(|m| m.0))
        }
    };
}

// duplicate functions for each Jewel vector type
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_search_simd_core!(levenshtein_search_simd_core_avx_1x32x8, Avx1x32x8, "avx2");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_search_simd_core!(levenshtein_search_simd_core_avx_2x32x8, Avx2x32x8, "avx2");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_search_simd_core!(levenshtein_search_simd_core_avx_4x32x8, Avx4x32x8, "avx2");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_search_simd_core!(levenshtein_search_simd_core_avx_8x32x8, Avx8x32x8, "avx2");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_search_simd_core!(levenshtein_search_simd_core_avx_nx16x16, AvxNx16x16, "avx2");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_search_simd_core!(levenshtein_search_simd_core_avx_nx8x32, AvxNx8x32, "avx2");

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_search_simd_core!(levenshtein_search_simd_core_sse_1x16x8, Sse1x16x8, "sse4.1");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_search_simd_core!(levenshtein_search_simd_core_sse_2x16x8, Sse2x16x8, "sse4.1");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_search_simd_core!(levenshtein_search_simd_core_sse_4x16x8, Sse4x16x8, "sse4.1");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_search_simd_core!(levenshtein_search_simd_core_sse_8x16x8, Sse8x16x8, "sse4.1");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_search_simd_core!(levenshtein_search_simd_core_sse_16x16x8, Sse16x16x8, "sse4.1");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_search_simd_core!(levenshtein_search_simd_core_sse_nx8x16, SseNx8x16, "sse4.1");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_levenshtein_search_simd_core!(levenshtein_search_simd_core_sse_nx4x32, SseNx4x32, "sse4.1");

/// Returns an iterator over best `Match`s by searching through the text `haystack` for the
/// pattern `needle` using SIMD acceleration.
///
/// The best matches are the matches with the lowest Levenshtein distance.
/// If multiple best matches end at the same position or fully overlap, then the longest match is chosen.
/// If `needle` is empty, then no `Match`es are returned.
/// Each returned `Match` requires at least half or more bytes of the `needle` to match
/// somewhere in the `haystack`.
/// Internally, this will call `levenshtein_search_simd`.
/// If AVX2 or SSE4.1 is not supported, then this will automatically fall back to a scalar alternative.
///
/// # Arguments
/// * `needle` - pattern string (slice)
/// * `haystack` - text string (slice)
///
/// # Example
/// ```
/// # use triple_accel::*;
/// let matches: Vec<Match> = levenshtein_search(b"abc", b"  abd").collect();
///
/// assert!(matches == vec![Match{start: 2, end: 5, k: 1}]);
/// ```
pub fn levenshtein_search<'a>(needle: &'a [u8], haystack: &'a [u8]) -> Box<dyn Iterator<Item = Match> + 'a> {
    levenshtein_search_simd(needle, haystack)
}

