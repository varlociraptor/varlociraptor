//! # triple_accel
//!
//! Rust edit distance routines accelerated using SIMD. Supports fast Hamming, Levenshtein,
//! restricted Damerau-Levenshtein, etc. distance calculations and string search.
//!
//! Although vectorized SIMD code allows for up to 20-30x speedups over their scalar counterparts,
//! the difficulty of handling platform-dependent SIMD code makes SIMD routines less attractive.
//! The goal of this library is to provide an easy-to-use abstraction over SIMD edit distance routines
//! that fall back to scalar routines if the target CPU architecture is not supported.
//! Additionally, all limitations and tradeoffs of the edit distance routines should be provided upfront
//! so the user knows exactly what to expect.
//! Finally, this library should lead to performance boosts on both short and longer strings, so it
//! can be used for a variety of tasks, from bioinformatics to natural language processing.
//! `triple_accel` is very lightweight: it only has dependencies on other crates for benchmarking.
//! It can be built on machines without CPUs that have AVX2 or SSE4.1 support. It can also run on
//! machines without SIMD support by automatically using scalar alternatives.
//!
//! ## Features
//!
//! This library provides routines for both searching for some needle string in a haystack string
//! and calculating the edit distance between two strings. Hamming distance (mismatches only),
//! Levenshtein distance (mismatches + gaps), and restricted Damerau-Levenshtein distance
//! (transpositions + mismatches + gaps) are supported, along with arbitrary edit costs. This
//! library provides a simple interface, in addition to powerful lower-level control over the edit
//! distance calculations.
//!
//! At runtime, the implementation for a certain algorithm is selected based on CPU support, going
//! down the list:
//!
//! 1. Vectorized implementation with 256-bit AVX vectors, if AVX2 is supported.
//! 2. Vectorized implementation with 128-bit SSE vectors, if SSE4.1 is supported.
//! 3. Scalar implementation.
//!
//! Currently, vectorized SIMD implementations are only available for x86 or x86-64 CPUs. However,
//! after compiling this library on a machine that supports those SIMD intrinsics, the library can
//! be used on other machines.
//! Additionally, the internal data structure for storing vectors and the bit width of the values
//! in the vectors are selected at runtime for maximum efficiency and accuracy, given the lengths
//! of the input strings.
//!
//! ## Limitations
//!
//! Due to the use of SIMD intrinsics, only binary strings that are represented with `u8` bytes
//! are supported. Unicode strings are not currently supported.
//!
//! ## Notation
//!
//! Quick notation notes that will often appear in the code/documentation:
//!
//! * `k` - the number of edits that are allowed
//! * `a` and `b` - any two strings; this is usually used for edit distance routines
//! * `needle` and `haystack` - any two strings; we want to search for where needle appears in
//! haystack
//!
//! ## Examples
//!
//! Calculating the Hamming distance (number of mismatches) between two strings is extremely simple:
//! ```
//! use triple_accel::*;
//!
//! let a = b"abcd";
//! let b = b"abcc";
//!
//! let dist = hamming(a, b);
//! assert!(dist == 1);
//! ```
//! By default, SIMD will be used if possible. Similarly, we can easily calculate the Levenshtein
//! distance (character mismatches and gaps all have a cost of 1) between two strings with the
//! following code:
//! ```
//! use triple_accel::*;
//!
//! let a = b"abc";
//! let b = b"abcd";
//!
//! let dist = levenshtein_exp(a, b);
//! assert!(dist == 1);
//! ```
//! This uses exponential search to estimate the number of edits between `a` and `b`, which makes it
//! more efficient than the alternative `levenshtein` function when the number of edits between `a`
//! and `b` is low.
//!
//! In addition to edit distance routines, `triple_accel` also provides search routines. These
//! routines return an iterator over matches that indicate where the `needle` string matches the `haystack`
//! string. `triple_accel` will attempt to maximize the length of matches that end at the same position and
//! remove shorter matches when some matches fully overlap.
//! ```
//! use triple_accel::*;
//!
//! let needle = b"helllo";
//! let haystack = b"hello world";
//!
//! let matches: Vec<Match> = levenshtein_search(needle, haystack).collect();
//! // note: start index is inclusive, end index is exclusive!
//! assert!(matches == vec![Match{start: 0, end: 5, k: 1}]);
//! ```
//! Sometimes, it is necessary to use the slightly lower level, but also more powerful routines
//! that `triple_accel` provides. For example, it is possible to allow transpositions (character swaps)
//! that have a cost of 1, in addition to mismatches and gaps:
//! ```
//! use triple_accel::levenshtein::*;
//!
//! let a = b"abcd";
//! let b = b"abdc";
//! let k = 2; // upper bound on allowed cost
//! let trace_on = false; // return edit traceback?
//!
//! let dist = levenshtein_simd_k_with_opts(a, b, k, trace_on, RDAMERAU_COSTS);
//! // note: dist may be None if a and b do not match within a cost of k
//! assert!(dist.unwrap().0 == 1);
//! ```
//! Don't let the name of the function fool you! `levenshtein_simd_k_with_opts` will still fall back to
//! the scalar implementation if AVX2 or SSE4.1 support is not available. It just prefers to use SIMD
//! where possible.
//!
//! For most common cases, the re-exported functions are enough, and the low level functions do not
//! have to be used directly.

use std::*;

mod jewel;
pub mod hamming;
pub mod levenshtein;

// re-export common functions
pub use hamming::{hamming, hamming_search};
pub use levenshtein::{levenshtein, rdamerau, levenshtein_exp, rdamerau_exp, levenshtein_search};

// some shared utility stuff below

/// A struct that describes a single matching location.
///
/// This is usually returned as part of searching routines.
#[derive(Debug, PartialEq)]
pub struct Match {
    /// The start index of the match (inclusive).
    pub start: usize,
    /// The end index of the match (exclusive).
    pub end: usize,
    /// Number of edits for the match.
    pub k: u32
}

/// An enum describing possible edit operations.
///
/// This is usually returned as part of the traceback for edit distance routines.
#[derive(Debug, PartialEq)]
pub enum EditType {
    Match,
    Mismatch,
    AGap,
    BGap,
    Transpose
}

/// A struct representing a sequence of edits of the same type.
///
/// This is returned in the run-length encoded traceback of edit distance routines.
#[derive(Debug, PartialEq)]
pub struct Edit {
    /// The type of edit operation.
    pub edit: EditType,
    /// The number of consecutive edit operations of the same type.
    pub count: usize
}

/// An enum representing whether to return all matches or just the best matches.
///
/// This is used as an argument for searching routines.
#[derive(Debug, PartialEq, Copy, Clone)]
pub enum SearchType {
    All,
    Best
}

/// This creates a vector with the alignment and padding for `u128` values, and
/// then convert it to a vector of `u8` values that is returned.
///
/// This is possible because u8 has looser alignment requirements than `u128`.
/// This vector can be easily converted back to `u128` or `u64` later, for Hamming
/// distance routines.
/// The returned vector can be edited by copying `u8` values into it.
/// However, do not do any operation (like `push`) that may cause the the vector to be
/// reallocated.
///
/// # Arguments
/// * `len` - the length of the resulting array of u8 values
///
/// # Example
/// ```
/// # use triple_accel::*;
/// let s = alloc_str(10);
///
/// assert!(s.len() == 10);
/// ```
#[inline]
pub fn alloc_str(len: usize) -> Vec<u8> {
    let words_len = (len >> 4) + (if (len & 15) > 0 {1} else {0});
    let words = vec![0u128; words_len];
    let mut words = mem::ManuallyDrop::new(words);

    unsafe {
        Vec::from_raw_parts(words.as_mut_ptr() as *mut u8, len, words_len << 4)
    }
}

/// Directly copy from the a source `u8` slice to a destination `u8` slice.
///
/// Can be used to copy string data after allocating a vector using `alloc_str`.
///
/// # Arguments
/// * `dest` - the destination slice
/// * `src` - the source slice
///
/// # Panics
/// * If the length of `src` is greater than the length of `dest`.
///
/// # Example
/// ```
/// # use triple_accel::*;
/// let mut a = vec![0u8; 5];
/// let b = vec![1u8, 2u8, 3u8, 4u8];
///
/// fill_str(&mut a, &b);
///
/// assert!(a == vec![1u8, 2u8, 3u8, 4u8, 0u8]);
/// ```
#[inline]
pub fn fill_str(dest: &mut [u8], src: &[u8]) {
    assert!(dest.len() >= src.len());

    unsafe {
        ptr::copy_nonoverlapping(src.as_ptr(), dest.as_mut_ptr(), src.len());
    }
}

fn check_no_null_bytes(s: &[u8]) {
    for i in 0..s.len() {
        if s[i] == 0u8 {
            panic!("No zero/null bytes allowed in the string!");
        }
    }
}

