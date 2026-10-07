//! This module provides many Hamming distance routines.
//!
//! These distance functions share the same efficient underlying SIMD-accelerated implementation:
//! * `hamming`
//! * `hamming_simd_parallel`
//!
//! These search functions share the same efficient underlying SIMD-accelerated implementation:
//! * `hamming_search`
//! * `hamming_search_simd`
//! * `hamming_search_simd_with_opts`

use std::*;

use super::*;
use super::jewel::*;

/// Returns the hamming distance between two strings by naively counting mismatches.
///
/// The length of `a` and `b` must be the same.
///
/// # Arguments
/// * `a` - first string (slice)
/// * `b` - second string (slice)
///
/// # Panics
/// * If the length of `a` does not equal the length of `b`.
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::hamming::*;
/// let dist = hamming_naive(b"abc", b"abd");
///
/// assert!(dist == 1);
/// ```
pub fn hamming_naive(a: &[u8], b: &[u8]) -> u32 {
    let len = a.len();
    assert!(len == b.len());

    let mut res = 0u32;

    for i in 0..len {
        res += (a[i] != b[i]) as u32;
    }

    res
}

/// Returns an iterator over best `Match`s by naively searching through the text `haystack`
/// for the pattern `needle`.
///
/// This is done by naively counting mismatches at every position in `haystack`.
/// Only the matches with the lowest Hamming distance are returned.
/// Each returned `Match` requires at least half or more bytes of the `needle` to match
/// somewhere in the `haystack`.
/// The length of `needle` must be less than or equal to the length of `haystack`.
///
/// # Arguments
/// * `needle` - pattern string (slice)
/// * `haystack` - text string (slice)
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::hamming::*;
/// let matches: Vec<Match> = hamming_search_naive(b"abc", b"  abd").collect();
///
/// assert!(matches == vec![Match{start: 2, end: 5, k: 1}]);
/// ```
pub fn hamming_search_naive<'a>(needle: &'a [u8], haystack: &'a [u8]) -> Box<dyn Iterator<Item = Match> + 'a> {
    hamming_search_naive_with_opts(needle, haystack, ((needle.len() as u32) >> 1) + ((needle.len() as u32) & 1), SearchType::Best)
}

/// Returns an iterator over `Match`s by naively searching through the text `haystack`
/// for the pattern `needle`, with extra options.
///
/// Only matches with less than `k` mismatches are returned.
/// This is done by naively counting mismatches at every position in `haystack`.
/// The length of `needle` must be less than or equal to the length of `haystack`.
///
/// # Arguments
/// * `needle` - pattern string (slice)
/// * `haystack` - text string (slice)
/// * `k` - number of mismatches allowed
/// * `search_type` - whether to only return the "best" matches with the lowest Hamming distance, or
/// all matches
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::hamming::*;
/// let matches: Vec<Match> = hamming_search_naive_with_opts(b"abc", b"  abd", 1, SearchType::All).collect();
///
/// assert!(matches == vec![Match{start: 2, end: 5, k: 1}]);
/// ```
pub fn hamming_search_naive_with_opts<'a>(needle: &'a [u8], haystack: &'a [u8], k: u32, search_type: SearchType) -> Box<dyn Iterator<Item = Match> + 'a> {
    let needle_len = needle.len();
    let haystack_len = haystack.len();

    if needle_len > haystack_len {
        return Box::new(iter::empty());
    }

    let len = haystack_len + 1 - needle_len;
    let mut curr_k = k;
    let mut i = 0;

    let res = iter::from_fn(move || {
        'outer: while i < len {
            let mut final_res = 0u32;

            for j in 0..needle_len {
                final_res += (needle[j] != haystack[i + j]) as u32;

                // early stop
                if final_res > curr_k {
                    i += 1;
                    continue 'outer;
                }
            }

            match search_type {
                SearchType::Best => curr_k = final_res,
                _ => ()
            }

            i += 1;

            return Some((Match{start: i - 1, end: i + needle_len - 1, k: final_res}, curr_k));
        }

        None
    });

    if search_type == SearchType::Best {
        let mut res_vec = Vec::with_capacity(haystack_len / needle_len);
        res.for_each(|m| {
            res_vec.push(m.0);
            curr_k = m.1;
        });

        return Box::new(res_vec.into_iter().filter(move |m| m.k == curr_k));
    }

    Box::new(res.map(|m| m.0))
}

/// Returns the hamming distance between two strings by efficiently counting mismatches in chunks of 64 bits.
///
/// The length of `a` and `b` must be the same.
/// Both `a` and `b` must be aligned and padded so they can be directly casted to chunks of `u64`.
/// Use `alloc_str` to create aligned and padded strings.
/// This should be faster than `hamming_naive` and maybe even `hamming_words_128`. This should be slower
/// than `hamming_simd_parallel/movemask`.
///
/// # Arguments
/// * `a` - first string (slice)
/// * `b` - second string (slice)
///
/// # Panics
/// * If the length of `a` does not equal the length of `b`.
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::hamming::*;
/// let mut a = alloc_str(3);
/// let mut b = alloc_str(3);
/// fill_str(&mut a, b"abc");
/// fill_str(&mut b, b"abd");
///
/// let dist = hamming_words_64(&a, &b);
///
/// assert!(dist == 1);
/// ```
pub fn hamming_words_64(a: &[u8], b: &[u8]) -> u32 {
    assert!(a.len() == b.len());

    unsafe {
        let mut res = 0u32;
        // the pointer address better be aligned for u64
        // may not be in little endian
        let a_ptr = a.as_ptr() as *const u64;
        let b_ptr = b.as_ptr() as *const u64;
        let words_len = (a.len() >> 3) as isize;

        for i in 0..words_len {
            // change to little endian omitted because it is not necessary in this case
            let mut r = (*a_ptr.offset(i)) ^ (*b_ptr.offset(i));
            // reduce or by "folding" one half of each byte onto the other multiple times
            r |= r >> 4;
            // ...00001111
            r &= 0x0f0f0f0f0f0f0f0fu64;
            r |= r >> 2;
            // ...00110011
            r &= 0x3333333333333333u64;
            r |= r >> 1;
            // ...01010101
            r &= 0x5555555555555555u64;
            res += r.count_ones();
        }

        let words_rem = a.len() & 7;

        if words_rem > 0 {
            let mut r = (*a_ptr.offset(words_len)) ^ (*b_ptr.offset(words_len));
            r |= r >> 4;
            r &= 0x0f0f0f0f0f0f0f0fu64;
            r |= r >> 2;
            r &= 0x3333333333333333u64;
            r |= r >> 1;
            r &= 0x5555555555555555u64;
            // make sure to mask out bits outside the string lengths
            res += (r & ((1u64 << ((words_rem as u64) << 3u64)) - 1u64)).count_ones();
        }

        res
    }
}

/// Returns the hamming distance between two strings by counting mismatches in chunks of 128 bits.
///
/// The length of `a` and `b` must be the same.
/// Both `a` and `b` must be aligned and padded so they can be directly casted to chunks of `u128`.
/// Use `alloc_str` to create aligned and padded strings.
/// This may be slower than `hamming_words_64` in practice, probably since Rust `u128` is not as
/// optimized. This should be slower than `hamming_simd_parallel/movemask`.
///
/// # Arguments
/// * `a` - first string (slice)
/// * `b` - second string (slice)
///
/// # Panics
/// * If the length of `a` does not equal the length of `b`.
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::hamming::*;
/// let mut a = alloc_str(3);
/// let mut b = alloc_str(3);
/// fill_str(&mut a, b"abc");
/// fill_str(&mut b, b"abd");
///
/// let dist = hamming_words_128(&a, &b);
///
/// assert!(dist == 1);
/// ```
pub fn hamming_words_128(a: &[u8], b: &[u8]) -> u32 {
    assert!(a.len() == b.len());

    unsafe {
        let mut res = 0u32;
        // the pointer address better be aligned for u128
        // may not be in little endian
        let a_ptr = a.as_ptr() as *const u128;
        let b_ptr = b.as_ptr() as *const u128;
        let words_len = (a.len() >> 4) as isize;

        for i in 0..words_len {
            // change to little endian omitted because it is not necessary in this case
            let mut r = (*a_ptr.offset(i)) ^ (*b_ptr.offset(i));
            // reduce or by "folding" one half of each byte onto the other multiple times
            r |= r >> 4;
            // ...00001111
            r &= 0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0fu128;
            r |= r >> 2;
            // ...00110011
            r &= 0x33333333333333333333333333333333u128;
            r |= r >> 1;
            // ...01010101
            r &= 0x55555555555555555555555555555555u128;
            res += r.count_ones();
        }

        let words_rem = a.len() & 15;

        if words_rem > 0 {
            let mut r = (*a_ptr.offset(words_len)) ^ (*b_ptr.offset(words_len));
            r |= r >> 4;
            r &= 0x0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f0fu128;
            r |= r >> 2;
            r &= 0x33333333333333333333333333333333u128;
            r |= r >> 1;
            r &= 0x55555555555555555555555555555555u128;
            // make sure to mask out bits outside the string lengths
            res += (r & ((1u128 << ((words_rem as u128) << 3u128)) - 1u128)).count_ones();
        }

        res
    }
}

/// Returns the hamming distance between two strings by counting mismatches using SIMD vectors to
/// increment multiple counters in parallel.
///
/// The length of `a` and `b` must be the same.
/// There are no constraints on how `a` and `b` are aligned and padded.
/// This will automatically fall back to `hamming_naive`, if AVX2 and SSE4.1 are not supported.
/// This should be faster than both `hamming_word_64/128` and `hamming_simd_movemask`.
///
/// # Arguments
/// * `a` - first string (slice)
/// * `b` - second string (slice)
///
/// # Panics
/// * If the length of `a` does not equal the length of `b`.
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::hamming::*;
/// let dist = hamming_simd_parallel(b"abc", b"abd");
///
/// assert!(dist == 1);
/// ```
pub fn hamming_simd_parallel(a: &[u8], b: &[u8]) -> u32 {
    assert!(a.len() == b.len());

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if cfg!(feature = "jewel-avx") && is_x86_feature_detected!("avx2") {
            return unsafe {Avx::count_mismatches(a.as_ptr(), b.as_ptr(), a.len())};
        }else if cfg!(feature = "jewel-sse") && is_x86_feature_detected!("sse4.1") {
            return unsafe {Sse::count_mismatches(a.as_ptr(), b.as_ptr(), a.len())};
        }
    }

    hamming_naive(a, b)
}

/// Returns the hamming distance between two strings by counting mismatches using the SIMD movemask intrinsic.
///
/// The length of `a` and `b` must be the same.
/// There are no constraints on how `a` and `b` are aligned and padded.
/// This will automatically fall back to `hamming_naive`, if AVX2 and SSE4.1 are not supported.
/// This should be faster than `hamming_word_64/128`, but slower than `hamming_simd_parallel`.
///
/// # Arguments
/// * `a` - first string (slice)
/// * `b` - second string (slice)
///
/// # Panics
/// * If the length of `a` does not equal the length of `b`.
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::hamming::*;
/// let dist = hamming_simd_movemask(b"abc", b"abd");
///
/// assert!(dist == 1);
/// ```
pub fn hamming_simd_movemask(a: &[u8], b: &[u8]) -> u32 {
    assert!(a.len() == b.len());

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if cfg!(feature = "jewel-avx") && is_x86_feature_detected!("avx2") {
            return unsafe {Avx::mm_count_mismatches(a.as_ptr(), b.as_ptr(), a.len())};
        }else if cfg!(feature = "jewel-sse") && is_x86_feature_detected!("sse4.1") {
            return unsafe {Sse::mm_count_mismatches(a.as_ptr(), b.as_ptr(), a.len())};
        }
    }

    hamming_naive(a, b)
}

/// Returns the hamming distance between two strings using the best method.
///
/// The length of `a` and `b` must be the same.
/// This will automatically fall back to a scalar alternative if AVX2 and
/// SSE4.1 are not supported.
/// Internally, this calls `hamming_simd_parallel`.
///
/// # Arguments
/// * `a` - first string (slice)
/// * `b` - second string (slice)
///
/// # Panics
/// * If the length of `a` does not equal the length of `b`.
///
/// # Example
/// ```
/// # use triple_accel::*;
/// let dist = hamming(b"abc", b"abd");
///
/// assert!(dist == 1);
/// ```
pub fn hamming(a: &[u8], b: &[u8]) -> u32 {
    hamming_simd_parallel(a, b)
}

/// Returns an iterator over best `Match`s by searching through the text `haystack`
/// for the pattern `needle` using SIMD.
///
/// This is done by counting mismatches at every position in `haystack`.
/// This will automatically fall back to `hamming_search_naive_with_opts` if AVX2 and SSE4.1
/// are not supported.
/// Null bytes/characters are not supported.
/// The length of `needle` must be less than or equal to the length of `haystack`.
/// Each returned `Match` requires at least half or more bytes of the `needle` to match
/// somwhere in the `haystack`.
/// Only the matches with the lowest Hamming distance are returned.
/// This should be faster than `hamming_search_naive`.
///
/// # Arguments
/// * `needle` - pattern string (slice)
/// * `haystack` - text string (slice)
///
/// # Panics
/// * When there are zero/null bytes in the `haystack` string.
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::hamming::*;
/// let matches: Vec<Match> = hamming_search_simd(b"abc", b"  abd").collect();
///
/// assert!(matches == vec![Match{start: 2, end: 5, k: 1}]);
/// ```
pub fn hamming_search_simd<'a>(needle: &'a [u8], haystack: &'a [u8]) -> Box<dyn Iterator<Item = Match> + 'a> {
    hamming_search_simd_with_opts(needle, haystack, ((needle.len() as u32) >> 1) + ((needle.len() as u32) & 1), SearchType::Best)
}

/// Returns an iterator over `Match`s by searching through the text `haystack` for the
/// pattern `needle` using SIMD, with extra options.
///
/// This is done by using SIMD to count mismatches at every position in `haystack`.
/// This will automatically fall back to `hamming_search_naive_with_opts` if AVX2 and SSE4.1
/// are not supported.
/// Null bytes/characters are not supported.
/// The length of `needle` must be less than or equal to the length of `haystack`.
/// This should be faster than `hamming_search_naive_with_opts`.
///
/// # Arguments
/// * `needle` - pattern string (slice)
/// * `haystack` - text string (slice)
/// * `k` - number of mismatches allowed
/// * `search_type` - whether to only return the "best" matches with the lowest Hamming distance, or
/// all matches
///
/// # Panics
/// * When there are zero/null bytes in the `haystack` string.
///
/// # Example
/// ```
/// # use triple_accel::*;
/// # use triple_accel::hamming::*;
/// let matches: Vec<Match> = hamming_search_simd_with_opts(b"abc", b"  abd", 1, SearchType::All).collect();
///
/// assert!(matches == vec![Match{start: 2, end: 5, k: 1}]);
/// ```
pub fn hamming_search_simd_with_opts<'a>(needle: &'a [u8], haystack: &'a [u8], k: u32, search_type: SearchType) -> Box<dyn Iterator<Item = Match> + 'a> {
    if needle.len() > haystack.len() {
        return Box::new(iter::empty());
    }

    if needle.len() == 0 {
        return Box::new(iter::empty());
    }

    check_no_null_bytes(haystack);

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if cfg!(feature = "jewel-avx") && is_x86_feature_detected!("avx2") {
            return unsafe {hamming_search_simd_core_avx(needle, haystack, k, search_type)};
        }else if cfg!(feature = "jewel-sse") && is_x86_feature_detected!("sse4.1") {
            return unsafe {hamming_search_simd_core_sse(needle, haystack, k, search_type)};
        }
    }

    hamming_search_naive_with_opts(needle, haystack, k, search_type)
}

macro_rules! create_hamming_search_simd_core {
    ($name:ident, $jewel:ty, $target:literal) => {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        #[target_feature(enable = $target)]
        unsafe fn $name<'a>(needle: &'a [u8], haystack: &'a [u8], k: u32, search_type: SearchType) -> Box<dyn Iterator<Item = Match> + 'a> {
            #[cfg(feature = "debug")]
            {
                println!("Debug: Hamming search Jewel vector type {} for target {}.", stringify!($jewel), stringify!($target));
            }

            let needle_len = needle.len();
            let haystack_len = haystack.len();
            let needle_vector = <$jewel>::loadu(needle.as_ptr(), needle_len);
            // calculate len using the unused bytes in the needle Jewel vector, for speed
            // there may be leftover positions in haystack that need to be calculated using a
            // scalar search afterwards
            // there should be no null bytes in the strings
            let len = if needle_vector.upper_bound() > haystack_len {0} else {haystack_len + 1 - needle_vector.upper_bound()};
            let real_len = haystack_len + 1 - needle_len;
            let haystack_ptr = haystack.as_ptr();
            let mut curr_k = k;
            let mut i = 0;

            let res = iter::from_fn(move || {
                while i < len {
                    let final_res = <$jewel>::vector_count_mismatches(&needle_vector, haystack_ptr.offset(i as isize), needle_len);
                    i += 1;

                    if final_res <= curr_k {
                        match search_type {
                            SearchType::Best => curr_k = final_res,
                            _ => ()
                        }

                        return Some((Match{start: i - 1, end: i + needle_len - 1, k: final_res}, curr_k));
                    }
                }

                // scalar search
                'outer: while i < real_len {
                    let mut final_res = 0u32;

                    for j in 0..needle_len {
                        final_res += (*needle.get_unchecked(j) != *haystack.get_unchecked(i + j)) as u32;

                        if final_res > curr_k {
                            i += 1;
                            continue 'outer;
                        }
                    }

                    match search_type {
                        SearchType::Best => curr_k = final_res,
                        _ => ()
                    }

                    i += 1;

                    return Some((Match{start: i - 1, end: i + needle_len - 1, k: final_res}, curr_k));
                }

                None
            });

            if search_type == SearchType::Best {
                let mut res_vec = Vec::with_capacity(haystack_len / needle_len);
                res.for_each(|m| {
                    res_vec.push(m.0);
                    curr_k = m.1;
                });

                return Box::new(res_vec.into_iter().filter(move |m| m.k == curr_k));
            }

            Box::new(res.map(|m| m.0))
        }
    };
}

// generate different versions for different intrinsics
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_hamming_search_simd_core!(hamming_search_simd_core_avx, Avx, "avx2");
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
create_hamming_search_simd_core!(hamming_search_simd_core_sse, Sse, "sse4.1");

/// Returns an iterator over best `Match`s by searching through the text `haystack`
/// for the pattern `needle` using SIMD.
///
/// This will automatically fall back to a scalar alternative if AVX2 and SSE4.1
/// are not supported.
/// Null bytes/characters are not supported.
/// The length of `needle` must be less than or equal to the length of `haystack`.
/// Each returned `Match` requires at least half or more bytes of the `needle` to match
/// somewhere in the `haystack`.
/// Only the matches with the lowest Hamming distance are returned.
/// Internally, this calls `hamming_search_simd`.
///
/// # Arguments
/// * `needle` - pattern string (slice)
/// * `haystack` - text string (slice)
///
/// # Panics
/// * When there are zero/null bytes in the `haystack` string.
///
/// # Example
/// ```
/// # use triple_accel::*;
/// let matches: Vec<Match> = hamming_search(b"abc", b"  abd").collect();
///
/// assert!(matches == vec![Match{start: 2, end: 5, k: 1}]);
/// ```
pub fn hamming_search<'a>(needle: &'a [u8], haystack: &'a [u8]) -> Box<dyn Iterator<Item = Match> + 'a> {
    hamming_search_simd(needle, haystack)
}

