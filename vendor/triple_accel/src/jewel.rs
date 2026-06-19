//! This module provides wrappers for SIMD intrinsics, so they can be used on multiple platforms.

use std::*;

#[cfg(target_arch = "x86")]
use std::arch::x86::*;

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

/// Jewel provides a uniform interface for SIMD operations.
///
/// To save space, most operations are modify in place.
/// Jewel vectors can be easily printed for debugging purposes.
/// Additionally, the functions should be inlined to the caller for
/// maximum efficiency.
pub trait Jewel: fmt::Display {
    /// Functions for allocating memory and creating a new Jewel vector.
    unsafe fn repeating(val: u32, len: usize) -> Self;
    unsafe fn repeating_max(len: usize) -> Self;

    /// Figure out the length of the created vector, which may
    /// be longer than the length given by the caller.
    fn upper_bound(&self) -> usize;
    /// Figure out the length if it is static.
    fn static_upper_bound() -> usize;

    /// These operations do not have to be very efficient.
    unsafe fn slow_loadu(&mut self, idx: usize, ptr: *const u8, len: usize, reverse: bool);
    unsafe fn slow_extract(&self, i: usize) -> u32;
    unsafe fn slow_insert(&mut self, i: usize, val: u32);

    /// These operations modify the Jewel struct where the function is called.
    /// last_0 is the last element, last_1 is the second to last, etc.
    unsafe fn insert_last_0(&mut self, val: u32);
    unsafe fn insert_last_1(&mut self, val: u32);
    unsafe fn insert_last_2(&mut self, val: u32);
    unsafe fn insert_last_max(&mut self);
    unsafe fn insert_first(&mut self, val: u32);
    unsafe fn insert_first_max(&mut self);

    unsafe fn add_mut(&mut self, b: &Self);
    unsafe fn adds_mut(&mut self, b: &Self);
    unsafe fn and_mut(&mut self, b: &Self);
    unsafe fn andnot_mut(&mut self, b: &Self);
    unsafe fn cmpeq_mut(&mut self, b: &Self);
    unsafe fn min_mut(&mut self, b: &Self);
    unsafe fn max_mut(&mut self, b: &Self);
    unsafe fn blendv_mut(&mut self, b: &Self, mask: &Self);
    unsafe fn shift_left_1_mut(&mut self);
    unsafe fn shift_left_2_mut(&mut self);
    unsafe fn shift_right_1_mut(&mut self);

    /// These operations overwrite a res vector to reduce memory allocations.
    unsafe fn add(a: &Self, b: &Self, res: &mut Self);
    unsafe fn adds(a: &Self, b: &Self, res: &mut Self);
    unsafe fn andnot(a: &Self, b: &Self, res: &mut Self);
    unsafe fn cmpeq(a: &Self, b: &Self, res: &mut Self);
    unsafe fn min(a: &Self, b: &Self, res: &mut Self);
    unsafe fn max(a: &Self, b: &Self, res: &mut Self);
    unsafe fn shift_left_1(a: &Self, res: &mut Self);
    unsafe fn shift_right_1(a: &Self, res: &mut Self);

    /// `triple_argmin` will allocate memory and create a new Jewel vector.
    unsafe fn triple_argmin(sub: &Self, a_gap: &Self, b_gap: &Self, res_min: &mut Self) -> Self;
    unsafe fn triple_min_length(sub: &Self, a_gap: &Self, b_gap: &Self, sub_length: &Self,
                                a_gap_length: &Self, b_gap_length: &Self, res_min: &mut Self, res_length: &mut Self);
    unsafe fn double_min_length(new_gap: &Self, res_cont_gap: &mut Self, new_gap_length: &Self, res_cont_gap_length: &mut Self);
}

// macros to help generate implementations for some of the Jewel vector functions
macro_rules! operation_param2 {
    ($target:literal, $fn_name:ident, $intrinsic:ident) => {
        #[target_feature(enable = $target)]
        #[inline]
        unsafe fn $fn_name(a: &Self, b: &Self, res: &mut Self) {
            for i in 0..a.v.len() {
                *res.v.get_unchecked_mut(i) = $intrinsic(*a.v.get_unchecked(i), *b.v.get_unchecked(i));
            }
        }
    };
}

macro_rules! operation_mut_param2 {
    ($target:literal, $fn_name:ident, $intrinsic:ident) => {
        #[target_feature(enable = $target)]
        #[inline]
        unsafe fn $fn_name(&mut self, b: &Self) {
            for i in 0..self.v.len() {
                *self.v.get_unchecked_mut(i) = $intrinsic(*self.v.get_unchecked(i), *b.v.get_unchecked(i));
            }
        }
    };
}

/// N x 32 x 8 vector backed with 256-bit AVX vectors.
macro_rules! create_avx_nx32x8 {
    ($name:ident, $num:literal) => {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        pub struct $name {
            v: [__m256i; $num]
        }

        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        impl Jewel for $name {
            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn repeating(val: u32, _len: usize) -> Self {
                let v = [_mm256_set1_epi8(val as i8); $num];

                Self{
                    v: v
                }
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn repeating_max(_len: usize) -> Self {
                let v = [_mm256_set1_epi8(-1i8); $num];

                Self{
                    v: v
                }
            }

            #[inline]
            fn upper_bound(&self) -> usize {
                self.v.len() << 5
            }

            #[inline]
            fn static_upper_bound() -> usize {
                $num << 5
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn slow_loadu(&mut self, idx: usize, ptr: *const u8, len: usize, reverse: bool) {
                if len == 0 {
                    return;
                }

                let mut arr = [0u8; 32];
                let arr_ptr = arr.as_mut_ptr() as *mut __m256i;
                let store_idx = if reverse {31} else {0};
                let load_idx = if reverse {0} else {31};

                for i in 0..len {
                    let curr_idx = if reverse {idx - i} else {idx + i};
                    let arr_idx = curr_idx & 31;

                    if arr_idx == store_idx || i == 0 {
                        // use part of the original array
                        _mm256_storeu_si256(arr_ptr, *self.v.get_unchecked(curr_idx >> 5));
                    }

                    *arr.get_unchecked_mut(arr_idx) = *ptr.offset(i as isize);

                    if arr_idx == load_idx || i == len - 1 {
                        *self.v.get_unchecked_mut(curr_idx >> 5) = _mm256_loadu_si256(arr_ptr);
                    }
                }
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn slow_extract(&self, i: usize) -> u32 {
                let idx = i >> 5;
                let j = i & 31;
                let mut arr = [0u8; 32];
                _mm256_storeu_si256(arr.as_mut_ptr() as *mut __m256i, *self.v.get_unchecked(idx));
                *arr.get_unchecked(j) as u32
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn slow_insert(&mut self, i: usize, val: u32) {
                let idx = i >> 5;
                let j = i & 31;
                let mut arr = [0u8; 32];
                let arr_ptr = arr.as_mut_ptr() as *mut __m256i;
                _mm256_storeu_si256(arr_ptr, *self.v.get_unchecked(idx));
                *arr.get_unchecked_mut(j) = val as u8;
                *self.v.get_unchecked_mut(idx) = _mm256_loadu_si256(arr_ptr);
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn insert_last_0(&mut self, val: u32) {
                let last = self.v.len() - 1;
                *self.v.get_unchecked_mut(last) = _mm256_insert_epi8(*self.v.get_unchecked(last), val as i8, 31i32);
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn insert_last_1(&mut self, val: u32) {
                let last = self.v.len() - 1;
                *self.v.get_unchecked_mut(last) = _mm256_insert_epi8(*self.v.get_unchecked(last), val as i8, 30i32);
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn insert_last_2(&mut self, val: u32) {
                let last = self.v.len() - 1;
                *self.v.get_unchecked_mut(last) = _mm256_insert_epi8(*self.v.get_unchecked(last), val as i8, 29i32);
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn insert_last_max(&mut self) {
                let last = self.v.len() - 1;
                *self.v.get_unchecked_mut(last) = _mm256_insert_epi8(*self.v.get_unchecked(last), -1i8, 31i32);
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn insert_first(&mut self, val: u32) {
                *self.v.get_unchecked_mut(0) = _mm256_insert_epi8(*self.v.get_unchecked(0), val as i8, 0i32);
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn insert_first_max(&mut self) {
                *self.v.get_unchecked_mut(0) = _mm256_insert_epi8(*self.v.get_unchecked(0), -1i8, 0i32);
            }

            operation_mut_param2!("avx2", add_mut, _mm256_add_epi8);
            operation_mut_param2!("avx2", adds_mut, _mm256_adds_epu8);
            operation_mut_param2!("avx2", and_mut, _mm256_and_si256);
            operation_mut_param2!("avx2", andnot_mut, _mm256_andnot_si256);
            operation_mut_param2!("avx2", cmpeq_mut, _mm256_cmpeq_epi8);
            operation_mut_param2!("avx2", min_mut, _mm256_min_epu8);
            operation_mut_param2!("avx2", max_mut, _mm256_max_epu8);

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn blendv_mut(&mut self, b: &Self, mask: &Self) {
                for i in 0..self.v.len() {
                    *self.v.get_unchecked_mut(i) = _mm256_blendv_epi8(*self.v.get_unchecked(i), *b.v.get_unchecked(i), *mask.v.get_unchecked(i));
                }
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn shift_left_1_mut(&mut self) {
                for i in 0..(self.v.len() - 1) {
                    let curr = *self.v.get_unchecked(i);
                    // permute concatenates the second half of the current vector and the first half of the next vector
                    *self.v.get_unchecked_mut(i) = _mm256_alignr_epi8(
                        _mm256_permute2x128_si256(curr, *self.v.get_unchecked(i + 1), 0b00100001i32), curr, 1i32);
                }

                // last one gets to shift in zeros
                let last = self.v.len() - 1;
                let curr = *self.v.get_unchecked(last);
                // permute concatenates the second half of the last vector and a vector of zeros
                *self.v.get_unchecked_mut(last) = _mm256_alignr_epi8(_mm256_permute2x128_si256(curr, curr, 0b10000001i32), curr, 1i32);
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn shift_left_2_mut(&mut self) {
                for i in 0..(self.v.len() - 1) {
                    let curr = *self.v.get_unchecked(i);
                    // permute concatenates the second half of the current vector and the first half of the next vector
                    *self.v.get_unchecked_mut(i) = _mm256_alignr_epi8(
                        _mm256_permute2x128_si256(curr, *self.v.get_unchecked(i + 1), 0b00100001i32), curr, 2i32);
                }

                // last one gets to shift in zeros
                let last = self.v.len() - 1;
                let curr = *self.v.get_unchecked(last);
                // permute concatenates the second half of the last vector and a vector of zeros
                *self.v.get_unchecked_mut(last) = _mm256_alignr_epi8(_mm256_permute2x128_si256(curr, curr, 0b10000001i32), curr, 2i32);
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn shift_right_1_mut(&mut self) {
                for i in (1..self.v.len()).rev() {
                    let curr = *self.v.get_unchecked(i);
                    // permute concatenates the second half of the previous vector and the first half of the current vector
                    *self.v.get_unchecked_mut(i) = _mm256_alignr_epi8(
                        curr, _mm256_permute2x128_si256(curr, *self.v.get_unchecked(i - 1), 0b00000011i32), 15i32);
                }

                // first one gets to shift in zeros
                let curr = *self.v.get_unchecked(0);
                // permute concatenates a vector of zeros and the first half of the first vector
                *self.v.get_unchecked_mut(0) = _mm256_alignr_epi8(curr, _mm256_permute2x128_si256(curr, curr, 0b00001000i32), 15i32);
            }

            operation_param2!("avx2", add, _mm256_add_epi8);
            operation_param2!("avx2", adds, _mm256_adds_epu8);
            operation_param2!("avx2", andnot, _mm256_andnot_si256);
            operation_param2!("avx2", cmpeq, _mm256_cmpeq_epi8);
            operation_param2!("avx2", min, _mm256_min_epu8);
            operation_param2!("avx2", max, _mm256_max_epu8);

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn shift_left_1(a: &Self, res: &mut Self) {
                for i in 0..(a.v.len() - 1) {
                    let curr = *a.v.get_unchecked(i);
                    // permute concatenates the second half of the current vector and the first half of the next vector
                    *res.v.get_unchecked_mut(i) = _mm256_alignr_epi8(
                        _mm256_permute2x128_si256(curr, *a.v.get_unchecked(i + 1), 0b00100001i32), curr, 1i32);
                }

                // last one gets to shift in zeros
                let last = a.v.len() - 1;
                let curr = *a.v.get_unchecked(last);
                // permute concatenates the second half of the last vector and a vector of zeros
                *res.v.get_unchecked_mut(last) = _mm256_alignr_epi8(_mm256_permute2x128_si256(curr, curr, 0b10000001i32), curr, 1i32);
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn shift_right_1(a: &Self, res: &mut Self) {
                for i in (1..a.v.len()).rev() {
                    let curr = *a.v.get_unchecked(i);
                    // permute concatenates the second half of the previous vector and the first half of the current vector
                    *res.v.get_unchecked_mut(i) = _mm256_alignr_epi8(
                        curr, _mm256_permute2x128_si256(curr, *a.v.get_unchecked(i - 1), 0b00000011i32), 15i32);
                }

                // first one gets to shift in zeros
                let curr = *a.v.get_unchecked(0);
                // permute concatenates a vector of zeros and the first half of the first vector
                *res.v.get_unchecked_mut(0) = _mm256_alignr_epi8(curr, _mm256_permute2x128_si256(curr, curr, 0b00001000i32), 15i32);
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn triple_argmin(sub: &Self, a_gap: &Self, b_gap: &Self, res_min: &mut Self) -> Self {
                // return the edit used in addition to doing a min operation
                let mut v = [_mm256_undefined_si256(); $num];
                let twos = _mm256_set1_epi8(2);

                for i in 0..sub.v.len() {
                    let sub = *sub.v.get_unchecked(i);
                    let a_gap = *a_gap.v.get_unchecked(i);
                    let b_gap = *b_gap.v.get_unchecked(i);

                    let res_min1 = _mm256_min_epu8(a_gap, b_gap);
                    // a gap: 2 + -1 = 1, b gap: 2 + 0 = 2
                    let res_arg1 = _mm256_add_epi8(twos, _mm256_cmpeq_epi8(a_gap, res_min1));

                    let res_min2 = _mm256_min_epu8(sub, res_min1);
                    // sub: 0
                    let res_arg2 = _mm256_andnot_si256(_mm256_cmpeq_epi8(sub, res_min2), res_arg1);

                    *res_min.v.get_unchecked_mut(i) = res_min2;
                    *v.get_unchecked_mut(i) = res_arg2;
                }

                Self{
                    v: v
                }
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn triple_min_length(sub: &Self, a_gap: &Self,
                                        b_gap: &Self, sub_length: &Self, a_gap_length: &Self,
                                        b_gap_length: &Self, res_min: &mut Self, res_length: &mut Self) {
                // choose the length based on which edit is chosen during the min operation
                // secondary objective of maximizing length if edit costs equal
                for i in 0..sub.v.len() {
                    let sub = *sub.v.get_unchecked(i);
                    let a_gap = *a_gap.v.get_unchecked(i);
                    let b_gap = *b_gap.v.get_unchecked(i);
                    let sub_length = *sub_length.v.get_unchecked(i);
                    let a_gap_length = *a_gap_length.v.get_unchecked(i);
                    let b_gap_length = *b_gap_length.v.get_unchecked(i);

                    let res_min1 = _mm256_min_epu8(a_gap, b_gap);
                    let a_b_gt_mask = _mm256_cmpeq_epi8(a_gap, res_min1); // a gap: -1, b gap: 0
                    let mut res_length1 = _mm256_blendv_epi8(b_gap_length, a_gap_length, a_b_gt_mask); // lengths based on edits
                    let a_b_eq_mask = _mm256_cmpeq_epi8(a_gap, b_gap); // equal: -1
                    let a_b_max_len = _mm256_max_epu8(a_gap_length, b_gap_length);
                    res_length1 = _mm256_blendv_epi8(res_length1, a_b_max_len, a_b_eq_mask); // maximize length if edits equal

                    let res_min2 = _mm256_min_epu8(sub, res_min1);
                    let sub_gt_mask = _mm256_cmpeq_epi8(sub, res_min2); // sub: -1, prev a or b gap: 0
                    let mut res_length2 = _mm256_blendv_epi8(res_length1, sub_length, sub_gt_mask); // length based on edits
                    let sub_eq_mask = _mm256_cmpeq_epi8(sub, res_min1);
                    let sub_max_len = _mm256_max_epu8(sub_length, res_length1);
                    res_length2 = _mm256_blendv_epi8(res_length2, sub_max_len, sub_eq_mask); // maximize length if edits equal

                    *res_min.v.get_unchecked_mut(i) = res_min2;
                    *res_length.v.get_unchecked_mut(i) = res_length2;
                }
            }

            #[target_feature(enable = "avx2")]
            #[inline]
            unsafe fn double_min_length(new_gap: &Self, res_cont_gap: &mut Self, new_gap_length: &Self, res_cont_gap_length: &mut Self) {
                // choose the length based on which gap type is chosen during the min operation
                // secondary objective of maximizing length if edit costs equal
                for i in 0..new_gap.v.len() {
                    let new_gap = *new_gap.v.get_unchecked(i);
                    let cont_gap = *res_cont_gap.v.get_unchecked(i);
                    let new_gap_length = *new_gap_length.v.get_unchecked(i);
                    let cont_gap_length = *res_cont_gap_length.v.get_unchecked(i);

                    let res_min = _mm256_min_epu8(new_gap, cont_gap);
                    let new_cont_gt_mask = _mm256_cmpeq_epi8(new_gap, res_min); // new gap: -1, continue gap: 0
                    let mut res_length = _mm256_blendv_epi8(cont_gap_length, new_gap_length, new_cont_gt_mask); // lengths based on edits
                    let new_cont_eq_mask = _mm256_cmpeq_epi8(new_gap, cont_gap); // equal: -1
                    let new_cont_max_len = _mm256_max_epu8(new_gap_length, cont_gap_length);
                    res_length = _mm256_blendv_epi8(res_length, new_cont_max_len, new_cont_eq_mask); // maximize length if edits equal

                    *res_cont_gap.v.get_unchecked_mut(i) = res_min;
                    *res_cont_gap_length.v.get_unchecked_mut(i) = res_length;
                }
            }
        }

        // this implementation will probably only be used for debugging
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
                #[target_feature(enable = "avx2")]
                #[inline]
                unsafe fn fmt_internal(s: &$name, f: &mut fmt::Formatter) -> fmt::Result {
                    write!(f, "[")?;

                    let mut arr = [0u8; 32];
                    let arr_ptr = arr.as_mut_ptr() as *mut __m256i;

                    for i in 0..(s.v.len() - 1) {
                        _mm256_storeu_si256(arr_ptr, *s.v.get_unchecked(i));

                        for j in 0..32 {
                            write!(f, "{:>3}, ", *arr.get_unchecked(j))?;
                        }
                    }

                    // leftover elements

                    _mm256_storeu_si256(arr_ptr, *s.v.get_unchecked(s.v.len() - 1));

                    let start = (s.v.len() - 1) << 5;

                    for i in 0..(s.upper_bound() - start) {
                        if i == s.upper_bound() - start - 1 {
                            write!(f, "{:>3}", *arr.get_unchecked(i))?;
                        }else{
                            write!(f, "{:>3}, ", *arr.get_unchecked(i))?;
                        }
                    }

                    write!(f, "]")
                }

                unsafe { fmt_internal(self, f) }
            }
        }
    };
}

// constant array size, so the compiler should unroll the loops
create_avx_nx32x8!(Avx1x32x8, 1);
create_avx_nx32x8!(Avx2x32x8, 2);
create_avx_nx32x8!(Avx4x32x8, 4);
create_avx_nx32x8!(Avx8x32x8, 8);

/// N x 16 x 16 vector backed with 256-bit AVX vectors.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
pub struct AvxNx16x16 {
    v: Vec<__m256i>
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
impl Jewel for AvxNx16x16 {
    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn repeating(val: u32, len: usize) -> Self {
        let v = vec![_mm256_set1_epi16(val as i16); (len >> 4) + if (len & 15) > 0 {1} else {0}];

        Self{
            v: v
        }
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn repeating_max(len: usize) -> Self {
        let v = vec![_mm256_set1_epi16(-1i16); (len >> 4) + if (len & 15) > 0 {1} else {0}];

        Self{
            v: v
        }
    }

    #[inline]
    fn upper_bound(&self) -> usize {
        self.v.len() << 4
    }

    #[inline]
    fn static_upper_bound() -> usize {
        unimplemented!()
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn slow_loadu(&mut self, idx: usize, ptr: *const u8, len: usize, reverse: bool) {
        if len == 0 {
            return;
        }

        let mut arr = [0u16; 16];
        let arr_ptr = arr.as_mut_ptr() as *mut __m256i;
        let store_idx = if reverse {15} else {0};
        let load_idx = if reverse {0} else {15};

        for i in 0..len {
            let curr_idx = if reverse {idx - i} else {idx + i};
            let arr_idx = curr_idx & 15;

            if arr_idx == store_idx || i == 0 {
                _mm256_storeu_si256(arr_ptr, *self.v.get_unchecked(curr_idx >> 4));
            }

            *arr.get_unchecked_mut(arr_idx) = *ptr.offset(i as isize) as u16;

            if arr_idx == load_idx || i == len - 1 {
                *self.v.get_unchecked_mut(curr_idx >> 4) = _mm256_loadu_si256(arr_ptr);
            }
        }
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn slow_extract(&self, i: usize) -> u32 {
        let idx = i >> 4;
        let j = i & 15;
        let mut arr = [0u16; 16];
        _mm256_storeu_si256(arr.as_mut_ptr() as *mut __m256i, *self.v.get_unchecked(idx));
        *arr.get_unchecked(j) as u32
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn slow_insert(&mut self, i: usize, val: u32) {
        let idx = i >> 4;
        let j = i & 15;
        let mut arr = [0u16; 16];
        let arr_ptr = arr.as_mut_ptr() as *mut __m256i;
        _mm256_storeu_si256(arr_ptr, *self.v.get_unchecked(idx));
        *arr.get_unchecked_mut(j) = val as u16;
        *self.v.get_unchecked_mut(idx) = _mm256_loadu_si256(arr_ptr);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn insert_last_0(&mut self, val: u32) {
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm256_insert_epi16(*self.v.get_unchecked(last), val as i16, 15i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn insert_last_1(&mut self, val: u32) {
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm256_insert_epi16(*self.v.get_unchecked(last), val as i16, 14i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn insert_last_2(&mut self, val: u32) {
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm256_insert_epi16(*self.v.get_unchecked(last), val as i16, 13i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn insert_last_max(&mut self) {
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm256_insert_epi16(*self.v.get_unchecked(last), -1i16, 15i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn insert_first(&mut self, val: u32) {
        *self.v.get_unchecked_mut(0) = _mm256_insert_epi16(*self.v.get_unchecked(0), val as i16, 0i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn insert_first_max(&mut self) {
        *self.v.get_unchecked_mut(0) = _mm256_insert_epi16(*self.v.get_unchecked(0), -1i16, 0i32);
    }

    operation_mut_param2!("avx2", add_mut, _mm256_add_epi16);
    operation_mut_param2!("avx2", adds_mut, _mm256_adds_epu16);
    operation_mut_param2!("avx2", and_mut, _mm256_and_si256);
    operation_mut_param2!("avx2", andnot_mut, _mm256_andnot_si256);
    operation_mut_param2!("avx2", cmpeq_mut, _mm256_cmpeq_epi16);
    operation_mut_param2!("avx2", min_mut, _mm256_min_epu16);
    operation_mut_param2!("avx2", max_mut, _mm256_max_epu16);

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn blendv_mut(&mut self, b: &Self, mask: &Self) {
        for i in 0..self.v.len() {
            *self.v.get_unchecked_mut(i) = _mm256_blendv_epi8(*self.v.get_unchecked(i), *b.v.get_unchecked(i), *mask.v.get_unchecked(i));
        }
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn shift_left_1_mut(&mut self) {
        for i in 0..(self.v.len() - 1) {
            let curr = *self.v.get_unchecked(i);
            // permute concatenates the second half of the current vector and the first half of the next vector
            *self.v.get_unchecked_mut(i) = _mm256_alignr_epi8(
                _mm256_permute2x128_si256(curr, *self.v.get_unchecked(i + 1), 0b00100001i32), curr, 2i32);
        }

        // last one gets to shift in zeros
        let last = self.v.len() - 1;
        let curr = *self.v.get_unchecked(last);
        // permute concatenates the second half of the last vector and a vector of zeros
        *self.v.get_unchecked_mut(last) = _mm256_alignr_epi8(_mm256_permute2x128_si256(curr, curr, 0b10000001i32), curr, 2i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn shift_left_2_mut(&mut self) {
        for i in 0..(self.v.len() - 1) {
            let curr = *self.v.get_unchecked(i);
            // permute concatenates the second half of the current vector and the first half of the next vector
            *self.v.get_unchecked_mut(i) = _mm256_alignr_epi8(
                _mm256_permute2x128_si256(curr, *self.v.get_unchecked(i + 1), 0b00100001i32), curr, 4i32);
        }

        // last one gets to shift in zeros
        let last = self.v.len() - 1;
        let curr = *self.v.get_unchecked(last);
        // permute concatenates the second half of the last vector and a vector of zeros
        *self.v.get_unchecked_mut(last) = _mm256_alignr_epi8(_mm256_permute2x128_si256(curr, curr, 0b10000001i32), curr, 4i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn shift_right_1_mut(&mut self) {
        for i in (1..self.v.len()).rev() {
            let curr = *self.v.get_unchecked(i);
            // permute concatenates the second half of the previous vector and the first half of the current vector
            *self.v.get_unchecked_mut(i) = _mm256_alignr_epi8(
                curr, _mm256_permute2x128_si256(curr, *self.v.get_unchecked(i - 1), 0b00000011i32), 14i32);
        }

        // first one gets to shift in zeros
        let curr = *self.v.get_unchecked(0);
        // permute concatenates a vector of zeros and the first half of the first vector
        *self.v.get_unchecked_mut(0) = _mm256_alignr_epi8(curr, _mm256_permute2x128_si256(curr, curr, 0b00001000i32), 14i32);
    }

    operation_param2!("avx2", add, _mm256_add_epi16);
    operation_param2!("avx2", adds, _mm256_adds_epu16);
    operation_param2!("avx2", andnot, _mm256_andnot_si256);
    operation_param2!("avx2", cmpeq, _mm256_cmpeq_epi16);
    operation_param2!("avx2", min, _mm256_min_epu16);
    operation_param2!("avx2", max, _mm256_max_epu16);

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn shift_left_1(a: &Self, res: &mut Self) {
        for i in 0..(a.v.len() - 1) {
            let curr = *a.v.get_unchecked(i);
            // permute concatenates the second half of the current vector and the first half of the next vector
            *res.v.get_unchecked_mut(i) = _mm256_alignr_epi8(
                _mm256_permute2x128_si256(curr, *a.v.get_unchecked(i + 1), 0b00100001i32), curr, 2i32);
        }

        // last one gets to shift in zeros
        let last = a.v.len() - 1;
        let curr = *a.v.get_unchecked(last);
        // permute concatenates the second half of the last vector and a vector of zeros
        *res.v.get_unchecked_mut(last) = _mm256_alignr_epi8(_mm256_permute2x128_si256(curr, curr, 0b10000001i32), curr, 2i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn shift_right_1(a: &Self, res: &mut Self) {
        for i in (1..a.v.len()).rev() {
            let curr = *a.v.get_unchecked(i);
            // permute concatenates the second half of the previous vector and the first half of the current vector
            *res.v.get_unchecked_mut(i) = _mm256_alignr_epi8(
                curr, _mm256_permute2x128_si256(curr, *a.v.get_unchecked(i - 1), 0b00000011i32), 14i32);
        }

        // first one gets to shift in zeros
        let curr = *a.v.get_unchecked(0);
        // permute concatenates a vector of zeros and the first half of the first vector
        *res.v.get_unchecked_mut(0) = _mm256_alignr_epi8(curr, _mm256_permute2x128_si256(curr, curr, 0b00001000i32), 14i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn triple_argmin(sub: &Self, a_gap: &Self, b_gap: &Self, res_min: &mut Self) -> Self {
        // return the edit used in addition to doing a min operation
        let mut v = Vec::with_capacity(sub.v.len());
        let twos = _mm256_set1_epi16(2);

        for i in 0..sub.v.len() {
            let sub = *sub.v.get_unchecked(i);
            let a_gap = *a_gap.v.get_unchecked(i);
            let b_gap = *b_gap.v.get_unchecked(i);

            let res_min1 = _mm256_min_epu16(a_gap, b_gap);
            // a gap: 2 + -1 = 1, b gap: 2 + 0 = 2
            let res_arg1 = _mm256_add_epi16(twos, _mm256_cmpeq_epi16(a_gap, res_min1));

            let res_min2 = _mm256_min_epu16(sub, res_min1);
            // sub: 0
            let res_arg2 = _mm256_andnot_si256(_mm256_cmpeq_epi16(sub, res_min2), res_arg1);

            *res_min.v.get_unchecked_mut(i) = res_min2;
            v.push(res_arg2);
        }

        Self{
            v: v
        }
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn triple_min_length(sub: &Self, a_gap: &Self,
                                b_gap: &Self, sub_length: &Self, a_gap_length: &Self,
                                b_gap_length: &Self, res_min: &mut Self, res_length: &mut Self) {
        // choose the length based on which edit is chosen during the min operation
        // secondary objective of maximizing length if edit costs equal
        for i in 0..sub.v.len() {
            let sub = *sub.v.get_unchecked(i);
            let a_gap = *a_gap.v.get_unchecked(i);
            let b_gap = *b_gap.v.get_unchecked(i);
            let sub_length = *sub_length.v.get_unchecked(i);
            let a_gap_length = *a_gap_length.v.get_unchecked(i);
            let b_gap_length = *b_gap_length.v.get_unchecked(i);

            let res_min1 = _mm256_min_epu16(a_gap, b_gap);
            let a_b_gt_mask = _mm256_cmpeq_epi16(a_gap, res_min1); // a gap: -1, b gap: 0
            let mut res_length1 = _mm256_blendv_epi8(b_gap_length, a_gap_length, a_b_gt_mask); // lengths based on edits
            let a_b_eq_mask = _mm256_cmpeq_epi16(a_gap, b_gap); // equal: -1
            let a_b_max_len = _mm256_max_epu16(a_gap_length, b_gap_length);
            res_length1 = _mm256_blendv_epi8(res_length1, a_b_max_len, a_b_eq_mask); // maximize length if edits equal

            let res_min2 = _mm256_min_epu16(sub, res_min1);
            let sub_gt_mask = _mm256_cmpeq_epi16(sub, res_min2); // sub: -1, prev a or b gap: 0
            let mut res_length2 = _mm256_blendv_epi8(res_length1, sub_length, sub_gt_mask); // length based on edits
            let sub_eq_mask = _mm256_cmpeq_epi16(sub, res_min1);
            let sub_max_len = _mm256_max_epu16(sub_length, res_length1);
            res_length2 = _mm256_blendv_epi8(res_length2, sub_max_len, sub_eq_mask); // maximize length if edits equal

            *res_min.v.get_unchecked_mut(i) = res_min2;
            *res_length.v.get_unchecked_mut(i) = res_length2;
        }
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn double_min_length(new_gap: &Self, res_cont_gap: &mut Self, new_gap_length: &Self, res_cont_gap_length: &mut Self) {
        // choose the length based on which gap type is chosen during the min operation
        // secondary objective of maximizing length if edit costs equal
        for i in 0..new_gap.v.len() {
            let new_gap = *new_gap.v.get_unchecked(i);
            let cont_gap = *res_cont_gap.v.get_unchecked(i);
            let new_gap_length = *new_gap_length.v.get_unchecked(i);
            let cont_gap_length = *res_cont_gap_length.v.get_unchecked(i);

            let res_min = _mm256_min_epu16(new_gap, cont_gap);
            let new_cont_gt_mask = _mm256_cmpeq_epi16(new_gap, res_min); // new gap: -1, continue gap: 0
            let mut res_length = _mm256_blendv_epi8(cont_gap_length, new_gap_length, new_cont_gt_mask); // lengths based on edits
            let new_cont_eq_mask = _mm256_cmpeq_epi16(new_gap, cont_gap); // equal: -1
            let new_cont_max_len = _mm256_max_epu16(new_gap_length, cont_gap_length);
            res_length = _mm256_blendv_epi8(res_length, new_cont_max_len, new_cont_eq_mask); // maximize length if edits equal

            *res_cont_gap.v.get_unchecked_mut(i) = res_min;
            *res_cont_gap_length.v.get_unchecked_mut(i) = res_length;
        }
    }
}

// this implementation will probably only be used for debugging
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
impl fmt::Display for AvxNx16x16 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        #[target_feature(enable = "avx2")]
        #[inline]
        unsafe fn fmt_internal(s: &AvxNx16x16, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "[")?;

            let mut arr = [0u16; 16];
            let arr_ptr = arr.as_mut_ptr() as *mut __m256i;

            for i in 0..(s.v.len() - 1) {
                _mm256_storeu_si256(arr_ptr, *s.v.get_unchecked(i));

                for j in 0..16 {
                    write!(f, "{:>3}, ", *arr.get_unchecked(j))?;
                }
            }

            // leftover elements

            _mm256_storeu_si256(arr_ptr, *s.v.get_unchecked(s.v.len() - 1));

            let start = (s.v.len() - 1) << 4;

            for i in 0..(s.upper_bound() - start) {
                if i == s.upper_bound() - start - 1 {
                    write!(f, "{:>3}", *arr.get_unchecked(i))?;
                }else{
                    write!(f, "{:>3}, ", *arr.get_unchecked(i))?;
                }
            }

            write!(f, "]")
        }

        unsafe { fmt_internal(self, f) }
    }
}

/// N x 8 x 32 vector backed with 256-bit AVX vectors.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
pub struct AvxNx8x32 {
    v: Vec<__m256i>
}

/// Workaround for the lack of the _mm256_adds_epu32 intrinsic.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn _mm256_adds_epu32(a: __m256i, b: __m256i) -> __m256i {
    let sum = _mm256_add_epi32(a, b);
    let min = _mm256_min_epu32(a, sum);
    let eq = _mm256_cmpeq_epi32(a, min);
    // if the sum is less than a, then saturate
    // note: sum is either greater than both a and b or less than both
    _mm256_blendv_epi8(_mm256_set1_epi32(-1i32), sum, eq)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
impl Jewel for AvxNx8x32 {
    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn repeating(val: u32, len: usize) -> Self {
        let v = vec![_mm256_set1_epi32(val as i32); (len >> 3) + if (len & 7) > 0 {1} else {0}];

        Self{
            v: v
        }
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn repeating_max(len: usize) -> Self {
        let v = vec![_mm256_set1_epi32(-1i32); (len >> 3) + if (len & 7) > 0 {1} else {0}];

        Self{
            v: v
        }
    }

    #[inline]
    fn upper_bound(&self) -> usize {
        self.v.len() << 3
    }

    #[inline]
    fn static_upper_bound() -> usize {
        unimplemented!()
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn slow_loadu(&mut self, idx: usize, ptr: *const u8, len: usize, reverse: bool) {
        if len == 0 {
            return;
        }

        let mut arr = [0u32; 8];
        let arr_ptr = arr.as_mut_ptr() as *mut __m256i;
        let store_idx = if reverse {7} else {0};
        let load_idx = if reverse {0} else {7};

        for i in 0..len {
            let curr_idx = if reverse {idx - i} else {idx + i};
            let arr_idx = curr_idx & 7;

            if arr_idx == store_idx || i == 0 {
                _mm256_storeu_si256(arr_ptr, *self.v.get_unchecked(curr_idx >> 3));
            }

            *arr.get_unchecked_mut(arr_idx) = *ptr.offset(i as isize) as u32;

            if arr_idx == load_idx || i == len - 1 {
                *self.v.get_unchecked_mut(curr_idx >> 3) = _mm256_loadu_si256(arr_ptr);
            }
        }
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn slow_extract(&self, i: usize) -> u32 {
        let idx = i >> 3;
        let j = i & 7;
        let mut arr = [0u32; 8];
        _mm256_storeu_si256(arr.as_mut_ptr() as *mut __m256i, *self.v.get_unchecked(idx));
        *arr.get_unchecked(j)
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn slow_insert(&mut self, i: usize, val: u32) {
        let idx = i >> 3;
        let j = i & 7;
        let mut arr = [0u32; 8];
        let arr_ptr = arr.as_mut_ptr() as *mut __m256i;
        _mm256_storeu_si256(arr_ptr, *self.v.get_unchecked(idx));
        *arr.get_unchecked_mut(j) = val;
        *self.v.get_unchecked_mut(idx) = _mm256_loadu_si256(arr_ptr);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn insert_last_0(&mut self, val: u32) {
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm256_insert_epi32(*self.v.get_unchecked(last), val as i32, 7i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn insert_last_1(&mut self, val: u32) {
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm256_insert_epi32(*self.v.get_unchecked(last), val as i32, 6i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn insert_last_2(&mut self, val: u32) {
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm256_insert_epi32(*self.v.get_unchecked(last), val as i32, 5i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn insert_last_max(&mut self) {
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm256_insert_epi32(*self.v.get_unchecked(last), -1i32, 7i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn insert_first(&mut self, val: u32) {
        *self.v.get_unchecked_mut(0) = _mm256_insert_epi32(*self.v.get_unchecked(0), val as i32, 0i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn insert_first_max(&mut self) {
        *self.v.get_unchecked_mut(0) = _mm256_insert_epi32(*self.v.get_unchecked(0), -1i32, 0i32);
    }

    operation_mut_param2!("avx2", add_mut, _mm256_add_epi32);
    operation_mut_param2!("avx2", adds_mut, _mm256_adds_epu32);
    operation_mut_param2!("avx2", and_mut, _mm256_and_si256);
    operation_mut_param2!("avx2", andnot_mut, _mm256_andnot_si256);
    operation_mut_param2!("avx2", cmpeq_mut, _mm256_cmpeq_epi32);
    operation_mut_param2!("avx2", min_mut, _mm256_min_epu32);
    operation_mut_param2!("avx2", max_mut, _mm256_max_epu32);

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn blendv_mut(&mut self, b: &Self, mask: &Self) {
        for i in 0..self.v.len() {
            *self.v.get_unchecked_mut(i) = _mm256_blendv_epi8(*self.v.get_unchecked(i), *b.v.get_unchecked(i), *mask.v.get_unchecked(i));
        }
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn shift_left_1_mut(&mut self) {
        for i in 0..(self.v.len() - 1) {
            let curr = *self.v.get_unchecked(i);
            // permute concatenates the second half of the current vector and the first half of the next vector
            *self.v.get_unchecked_mut(i) = _mm256_alignr_epi8(
                _mm256_permute2x128_si256(curr, *self.v.get_unchecked(i + 1), 0b00100001i32), curr, 4i32);
        }

        // last one gets to shift in zeros
        let last = self.v.len() - 1;
        let curr = *self.v.get_unchecked(last);
        // permute concatenates the second half of the last vector and a vector of zeros
        *self.v.get_unchecked_mut(last) = _mm256_alignr_epi8(_mm256_permute2x128_si256(curr, curr, 0b10000001i32), curr, 4i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn shift_left_2_mut(&mut self) {
        for i in 0..(self.v.len() - 1) {
            let curr = *self.v.get_unchecked(i);
            // permute concatenates the second half of the current vector and the first half of the next vector
            *self.v.get_unchecked_mut(i) = _mm256_alignr_epi8(
                _mm256_permute2x128_si256(curr, *self.v.get_unchecked(i + 1), 0b00100001i32), curr, 8i32);
        }

        // last one gets to shift in zeros
        let last = self.v.len() - 1;
        let curr = *self.v.get_unchecked(last);
        // permute concatenates the second half of the last vector and a vector of zeros
        *self.v.get_unchecked_mut(last) = _mm256_alignr_epi8(_mm256_permute2x128_si256(curr, curr, 0b10000001i32), curr, 8i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn shift_right_1_mut(&mut self) {
        for i in (1..self.v.len()).rev() {
            let curr = *self.v.get_unchecked(i);
            // permute concatenates the second half of the previous vector and the first half of the current vector
            *self.v.get_unchecked_mut(i) = _mm256_alignr_epi8(
                curr, _mm256_permute2x128_si256(curr, *self.v.get_unchecked(i - 1), 0b00000011i32), 12i32);
        }

        // first one gets to shift in zeros
        let curr = *self.v.get_unchecked(0);
        // permute concatenates a vector of zeros and the first half of the first vector
        *self.v.get_unchecked_mut(0) = _mm256_alignr_epi8(curr, _mm256_permute2x128_si256(curr, curr, 0b00001000i32), 12i32);
    }

    operation_param2!("avx2", add, _mm256_add_epi32);
    operation_param2!("avx2", adds, _mm256_adds_epu32);
    operation_param2!("avx2", andnot, _mm256_andnot_si256);
    operation_param2!("avx2", cmpeq, _mm256_cmpeq_epi32);
    operation_param2!("avx2", min, _mm256_min_epu32);
    operation_param2!("avx2", max, _mm256_max_epu32);

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn shift_left_1(a: &Self, res: &mut Self) {
        for i in 0..(a.v.len() - 1) {
            let curr = *a.v.get_unchecked(i);
            // permute concatenates the second half of the current vector and the first half of the next vector
            *res.v.get_unchecked_mut(i) = _mm256_alignr_epi8(
                _mm256_permute2x128_si256(curr, *a.v.get_unchecked(i + 1), 0b00100001i32), curr, 4i32);
        }

        // last one gets to shift in zeros
        let last = a.v.len() - 1;
        let curr = *a.v.get_unchecked(last);
        // permute concatenates the second half of the last vector and a vector of zeros
        *res.v.get_unchecked_mut(last) = _mm256_alignr_epi8(_mm256_permute2x128_si256(curr, curr, 0b10000001i32), curr, 4i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn shift_right_1(a: &Self, res: &mut Self) {
        for i in (1..a.v.len()).rev() {
            let curr = *a.v.get_unchecked(i);
            // permute concatenates the second half of the previous vector and the first half of the current vector
            *res.v.get_unchecked_mut(i) = _mm256_alignr_epi8(
                curr, _mm256_permute2x128_si256(curr, *a.v.get_unchecked(i - 1), 0b00000011i32), 12i32);
        }

        // first one gets to shift in zeros
        let curr = *a.v.get_unchecked(0);
        // permute concatenates a vector of zeros and the first half of the first vector
        *res.v.get_unchecked_mut(0) = _mm256_alignr_epi8(curr, _mm256_permute2x128_si256(curr, curr, 0b00001000i32), 12i32);
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn triple_argmin(sub: &Self, a_gap: &Self, b_gap: &Self, res_min: &mut Self) -> Self {
        // return the edit used in addition to doing a min operation
        let mut v = Vec::with_capacity(sub.v.len());
        let twos = _mm256_set1_epi32(2);

        for i in 0..sub.v.len() {
            let sub = *sub.v.get_unchecked(i);
            let a_gap = *a_gap.v.get_unchecked(i);
            let b_gap = *b_gap.v.get_unchecked(i);

            let res_min1 = _mm256_min_epu32(a_gap, b_gap);
            // a gap: 2 + -1 = 1, b gap: 2 + 0 = 2
            let res_arg1 = _mm256_add_epi32(twos, _mm256_cmpeq_epi32(a_gap, res_min1));

            let res_min2 = _mm256_min_epu32(sub, res_min1);
            // sub: 0
            let res_arg2 = _mm256_andnot_si256(_mm256_cmpeq_epi32(sub, res_min2), res_arg1);

            *res_min.v.get_unchecked_mut(i) = res_min2;
            v.push(res_arg2);
        }

        Self{
            v: v
        }
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn triple_min_length(sub: &Self, a_gap: &Self,
                                b_gap: &Self, sub_length: &Self, a_gap_length: &Self,
                                b_gap_length: &Self, res_min: &mut Self, res_length: &mut Self) {
        // choose the length based on which edit is chosen during the min operation
        // secondary objective of maximizing length if edit costs equal
        for i in 0..sub.v.len() {
            let sub = *sub.v.get_unchecked(i);
            let a_gap = *a_gap.v.get_unchecked(i);
            let b_gap = *b_gap.v.get_unchecked(i);
            let sub_length = *sub_length.v.get_unchecked(i);
            let a_gap_length = *a_gap_length.v.get_unchecked(i);
            let b_gap_length = *b_gap_length.v.get_unchecked(i);

            let res_min1 = _mm256_min_epu32(a_gap, b_gap);
            let a_b_gt_mask = _mm256_cmpeq_epi32(a_gap, res_min1); // a gap: -1, b gap: 0
            let mut res_length1 = _mm256_blendv_epi8(b_gap_length, a_gap_length, a_b_gt_mask); // lengths based on edits
            let a_b_eq_mask = _mm256_cmpeq_epi32(a_gap, b_gap); // equal: -1
            let a_b_max_len = _mm256_max_epu32(a_gap_length, b_gap_length);
            res_length1 = _mm256_blendv_epi8(res_length1, a_b_max_len, a_b_eq_mask); // maximize length if edits equal

            let res_min2 = _mm256_min_epu32(sub, res_min1);
            let sub_gt_mask = _mm256_cmpeq_epi32(sub, res_min2); // sub: -1, prev a or b gap: 0
            let mut res_length2 = _mm256_blendv_epi8(res_length1, sub_length, sub_gt_mask); // length based on edits
            let sub_eq_mask = _mm256_cmpeq_epi32(sub, res_min1);
            let sub_max_len = _mm256_max_epu32(sub_length, res_length1);
            res_length2 = _mm256_blendv_epi8(res_length2, sub_max_len, sub_eq_mask); // maximize length if edits equal

            *res_min.v.get_unchecked_mut(i) = res_min2;
            *res_length.v.get_unchecked_mut(i) = res_length2;
        }
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn double_min_length(new_gap: &Self, res_cont_gap: &mut Self, new_gap_length: &Self, res_cont_gap_length: &mut Self) {
        // choose the length based on which gap type is chosen during the min operation
        // secondary objective of maximizing length if edit costs equal
        for i in 0..new_gap.v.len() {
            let new_gap = *new_gap.v.get_unchecked(i);
            let cont_gap = *res_cont_gap.v.get_unchecked(i);
            let new_gap_length = *new_gap_length.v.get_unchecked(i);
            let cont_gap_length = *res_cont_gap_length.v.get_unchecked(i);

            let res_min = _mm256_min_epu32(new_gap, cont_gap);
            let new_cont_gt_mask = _mm256_cmpeq_epi32(new_gap, res_min); // new gap: -1, continue gap: 0
            let mut res_length = _mm256_blendv_epi8(cont_gap_length, new_gap_length, new_cont_gt_mask); // lengths based on edits
            let new_cont_eq_mask = _mm256_cmpeq_epi32(new_gap, cont_gap); // equal: -1
            let new_cont_max_len = _mm256_max_epu32(new_gap_length, cont_gap_length);
            res_length = _mm256_blendv_epi8(res_length, new_cont_max_len, new_cont_eq_mask); // maximize length if edits equal

            *res_cont_gap.v.get_unchecked_mut(i) = res_min;
            *res_cont_gap_length.v.get_unchecked_mut(i) = res_length;
        }
    }
}

// this implementation will probably only be used for debugging
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
impl fmt::Display for AvxNx8x32 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        #[target_feature(enable = "avx2")]
        #[inline]
        unsafe fn fmt_internal(s: &AvxNx8x32, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "[")?;

            let mut arr = [0u32; 8];
            let arr_ptr = arr.as_mut_ptr() as *mut __m256i;

            for i in 0..(s.v.len() - 1) {
                _mm256_storeu_si256(arr_ptr, *s.v.get_unchecked(i));

                for j in 0..8 {
                    write!(f, "{:>3}, ", *arr.get_unchecked(j))?;
                }
            }

            // leftover elements

            _mm256_storeu_si256(arr_ptr, *s.v.get_unchecked(s.v.len() - 1));

            let start = (s.v.len() - 1) << 3;

            for i in 0..(s.upper_bound() - start) {
                if i == s.upper_bound() - start - 1 {
                    write!(f, "{:>3}", *arr.get_unchecked(i))?;
                }else{
                    write!(f, "{:>3}, ", *arr.get_unchecked(i))?;
                }
            }

            return write!(f, "]");
        }

        unsafe { fmt_internal(self, f) }
    }
}

/// N x 16 x 8 vector backed with 128-bit SSE vectors.
macro_rules! create_sse_nx16x8 {
    ($name:ident, $num:literal) => {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        pub struct $name {
            v: [__m128i; $num]
        }

        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        impl Jewel for $name {
            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn repeating(val: u32, _len: usize) -> Self {
                let v = [_mm_set1_epi8(val as i8); $num];

                Self{
                    v: v
                }
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn repeating_max(_len: usize) -> Self {
                let v = [_mm_set1_epi8(-1i8); $num];

                Self{
                    v: v
                }
            }

            #[inline]
            fn upper_bound(&self) -> usize {
                self.v.len() << 4
            }

            #[inline]
            fn static_upper_bound() -> usize {
                $num << 4
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn slow_loadu(&mut self, idx: usize, ptr: *const u8, len: usize, reverse: bool) {
                if len == 0 {
                    return;
                }

                let mut arr = [0u8; 16];
                let arr_ptr = arr.as_mut_ptr() as *mut __m128i;
                let store_idx = if reverse {15} else {0};
                let load_idx = if reverse {0} else {15};

                for i in 0..len {
                    let curr_idx = if reverse {idx - i} else {idx + i};
                    let arr_idx = curr_idx & 15;

                    if arr_idx == store_idx || i == 0 {
                        _mm_storeu_si128(arr_ptr, *self.v.get_unchecked(curr_idx >> 4));
                    }

                    *arr.get_unchecked_mut(arr_idx) = *ptr.offset(i as isize);

                    if arr_idx == load_idx || i == len - 1 {
                        *self.v.get_unchecked_mut(curr_idx >> 4) = _mm_loadu_si128(arr_ptr);
                    }
                }
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn slow_extract(&self, i: usize) -> u32 {
                let idx = i >> 4;
                let j = i & 15;
                let mut arr = [0u8; 16];
                _mm_storeu_si128(arr.as_mut_ptr() as *mut __m128i, *self.v.get_unchecked(idx));
                *arr.get_unchecked(j) as u32
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn slow_insert(&mut self, i: usize, val: u32) {
                let idx = i >> 4;
                let j = i & 15;
                let mut arr = [0u8; 16];
                let arr_ptr = arr.as_mut_ptr() as *mut __m128i;
                _mm_storeu_si128(arr_ptr, *self.v.get_unchecked(idx));
                *arr.get_unchecked_mut(j) = val as u8;
                *self.v.get_unchecked_mut(idx) = _mm_loadu_si128(arr_ptr);
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn insert_last_0(&mut self, val: u32) {
                let last = self.v.len() - 1;
                *self.v.get_unchecked_mut(last) = _mm_insert_epi8(*self.v.get_unchecked(last), val as i32, 15i32);
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn insert_last_1(&mut self, val: u32) {
                let last = self.v.len() - 1;
                *self.v.get_unchecked_mut(last) = _mm_insert_epi8(*self.v.get_unchecked(last), val as i32, 14i32);
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn insert_last_2(&mut self, val: u32) {
                let last = self.v.len() - 1;
                *self.v.get_unchecked_mut(last) = _mm_insert_epi8(*self.v.get_unchecked(last), val as i32, 13i32);
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn insert_last_max(&mut self) {
                let last = self.v.len() - 1;
                *self.v.get_unchecked_mut(last) = _mm_insert_epi8(*self.v.get_unchecked(last), u8::MAX as i32, 15i32);
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn insert_first(&mut self, val: u32) {
                *self.v.get_unchecked_mut(0) = _mm_insert_epi8(*self.v.get_unchecked(0), val as i32, 0i32);
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn insert_first_max(&mut self) {
                *self.v.get_unchecked_mut(0) = _mm_insert_epi8(*self.v.get_unchecked(0), u8::MAX as i32, 0i32);
            }

            operation_mut_param2!("sse4.1", add_mut, _mm_add_epi8);
            operation_mut_param2!("sse4.1", adds_mut, _mm_adds_epu8);
            operation_mut_param2!("sse4.1", and_mut, _mm_and_si128);
            operation_mut_param2!("sse4.1", andnot_mut, _mm_andnot_si128);
            operation_mut_param2!("sse4.1", cmpeq_mut, _mm_cmpeq_epi8);
            operation_mut_param2!("sse4.1", min_mut, _mm_min_epu8);
            operation_mut_param2!("sse4.1", max_mut, _mm_max_epu8);

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn blendv_mut(&mut self, b: &Self, mask: &Self) {
                for i in 0..self.v.len() {
                    *self.v.get_unchecked_mut(i) = _mm_blendv_epi8(*self.v.get_unchecked(i), *b.v.get_unchecked(i), *mask.v.get_unchecked(i));
                }
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn shift_left_1_mut(&mut self) {
                for i in 0..(self.v.len() - 1) {
                    *self.v.get_unchecked_mut(i) = _mm_alignr_epi8(*self.v.get_unchecked(i + 1), *self.v.get_unchecked(i), 1i32);
                }

                // last one gets to shift in zeros
                let last = self.v.len() - 1;
                *self.v.get_unchecked_mut(last) = _mm_srli_si128(*self.v.get_unchecked(last), 1i32);
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn shift_left_2_mut(&mut self) {
                for i in 0..(self.v.len() - 1) {
                    *self.v.get_unchecked_mut(i) = _mm_alignr_epi8(*self.v.get_unchecked(i + 1), *self.v.get_unchecked(i), 2i32);
                }

                // last one gets to shift in zeros
                let last = self.v.len() - 1;
                *self.v.get_unchecked_mut(last) = _mm_srli_si128(*self.v.get_unchecked(last), 2i32);
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn shift_right_1_mut(&mut self) {
                for i in (1..self.v.len()).rev() {
                    *self.v.get_unchecked_mut(i) = _mm_alignr_epi8(*self.v.get_unchecked(i), *self.v.get_unchecked(i - 1), 15i32);
                }

                // first one gets to shift in zeros
                *self.v.get_unchecked_mut(0) = _mm_slli_si128(*self.v.get_unchecked(0), 1i32);
            }

            operation_param2!("sse4.1", add, _mm_add_epi8);
            operation_param2!("sse4.1", adds, _mm_adds_epu8);
            operation_param2!("sse4.1", andnot, _mm_andnot_si128);
            operation_param2!("sse4.1", cmpeq, _mm_cmpeq_epi8);
            operation_param2!("sse4.1", min, _mm_min_epu8);
            operation_param2!("sse4.1", max, _mm_max_epu8);

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn shift_left_1(a: &Self, res: &mut Self) {
                for i in 0..(a.v.len() - 1) {
                    *res.v.get_unchecked_mut(i) = _mm_alignr_epi8(*a.v.get_unchecked(i + 1), *a.v.get_unchecked(i), 1i32);
                }

                // last one gets to shift in zeros
                let last = a.v.len() - 1;
                *res.v.get_unchecked_mut(last) = _mm_srli_si128(*a.v.get_unchecked(last), 1i32);
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn shift_right_1(a: &Self, res: &mut Self) {
                for i in (1..a.v.len()).rev() {
                    *res.v.get_unchecked_mut(i) = _mm_alignr_epi8(*a.v.get_unchecked(i), *a.v.get_unchecked(i - 1), 15i32);
                }

                // first one gets to shift in zeros
                *res.v.get_unchecked_mut(0) = _mm_slli_si128(*a.v.get_unchecked(0), 1i32);
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn triple_argmin(sub: &Self, a_gap: &Self, b_gap: &Self, res_min: &mut Self) -> Self {
                // return the edit used in addition to doing a min operation
                let mut v = [_mm_undefined_si128(); $num];
                let twos = _mm_set1_epi8(2);

                for i in 0..sub.v.len() {
                    let sub = *sub.v.get_unchecked(i);
                    let a_gap = *a_gap.v.get_unchecked(i);
                    let b_gap = *b_gap.v.get_unchecked(i);

                    let res_min1 = _mm_min_epu8(a_gap, b_gap);
                    // a gap: 2 + -1 = 1, b gap: 2 + 0 = 2
                    let res_arg1 = _mm_add_epi8(twos, _mm_cmpeq_epi8(a_gap, res_min1));

                    let res_min2 = _mm_min_epu8(sub, res_min1);
                    // sub: 0
                    let res_arg2 = _mm_andnot_si128(_mm_cmpeq_epi8(sub, res_min2), res_arg1);

                    *res_min.v.get_unchecked_mut(i) = res_min2;
                    *v.get_unchecked_mut(i) = res_arg2;
                }

                Self{
                    v: v
                }
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn triple_min_length(sub: &Self, a_gap: &Self,
                                        b_gap: &Self, sub_length: &Self, a_gap_length: &Self,
                                        b_gap_length: &Self, res_min: &mut Self, res_length: &mut Self) {
                // choose the length based on which edit is chosen during the min operation
                // secondary objective of maximizing length if edit costs equal
                for i in 0..sub.v.len() {
                    let sub = *sub.v.get_unchecked(i);
                    let a_gap = *a_gap.v.get_unchecked(i);
                    let b_gap = *b_gap.v.get_unchecked(i);
                    let sub_length = *sub_length.v.get_unchecked(i);
                    let a_gap_length = *a_gap_length.v.get_unchecked(i);
                    let b_gap_length = *b_gap_length.v.get_unchecked(i);

                    let res_min1 = _mm_min_epu8(a_gap, b_gap);
                    let a_b_gt_mask = _mm_cmpeq_epi8(a_gap, res_min1); // a gap: -1, b gap: 0
                    let mut res_length1 = _mm_blendv_epi8(b_gap_length, a_gap_length, a_b_gt_mask); // lengths based on edits
                    let a_b_eq_mask = _mm_cmpeq_epi8(a_gap, b_gap); // equal: -1
                    let a_b_max_len = _mm_max_epu8(a_gap_length, b_gap_length);
                    res_length1 = _mm_blendv_epi8(res_length1, a_b_max_len, a_b_eq_mask); // maximize length if edits equal

                    let res_min2 = _mm_min_epu8(sub, res_min1);
                    let sub_gt_mask = _mm_cmpeq_epi8(sub, res_min2); // sub: -1, prev a or b gap: 0
                    let mut res_length2 = _mm_blendv_epi8(res_length1, sub_length, sub_gt_mask); // length based on edits
                    let sub_eq_mask = _mm_cmpeq_epi8(sub, res_min1);
                    let sub_max_len = _mm_max_epu8(sub_length, res_length1);
                    res_length2 = _mm_blendv_epi8(res_length2, sub_max_len, sub_eq_mask); // maximize length if edits equal

                    *res_min.v.get_unchecked_mut(i) = res_min2;
                    *res_length.v.get_unchecked_mut(i) = res_length2;
                }
            }

            #[target_feature(enable = "sse4.1")]
            #[inline]
            unsafe fn double_min_length(new_gap: &Self, res_cont_gap: &mut Self, new_gap_length: &Self, res_cont_gap_length: &mut Self) {
                // choose the length based on which gap type is chosen during the min operation
                // secondary objective of maximizing length if edit costs equal
                for i in 0..new_gap.v.len() {
                    let new_gap = *new_gap.v.get_unchecked(i);
                    let cont_gap = *res_cont_gap.v.get_unchecked(i);
                    let new_gap_length = *new_gap_length.v.get_unchecked(i);
                    let cont_gap_length = *res_cont_gap_length.v.get_unchecked(i);

                    let res_min = _mm_min_epu8(new_gap, cont_gap);
                    let new_cont_gt_mask = _mm_cmpeq_epi8(new_gap, res_min); // new gap: -1, continue gap: 0
                    let mut res_length = _mm_blendv_epi8(cont_gap_length, new_gap_length, new_cont_gt_mask); // lengths based on edits
                    let new_cont_eq_mask = _mm_cmpeq_epi8(new_gap, cont_gap); // equal: -1
                    let new_cont_max_len = _mm_max_epu8(new_gap_length, cont_gap_length);
                    res_length = _mm_blendv_epi8(res_length, new_cont_max_len, new_cont_eq_mask); // maximize length if edits equal

                    *res_cont_gap.v.get_unchecked_mut(i) = res_min;
                    *res_cont_gap_length.v.get_unchecked_mut(i) = res_length;
                }
            }
        }

        // this implementation will probably only be used for debugging
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
                #[target_feature(enable = "sse4.1")]
                #[inline]
                unsafe fn fmt_internal(s: &$name, f: &mut fmt::Formatter) -> fmt::Result {
                    write!(f, "[")?;

                    let mut arr = [0u8; 16];
                    let arr_ptr = arr.as_mut_ptr() as *mut __m128i;

                    for i in 0..(s.v.len() - 1) {
                        _mm_storeu_si128(arr_ptr, *s.v.get_unchecked(i));

                        for j in 0..16 {
                            write!(f, "{:>3}, ", *arr.get_unchecked(j))?;
                        }
                    }

                    // leftover elements

                    _mm_storeu_si128(arr_ptr, *s.v.get_unchecked(s.v.len() - 1));

                    let start = (s.v.len() - 1) << 4;

                    for i in 0..(s.upper_bound() - start) {
                        if i == s.upper_bound() - start - 1 {
                            write!(f, "{:>3}", *arr.get_unchecked(i))?;
                        }else{
                            write!(f, "{:>3}, ", *arr.get_unchecked(i))?;
                        }
                    }

                    write!(f, "]")
                }

                unsafe { fmt_internal(self, f) }
            }
        }
    };
}

// constant array size, so the compiler should unroll the loops
create_sse_nx16x8!(Sse1x16x8, 1);
create_sse_nx16x8!(Sse2x16x8, 2);
create_sse_nx16x8!(Sse4x16x8, 4);
create_sse_nx16x8!(Sse8x16x8, 8);
create_sse_nx16x8!(Sse16x16x8, 16);

/// N x 8 x 16 vector backed with 128-bit SSE vectors.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
pub struct SseNx8x16 {
    v: Vec<__m128i>
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
impl Jewel for SseNx8x16 {
    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn repeating(val: u32, len: usize) -> Self {
        let v = vec![_mm_set1_epi16(val as i16); (len >> 3) + if (len & 7) > 0 {1} else {0}];

        Self{
            v: v
        }
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn repeating_max(len: usize) -> Self {
        let v = vec![_mm_set1_epi16(-1i16); (len >> 3) + if (len & 7) > 0 {1} else {0}];

        Self{
            v: v
        }
    }

    #[inline]
    fn upper_bound(&self) -> usize {
        self.v.len() << 3
    }

    #[inline]
    fn static_upper_bound() -> usize {
        unimplemented!()
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn slow_loadu(&mut self, idx: usize, ptr: *const u8, len: usize, reverse: bool) {
        if len == 0 {
            return;
        }

        let mut arr = [0u16; 8];
        let arr_ptr = arr.as_mut_ptr() as *mut __m128i;
        let store_idx = if reverse {7} else {0};
        let load_idx = if reverse {0} else {7};

        for i in 0..len {
            let curr_idx = if reverse {idx - i} else {idx + i};
            let arr_idx = curr_idx & 7;

            if arr_idx == store_idx || i == 0 {
                _mm_storeu_si128(arr_ptr, *self.v.get_unchecked(curr_idx >> 3));
            }

            *arr.get_unchecked_mut(arr_idx) = *ptr.offset(i as isize) as u16;

            if arr_idx == load_idx || i == len - 1 {
                *self.v.get_unchecked_mut(curr_idx >> 3) = _mm_loadu_si128(arr_ptr);
            }
        }
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn slow_extract(&self, i: usize) -> u32 {
        let idx = i >> 3;
        let j = i & 7;
        let mut arr = [0u16; 8];
        _mm_storeu_si128(arr.as_mut_ptr() as *mut __m128i, *self.v.get_unchecked(idx));
        *arr.get_unchecked(j) as u32
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn slow_insert(&mut self, i: usize, val: u32) {
        let idx = i >> 3;
        let j = i & 7;
        let mut arr = [0u16; 8];
        let arr_ptr = arr.as_mut_ptr() as *mut __m128i;
        _mm_storeu_si128(arr_ptr, *self.v.get_unchecked(idx));
        *arr.get_unchecked_mut(j) = val as u16;
        *self.v.get_unchecked_mut(idx) = _mm_loadu_si128(arr_ptr);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn insert_last_0(&mut self, val: u32) {
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm_insert_epi16(*self.v.get_unchecked(last), val as i32, 7i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn insert_last_1(&mut self, val: u32) {
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm_insert_epi16(*self.v.get_unchecked(last), val as i32, 6i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn insert_last_2(&mut self, val: u32) {
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm_insert_epi16(*self.v.get_unchecked(last), val as i32, 5i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn insert_last_max(&mut self) {
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm_insert_epi16(*self.v.get_unchecked(last), u16::MAX as i32, 7i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn insert_first(&mut self, val: u32) {
        *self.v.get_unchecked_mut(0) = _mm_insert_epi16(*self.v.get_unchecked(0), val as i32, 0i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn insert_first_max(&mut self) {
        *self.v.get_unchecked_mut(0) = _mm_insert_epi16(*self.v.get_unchecked(0), u16::MAX as i32, 0i32);
    }

    operation_mut_param2!("sse4.1", add_mut, _mm_add_epi16);
    operation_mut_param2!("sse4.1", adds_mut, _mm_adds_epu16);
    operation_mut_param2!("sse4.1", and_mut, _mm_and_si128);
    operation_mut_param2!("sse4.1", andnot_mut, _mm_andnot_si128);
    operation_mut_param2!("sse4.1", cmpeq_mut, _mm_cmpeq_epi16);
    operation_mut_param2!("sse4.1", min_mut, _mm_min_epu16);
    operation_mut_param2!("sse4.1", max_mut, _mm_max_epu16);

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn blendv_mut(&mut self, b: &Self, mask: &Self) {
        for i in 0..self.v.len() {
            *self.v.get_unchecked_mut(i) = _mm_blendv_epi8(*self.v.get_unchecked(i), *b.v.get_unchecked(i), *mask.v.get_unchecked(i));
        }
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn shift_left_1_mut(&mut self) {
        for i in 0..(self.v.len() - 1) {
            *self.v.get_unchecked_mut(i) = _mm_alignr_epi8(*self.v.get_unchecked(i + 1), *self.v.get_unchecked(i), 2i32);
        }

        // last one gets to shift in zeros
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm_srli_si128(*self.v.get_unchecked(last), 2i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn shift_left_2_mut(&mut self) {
        for i in 0..(self.v.len() - 1) {
            *self.v.get_unchecked_mut(i) = _mm_alignr_epi8(*self.v.get_unchecked(i + 1), *self.v.get_unchecked(i), 4i32);
        }

        // last one gets to shift in zeros
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm_srli_si128(*self.v.get_unchecked(last), 4i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn shift_right_1_mut(&mut self) {
        for i in (1..self.v.len()).rev() {
            *self.v.get_unchecked_mut(i) = _mm_alignr_epi8(*self.v.get_unchecked(i), *self.v.get_unchecked(i - 1), 14i32);
        }

        // first one gets to shift in zeros
        *self.v.get_unchecked_mut(0) = _mm_slli_si128(*self.v.get_unchecked(0), 2i32);
    }

    operation_param2!("sse4.1", add, _mm_add_epi16);
    operation_param2!("sse4.1", adds, _mm_adds_epu16);
    operation_param2!("sse4.1", andnot, _mm_andnot_si128);
    operation_param2!("sse4.1", cmpeq, _mm_cmpeq_epi16);
    operation_param2!("sse4.1", min, _mm_min_epu16);
    operation_param2!("sse4.1", max, _mm_max_epu16);

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn shift_left_1(a: &Self, res: &mut Self) {
        for i in 0..(a.v.len() - 1) {
            *res.v.get_unchecked_mut(i) = _mm_alignr_epi8(*a.v.get_unchecked(i + 1), *a.v.get_unchecked(i), 2i32);
        }

        // last one gets to shift in zeros
        let last = a.v.len() - 1;
        *res.v.get_unchecked_mut(last) = _mm_srli_si128(*a.v.get_unchecked(last), 2i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn shift_right_1(a: &Self, res: &mut Self) {
        for i in (1..a.v.len()).rev() {
            *res.v.get_unchecked_mut(i) = _mm_alignr_epi8(*a.v.get_unchecked(i), *a.v.get_unchecked(i - 1), 14i32);
        }

        // first one gets to shift in zeros
        *res.v.get_unchecked_mut(0) = _mm_slli_si128(*a.v.get_unchecked(0), 2i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn triple_argmin(sub: &Self, a_gap: &Self, b_gap: &Self, res_min: &mut Self) -> Self {
        // return the edit used in addition to doing a min operation
        let mut v = Vec::with_capacity(sub.v.len());
        let twos = _mm_set1_epi16(2);

        for i in 0..sub.v.len() {
            let sub = *sub.v.get_unchecked(i);
            let a_gap = *a_gap.v.get_unchecked(i);
            let b_gap = *b_gap.v.get_unchecked(i);

            let res_min1 = _mm_min_epu16(a_gap, b_gap);
            // a gap: 2 + -1 = 1, b gap: 2 + 0 = 2
            let res_arg1 = _mm_add_epi16(twos, _mm_cmpeq_epi16(a_gap, res_min1));

            let res_min2 = _mm_min_epu16(sub, res_min1);
            // sub: 0
            let res_arg2 = _mm_andnot_si128(_mm_cmpeq_epi16(sub, res_min2), res_arg1);

            *res_min.v.get_unchecked_mut(i) = res_min2;
            v.push(res_arg2);
        }

        Self{
            v: v
        }
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn triple_min_length(sub: &Self, a_gap: &Self,
                                b_gap: &Self, sub_length: &Self, a_gap_length: &Self,
                                b_gap_length: &Self, res_min: &mut Self, res_length: &mut Self) {
        // choose the length based on which edit is chosen during the min operation
        // secondary objective of maximizing length if edit costs equal
        for i in 0..sub.v.len() {
            let sub = *sub.v.get_unchecked(i);
            let a_gap = *a_gap.v.get_unchecked(i);
            let b_gap = *b_gap.v.get_unchecked(i);
            let sub_length = *sub_length.v.get_unchecked(i);
            let a_gap_length = *a_gap_length.v.get_unchecked(i);
            let b_gap_length = *b_gap_length.v.get_unchecked(i);

            let res_min1 = _mm_min_epu16(a_gap, b_gap);
            let a_b_gt_mask = _mm_cmpeq_epi16(a_gap, res_min1); // a gap: -1, b gap: 0
            let mut res_length1 = _mm_blendv_epi8(b_gap_length, a_gap_length, a_b_gt_mask); // lengths based on edits
            let a_b_eq_mask = _mm_cmpeq_epi16(a_gap, b_gap); // equal: -1
            let a_b_max_len = _mm_max_epu16(a_gap_length, b_gap_length);
            res_length1 = _mm_blendv_epi8(res_length1, a_b_max_len, a_b_eq_mask); // maximize length if edits equal

            let res_min2 = _mm_min_epu16(sub, res_min1);
            let sub_gt_mask = _mm_cmpeq_epi16(sub, res_min2); // sub: -1, prev a or b gap: 0
            let mut res_length2 = _mm_blendv_epi8(res_length1, sub_length, sub_gt_mask); // length based on edits
            let sub_eq_mask = _mm_cmpeq_epi16(sub, res_min1);
            let sub_max_len = _mm_max_epu16(sub_length, res_length1);
            res_length2 = _mm_blendv_epi8(res_length2, sub_max_len, sub_eq_mask); // maximize length if edits equal

            *res_min.v.get_unchecked_mut(i) = res_min2;
            *res_length.v.get_unchecked_mut(i) = res_length2;
        }
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn double_min_length(new_gap: &Self, res_cont_gap: &mut Self, new_gap_length: &Self, res_cont_gap_length: &mut Self) {
        // choose the length based on which gap type is chosen during the min operation
        // secondary objective of maximizing length if edit costs equal
        for i in 0..new_gap.v.len() {
            let new_gap = *new_gap.v.get_unchecked(i);
            let cont_gap = *res_cont_gap.v.get_unchecked(i);
            let new_gap_length = *new_gap_length.v.get_unchecked(i);
            let cont_gap_length = *res_cont_gap_length.v.get_unchecked(i);

            let res_min = _mm_min_epu16(new_gap, cont_gap);
            let new_cont_gt_mask = _mm_cmpeq_epi16(new_gap, res_min); // new gap: -1, continue gap: 0
            let mut res_length = _mm_blendv_epi8(cont_gap_length, new_gap_length, new_cont_gt_mask); // lengths based on edits
            let new_cont_eq_mask = _mm_cmpeq_epi16(new_gap, cont_gap); // equal: -1
            let new_cont_max_len = _mm_max_epu16(new_gap_length, cont_gap_length);
            res_length = _mm_blendv_epi8(res_length, new_cont_max_len, new_cont_eq_mask); // maximize length if edits equal

            *res_cont_gap.v.get_unchecked_mut(i) = res_min;
            *res_cont_gap_length.v.get_unchecked_mut(i) = res_length;
        }
    }
}

// this implementation will probably only be used for debugging
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
impl fmt::Display for SseNx8x16 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        #[target_feature(enable = "sse4.1")]
        #[inline]
        unsafe fn fmt_internal(s: &SseNx8x16, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "[")?;

            let mut arr = [0u16; 8];
            let arr_ptr = arr.as_mut_ptr() as *mut __m128i;

            for i in 0..(s.v.len() - 1) {
                _mm_storeu_si128(arr_ptr, *s.v.get_unchecked(i));

                for j in 0..8 {
                    write!(f, "{:>3}, ", *arr.get_unchecked(j))?;
                }
            }

            // leftover elements

            _mm_storeu_si128(arr_ptr, *s.v.get_unchecked(s.v.len() - 1));

            let start = (s.v.len() - 1) << 3;

            for i in 0..(s.upper_bound() - start) {
                if i == s.upper_bound() - start - 1 {
                    write!(f, "{:>3}", *arr.get_unchecked(i))?;
                }else{
                    write!(f, "{:>3}, ", *arr.get_unchecked(i))?;
                }
            }

            write!(f, "]")
        }

        unsafe { fmt_internal(self, f) }
    }
}

/// N x 4 x 32 vector backed with 128-bit SSE vectors.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
pub struct SseNx4x32 {
    v: Vec<__m128i>
}

/// Workaround for the lack of the _mm_adds_epu32 intrinsic.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse4.1")]
#[inline]
unsafe fn _mm_adds_epu32(a: __m128i, b: __m128i) -> __m128i {
    let sum = _mm_add_epi32(a, b);
    let min = _mm_min_epu32(a, sum);
    let eq = _mm_cmpeq_epi32(a, min);
    // if the sum is less than a, then saturate
    // note: sum is either greater than both a and b or less than both
    _mm_blendv_epi8(_mm_set1_epi32(-1i32), sum, eq)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
impl Jewel for SseNx4x32 {
    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn repeating(val: u32, len: usize) -> Self {
        let v = vec![_mm_set1_epi32(val as i32); (len >> 2) + if (len & 3) > 0 {1} else {0}];

        Self{
            v: v
        }
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn repeating_max(len: usize) -> Self {
        let v = vec![_mm_set1_epi32(-1i32); (len >> 2) + if (len & 3) > 0 {1} else {0}];

        Self{
            v: v
        }
    }

    #[inline]
    fn upper_bound(&self) -> usize {
        self.v.len() << 2
    }

    #[inline]
    fn static_upper_bound() -> usize {
        unimplemented!()
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn slow_loadu(&mut self, idx: usize, ptr: *const u8, len: usize, reverse: bool) {
        if len == 0 {
            return;
        }

        let mut arr = [0u32; 4];
        let arr_ptr = arr.as_mut_ptr() as *mut __m128i;
        let store_idx = if reverse {3} else {0};
        let load_idx = if reverse {0} else {3};

        for i in 0..len {
            let curr_idx = if reverse {idx - i} else {idx + i};
            let arr_idx = curr_idx & 3;

            if arr_idx == store_idx || i == 0 {
                _mm_storeu_si128(arr_ptr, *self.v.get_unchecked(curr_idx >> 2));
            }

            *arr.get_unchecked_mut(arr_idx) = *ptr.offset(i as isize) as u32;

            if arr_idx == load_idx || i == len - 1 {
                *self.v.get_unchecked_mut(curr_idx >> 2) = _mm_loadu_si128(arr_ptr);
            }
        }
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn slow_extract(&self, i: usize) -> u32 {
        let idx = i >> 2;
        let j = i & 3;
        let mut arr = [0u32; 4];
        _mm_storeu_si128(arr.as_mut_ptr() as *mut __m128i, *self.v.get_unchecked(idx));
        *arr.get_unchecked(j)
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn slow_insert(&mut self, i: usize, val: u32) {
        let idx = i >> 2;
        let j = i & 3;
        let mut arr = [0u32; 4];
        let arr_ptr = arr.as_mut_ptr() as *mut __m128i;
        _mm_storeu_si128(arr_ptr, *self.v.get_unchecked(idx));
        *arr.get_unchecked_mut(j) = val;
        *self.v.get_unchecked_mut(idx) = _mm_loadu_si128(arr_ptr);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn insert_last_0(&mut self, val: u32) {
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm_insert_epi32(*self.v.get_unchecked(last), val as i32, 3i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn insert_last_1(&mut self, val: u32) {
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm_insert_epi32(*self.v.get_unchecked(last), val as i32, 2i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn insert_last_2(&mut self, val: u32) {
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm_insert_epi32(*self.v.get_unchecked(last), val as i32, 1i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn insert_last_max(&mut self) {
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm_insert_epi32(*self.v.get_unchecked(last), -1i32, 3i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn insert_first(&mut self, val: u32) {
        *self.v.get_unchecked_mut(0) = _mm_insert_epi32(*self.v.get_unchecked(0), val as i32, 0i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn insert_first_max(&mut self) {
        *self.v.get_unchecked_mut(0) = _mm_insert_epi32(*self.v.get_unchecked(0), -1i32, 0i32);
    }

    operation_mut_param2!("sse4.1", add_mut, _mm_add_epi32);
    operation_mut_param2!("sse4.1", adds_mut, _mm_adds_epu32);
    operation_mut_param2!("sse4.1", and_mut, _mm_and_si128);
    operation_mut_param2!("sse4.1", andnot_mut, _mm_andnot_si128);
    operation_mut_param2!("sse4.1", cmpeq_mut, _mm_cmpeq_epi32);
    operation_mut_param2!("sse4.1", min_mut, _mm_min_epu32);
    operation_mut_param2!("sse4.1", max_mut, _mm_max_epu32);

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn blendv_mut(&mut self, b: &Self, mask: &Self) {
        for i in 0..self.v.len() {
            *self.v.get_unchecked_mut(i) = _mm_blendv_epi8(*self.v.get_unchecked(i), *b.v.get_unchecked(i), *mask.v.get_unchecked(i));
        }
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn shift_left_1_mut(&mut self) {
        for i in 0..(self.v.len() - 1) {
            *self.v.get_unchecked_mut(i) = _mm_alignr_epi8(*self.v.get_unchecked(i + 1), *self.v.get_unchecked(i), 4i32);
        }

        // last one gets to shift in zeros
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm_srli_si128(*self.v.get_unchecked(last), 4i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn shift_left_2_mut(&mut self) {
        for i in 0..(self.v.len() - 1) {
            *self.v.get_unchecked_mut(i) = _mm_alignr_epi8(*self.v.get_unchecked(i + 1), *self.v.get_unchecked(i), 8i32);
        }

        // last one gets to shift in zeros
        let last = self.v.len() - 1;
        *self.v.get_unchecked_mut(last) = _mm_srli_si128(*self.v.get_unchecked(last), 8i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn shift_right_1_mut(&mut self) {
        for i in (1..self.v.len()).rev() {
            *self.v.get_unchecked_mut(i) = _mm_alignr_epi8(*self.v.get_unchecked(i), *self.v.get_unchecked(i - 1), 12i32);
        }

        // first one gets to shift in zeros
        *self.v.get_unchecked_mut(0) = _mm_slli_si128(*self.v.get_unchecked(0), 4i32);
    }

    operation_param2!("sse4.1", add, _mm_add_epi32);
    operation_param2!("sse4.1", adds, _mm_adds_epu32);
    operation_param2!("sse4.1", andnot, _mm_andnot_si128);
    operation_param2!("sse4.1", cmpeq, _mm_cmpeq_epi32);
    operation_param2!("sse4.1", min, _mm_min_epu32);
    operation_param2!("sse4.1", max, _mm_max_epu32);

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn shift_left_1(a: &Self, res: &mut Self) {
        for i in 0..(a.v.len() - 1) {
            *res.v.get_unchecked_mut(i) = _mm_alignr_epi8(*a.v.get_unchecked(i + 1), *a.v.get_unchecked(i), 4i32);
        }

        // last one gets to shift in zeros
        let last = a.v.len() - 1;
        *res.v.get_unchecked_mut(last) = _mm_srli_si128(*a.v.get_unchecked(last), 4i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn shift_right_1(a: &Self, res: &mut Self) {
        for i in (1..a.v.len()).rev() {
            *res.v.get_unchecked_mut(i) = _mm_alignr_epi8(*a.v.get_unchecked(i), *a.v.get_unchecked(i - 1), 12i32);
        }

        // first one gets to shift in zeros
        *res.v.get_unchecked_mut(0) = _mm_slli_si128(*a.v.get_unchecked(0), 4i32);
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn triple_argmin(sub: &Self, a_gap: &Self, b_gap: &Self, res_min: &mut Self) -> Self {
        // return the edit used in addition to doing a min operation
        let mut v = Vec::with_capacity(sub.v.len());
        let twos = _mm_set1_epi32(2);

        for i in 0..sub.v.len() {
            let sub = *sub.v.get_unchecked(i);
            let a_gap = *a_gap.v.get_unchecked(i);
            let b_gap = *b_gap.v.get_unchecked(i);

            let res_min1 = _mm_min_epu32(a_gap, b_gap);
            // a gap: 2 + -1 = 1, b gap: 2 + 0 = 2
            let res_arg1 = _mm_add_epi32(twos, _mm_cmpeq_epi32(a_gap, res_min1));

            let res_min2 = _mm_min_epu32(sub, res_min1);
            // sub: 0
            let res_arg2 = _mm_andnot_si128(_mm_cmpeq_epi32(sub, res_min2), res_arg1);

            *res_min.v.get_unchecked_mut(i) = res_min2;
            v.push(res_arg2);
        }

        Self{
            v: v
        }
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn triple_min_length(sub: &Self, a_gap: &Self,
                                b_gap: &Self, sub_length: &Self, a_gap_length: &Self,
                                b_gap_length: &Self, res_min: &mut Self, res_length: &mut Self) {
        // choose the length based on which edit is chosen during the min operation
        // secondary objective of maximizing length if edit costs equal
        for i in 0..sub.v.len() {
            let sub = *sub.v.get_unchecked(i);
            let a_gap = *a_gap.v.get_unchecked(i);
            let b_gap = *b_gap.v.get_unchecked(i);
            let sub_length = *sub_length.v.get_unchecked(i);
            let a_gap_length = *a_gap_length.v.get_unchecked(i);
            let b_gap_length = *b_gap_length.v.get_unchecked(i);

            let res_min1 = _mm_min_epu32(a_gap, b_gap);
            let a_b_gt_mask = _mm_cmpeq_epi32(a_gap, res_min1); // a gap: -1, b gap: 0
            let mut res_length1 = _mm_blendv_epi8(b_gap_length, a_gap_length, a_b_gt_mask); // lengths based on edits
            let a_b_eq_mask = _mm_cmpeq_epi32(a_gap, b_gap); // equal: -1
            let a_b_max_len = _mm_max_epu32(a_gap_length, b_gap_length);
            res_length1 = _mm_blendv_epi8(res_length1, a_b_max_len, a_b_eq_mask); // maximize length if edits equal

            let res_min2 = _mm_min_epu32(sub, res_min1);
            let sub_gt_mask = _mm_cmpeq_epi32(sub, res_min2); // sub: -1, prev a or b gap: 0
            let mut res_length2 = _mm_blendv_epi8(res_length1, sub_length, sub_gt_mask); // length based on edits
            let sub_eq_mask = _mm_cmpeq_epi32(sub, res_min1);
            let sub_max_len = _mm_max_epu32(sub_length, res_length1);
            res_length2 = _mm_blendv_epi8(res_length2, sub_max_len, sub_eq_mask); // maximize length if edits equal

            *res_min.v.get_unchecked_mut(i) = res_min2;
            *res_length.v.get_unchecked_mut(i) = res_length2;
        }
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn double_min_length(new_gap: &Self, res_cont_gap: &mut Self, new_gap_length: &Self, res_cont_gap_length: &mut Self) {
        // choose the length based on which gap type is chosen during the min operation
        // secondary objective of maximizing length if edit costs equal
        for i in 0..new_gap.v.len() {
            let new_gap = *new_gap.v.get_unchecked(i);
            let cont_gap = *res_cont_gap.v.get_unchecked(i);
            let new_gap_length = *new_gap_length.v.get_unchecked(i);
            let cont_gap_length = *res_cont_gap_length.v.get_unchecked(i);

            let res_min = _mm_min_epu32(new_gap, cont_gap);
            let new_cont_gt_mask = _mm_cmpeq_epi32(new_gap, res_min); // new gap: -1, continue gap: 0
            let mut res_length = _mm_blendv_epi8(cont_gap_length, new_gap_length, new_cont_gt_mask); // lengths based on edits
            let new_cont_eq_mask = _mm_cmpeq_epi32(new_gap, cont_gap); // equal: -1
            let new_cont_max_len = _mm_max_epu32(new_gap_length, cont_gap_length);
            res_length = _mm_blendv_epi8(res_length, new_cont_max_len, new_cont_eq_mask); // maximize length if edits equal

            *res_cont_gap.v.get_unchecked_mut(i) = res_min;
            *res_cont_gap_length.v.get_unchecked_mut(i) = res_length;
        }
    }
}

// this implementation will probably only be used for debugging
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
impl fmt::Display for SseNx4x32 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        #[target_feature(enable = "sse4.1")]
        #[inline]
        unsafe fn fmt_internal(s: &SseNx4x32, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "[")?;

            let mut arr = [0u32; 4];
            let arr_ptr = arr.as_mut_ptr() as *mut __m128i;

            for i in 0..(s.v.len() - 1) {
                _mm_storeu_si128(arr_ptr, *s.v.get_unchecked(i));

                for j in 0..4 {
                    write!(f, "{:>3}, ", *arr.get_unchecked(j))?;
                }
            }

            // leftover elements

            _mm_storeu_si128(arr_ptr, *s.v.get_unchecked(s.v.len() - 1));

            let start = (s.v.len() - 1) << 2;

            for i in 0..(s.upper_bound() - start) {
                if i == s.upper_bound() - start - 1 {
                    write!(f, "{:>3}", *arr.get_unchecked(i))?;
                }else{
                    write!(f, "{:>3}, ", *arr.get_unchecked(i))?;
                }
            }

            write!(f, "]")
        }

        unsafe { fmt_internal(self, f) }
    }
}

pub trait HammingJewel {
    unsafe fn loadu(ptr: *const u8, len: usize) -> Self;
    fn upper_bound(&self) -> usize;
    unsafe fn mm_count_mismatches(a_ptr: *const u8, b_ptr: *const u8, len: usize) -> u32;
    unsafe fn count_mismatches(a_ptr: *const u8, b_ptr: *const u8, len: usize) -> u32;
    unsafe fn vector_count_mismatches(a: &Self, b_ptr: *const u8, len: usize) -> u32;
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
pub struct Avx {
    v: Vec<__m256i>
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
impl HammingJewel for Avx {
    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn loadu(ptr: *const u8, len: usize) -> Self {
        let word_len = len >> 5;
        let word_rem = len & 31;
        let mut v = Vec::with_capacity(word_len + if word_rem > 0 {1} else {0});
        let avx2_ptr = ptr as *const __m256i;

        for i in 0..word_len {
            v.push(_mm256_loadu_si256(avx2_ptr.offset(i as isize)));
        }

        if word_rem > 0 {
            let mut arr = [0u8; 32];
            let end_ptr = ptr.offset((word_len << 5) as isize);

            for i in 0..word_rem {
                *arr.get_unchecked_mut(i) = *end_ptr.offset(i as isize);
            }

            v.push(_mm256_loadu_si256(arr.as_ptr() as *const __m256i));
        }

        Self{
            v: v
        }
    }

    #[inline]
    fn upper_bound(&self) -> usize {
        self.v.len() << 5
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn mm_count_mismatches(a_ptr: *const u8, b_ptr: *const u8, len: usize) -> u32 {
        let mut res = 0u32;
        let div_len = (len >> 5) as isize;
        let avx2_a_ptr = a_ptr as *const __m256i;
        let avx2_b_ptr = b_ptr as *const __m256i;

        for i in 0..div_len {
            let a = _mm256_loadu_si256(avx2_a_ptr.offset(i));
            let b = _mm256_loadu_si256(avx2_b_ptr.offset(i));
            let eq = _mm256_cmpeq_epi8(a, b);
            // basic movemask count equal bytes
            res += _mm256_movemask_epi8(eq).count_ones();
        }

        for i in (div_len << 5)..len as isize {
            res += (*a_ptr.offset(i) == *b_ptr.offset(i)) as u32;
        }

        len as u32 - res
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn count_mismatches(a_ptr: *const u8, b_ptr: *const u8, len: usize) -> u32 {
        let refresh_len = (len / (255 * 32)) as isize;
        let zeros = _mm256_setzero_si256();
        let mut sad = zeros;
        let avx2_a_ptr = a_ptr as *const __m256i;
        let avx2_b_ptr = b_ptr as *const __m256i;

        for i in 0..refresh_len {
            let mut curr = zeros;

            for j in (i * 255)..((i + 1) * 255) {
                let a = _mm256_loadu_si256(avx2_a_ptr.offset(j));
                let b = _mm256_loadu_si256(avx2_b_ptr.offset(j));
                let eq = _mm256_cmpeq_epi8(a, b);
                curr = _mm256_sub_epi8(curr, eq); // subtract -1 = add 1 when matching
                // counting matches instead of mismatches for speed
            }

            // subtract 0 and sum up 8 bytes at once horizontally into four 64 bit ints
            // accumulate those 64 bit ints
            sad = _mm256_add_epi64(sad, _mm256_sad_epu8(curr, zeros));
        }

        let word_len = (len >> 5) as isize;
        let mut curr = zeros;

        // leftover blocks of 32 bytes
        for i in (refresh_len * 255)..word_len {
            let a = _mm256_loadu_si256(avx2_a_ptr.offset(i));
            let b = _mm256_loadu_si256(avx2_b_ptr.offset(i));
            let eq = _mm256_cmpeq_epi8(a, b);
            curr = _mm256_sub_epi8(curr, eq); // subtract -1 = add 1 when matching
        }

        sad = _mm256_add_epi64(sad, _mm256_sad_epu8(curr, zeros));
        let mut sad_arr = [0u32; 8];
        _mm256_storeu_si256(sad_arr.as_mut_ptr() as *mut __m256i, sad);
        let mut res = *sad_arr.get_unchecked(0) + *sad_arr.get_unchecked(2)
            + *sad_arr.get_unchecked(4) + *sad_arr.get_unchecked(6);

        for i in (word_len << 5)..len as isize {
            res += (*a_ptr.offset(i) == *b_ptr.offset(i)) as u32;
        }

        len as u32 - res
    }

    #[target_feature(enable = "avx2")]
    #[inline]
    unsafe fn vector_count_mismatches(a: &Self, b_ptr: *const u8, len: usize) -> u32 {
        let refresh_len = (a.v.len() / 255) as isize;
        let zeros = _mm256_setzero_si256();
        let mut sad = zeros;
        let avx2_b_ptr = b_ptr as *const __m256i;

        for i in 0..refresh_len {
            let mut curr = zeros;

            for j in (i * 255)..((i + 1) * 255) {
                let a = *a.v.get_unchecked(j as usize);
                let b = _mm256_loadu_si256(avx2_b_ptr.offset(j));
                let eq = _mm256_cmpeq_epi8(a, b);
                curr = _mm256_sub_epi8(curr, eq); // subtract -1 = add 1 when matching
                // counting matches instead of mismatches for speed
            }

            // subtract 0 and sum up 8 bytes at once horizontally into four 64 bit ints
            // accumulate those 64 bit ints
            sad = _mm256_add_epi64(sad, _mm256_sad_epu8(curr, zeros));
        }

        let mut curr = zeros;

        // leftover blocks of 32 bytes
        for i in (refresh_len * 255)..a.v.len() as isize {
            let a = *a.v.get_unchecked(i as usize);
            let b = _mm256_loadu_si256(avx2_b_ptr.offset(i));
            let eq = _mm256_cmpeq_epi8(a, b);
            curr = _mm256_sub_epi8(curr, eq); // subtract -1 = add 1 when matching
        }

        sad = _mm256_add_epi64(sad, _mm256_sad_epu8(curr, zeros));
        let mut sad_arr = [0u32; 8];
        _mm256_storeu_si256(sad_arr.as_mut_ptr() as *mut __m256i, sad);
        let res = *sad_arr.get_unchecked(0) + *sad_arr.get_unchecked(2)
            + *sad_arr.get_unchecked(4) + *sad_arr.get_unchecked(6);

        len as u32 - res
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
pub struct Sse {
    v: Vec<__m128i>
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
impl HammingJewel for Sse {
    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn loadu(ptr: *const u8, len: usize) -> Self {
        let word_len = len >> 4;
        let word_rem = len & 15;
        let mut v = Vec::with_capacity(word_len + if word_rem > 0 {1} else {0});
        let sse_ptr = ptr as *const __m128i;

        for i in 0..word_len {
            v.push(_mm_loadu_si128(sse_ptr.offset(i as isize)));
        }

        if word_rem > 0 {
            let mut arr = [0u8; 16];
            let end_ptr = ptr.offset((word_len << 4) as isize);

            for i in 0..word_rem {
                *arr.get_unchecked_mut(i) = *end_ptr.offset(i as isize);
            }

            v.push(_mm_loadu_si128(arr.as_ptr() as *const __m128i));
        }

        Self{
            v: v
        }
    }

    #[inline]
    fn upper_bound(&self) -> usize {
        self.v.len() << 4
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn mm_count_mismatches(a_ptr: *const u8, b_ptr: *const u8, len: usize) -> u32 {
        let mut res = 0u32;
        let div_len = (len >> 4) as isize;
        let sse_a_ptr = a_ptr as *const __m128i;
        let sse_b_ptr = b_ptr as *const __m128i;

        for i in 0..div_len {
            let a = _mm_loadu_si128(sse_a_ptr.offset(i));
            let b = _mm_loadu_si128(sse_b_ptr.offset(i));
            let eq = _mm_cmpeq_epi8(a, b);
            // basic movemask count equal bytes
            res += _mm_movemask_epi8(eq).count_ones();
        }

        for i in (div_len << 4)..len as isize {
            res += (*a_ptr.offset(i) == *b_ptr.offset(i)) as u32;
        }

        len as u32 - res
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn count_mismatches(a_ptr: *const u8, b_ptr: *const u8, len: usize) -> u32 {
        let refresh_len = (len / (255 * 16)) as isize;
        let zeros = _mm_setzero_si128();
        let mut sad = zeros;
        let sse_a_ptr = a_ptr as *const __m128i;
        let sse_b_ptr = b_ptr as *const __m128i;

        for i in 0..refresh_len {
            let mut curr = zeros;

            for j in (i * 255)..((i + 1) * 255) {
                let a = _mm_loadu_si128(sse_a_ptr.offset(j));
                let b = _mm_loadu_si128(sse_b_ptr.offset(j));
                let eq = _mm_cmpeq_epi8(a, b);
                curr = _mm_sub_epi8(curr, eq); // subtract -1 = add 1 when matching
                // counting matches instead of mismatches for speed
            }

            // subtract 0 and sum up 8 bytes at once horizontally into two 64 bit ints
            // accumulate those 64 bit ints
            sad = _mm_add_epi64(sad, _mm_sad_epu8(curr, zeros));
        }

        let word_len = (len >> 4) as isize;
        let mut curr = zeros;

        // leftover blocks of 16 bytes
        for i in (refresh_len * 255)..word_len {
            let a = _mm_loadu_si128(sse_a_ptr.offset(i));
            let b = _mm_loadu_si128(sse_b_ptr.offset(i));
            let eq = _mm_cmpeq_epi8(a, b);
            curr = _mm_sub_epi8(curr, eq); // subtract -1 = add 1 when matching
        }

        sad = _mm_add_epi64(sad, _mm_sad_epu8(curr, zeros));
        let mut sad_arr = [0u32; 4];
        _mm_storeu_si128(sad_arr.as_mut_ptr() as *mut __m128i, sad);
        let mut res = *sad_arr.get_unchecked(0) + *sad_arr.get_unchecked(2);

        for i in (word_len << 4)..len as isize {
            res += (*a_ptr.offset(i) == *b_ptr.offset(i)) as u32;
        }

        len as u32 - res
    }

    #[target_feature(enable = "sse4.1")]
    #[inline]
    unsafe fn vector_count_mismatches(a: &Self, b_ptr: *const u8, len: usize) -> u32 {
        let refresh_len = (a.v.len() / 255) as isize;
        let zeros = _mm_setzero_si128();
        let mut sad = zeros;
        let sse_b_ptr = b_ptr as *const __m128i;

        for i in 0..refresh_len {
            let mut curr = zeros;

            for j in (i * 255)..((i + 1) * 255) {
                let a = *a.v.get_unchecked(j as usize);
                let b = _mm_loadu_si128(sse_b_ptr.offset(j));
                let eq = _mm_cmpeq_epi8(a, b);
                curr = _mm_sub_epi8(curr, eq); // subtract -1 = add 1 when matching
                // counting matches instead of mismatches for speed
            }

            // subtract 0 and sum up 8 bytes at once horizontally into two 64 bit ints
            // accumulate those 64 bit ints
            sad = _mm_add_epi64(sad, _mm_sad_epu8(curr, zeros));
        }

        let mut curr = zeros;

        // leftover blocks of 16 bytes
        for i in (refresh_len * 255)..a.v.len() as isize {
            let a = *a.v.get_unchecked(i as usize);
            let b = _mm_loadu_si128(sse_b_ptr.offset(i));
            let eq = _mm_cmpeq_epi8(a, b);
            curr = _mm_sub_epi8(curr, eq); // subtract -1 = add 1 when matching
        }

        sad = _mm_add_epi64(sad, _mm_sad_epu8(curr, zeros));
        let mut sad_arr = [0u32; 4];
        _mm_storeu_si128(sad_arr.as_mut_ptr() as *mut __m128i, sad);
        let res = *sad_arr.get_unchecked(0) + *sad_arr.get_unchecked(2);

        len as u32 - res
    }
}
