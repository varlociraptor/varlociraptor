use core::cmp::Ordering;

use Iter;

pub trait Bits: Eq + PartialEq + PartialOrd + Ord + Copy {
    fn as_u64(self) -> u64;
}
impl Bits for u32 {
    fn as_u64(self) -> u64 { self as u64 }
}
impl Bits for u64 {
    fn as_u64(self) -> u64 { self }
}

/// Types that are IEEE754 floating point numbers.
pub trait Ieee754: Copy + PartialEq + PartialOrd {
    /// Iterate over each value of `Self` in `[self, lim]`.
    ///
    /// The returned iterator will include subnormal numbers, and will
    /// only include one of `-0.0` and `0.0`.
    ///
    /// # Panics
    ///
    /// Panics if `self > lim`, or if either are NaN.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use ieee754::Ieee754;
    ///
    /// // there are 840 single-precision floats in between 1.0 and 1.0001
    /// // (inclusive).
    /// assert_eq!(1_f32.upto(1.0001).count(), 840);
    /// ```
    fn upto(self, lim: Self) -> Iter<Self>;

    /// A type that represents the raw bits of `Self`.
    type Bits: Bits;
    /// A type large enough to store the true exponent of `Self`.
    type Exponent;
    /// A type large enough to store the raw exponent (i.e. with the bias).
    type RawExponent;
    /// A type large enough to store the significand of `Self`.
    type Significand;

    /// Return the next value after `self`.
    ///
    /// Calling this on NaN or positive infinity will yield nonsense.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use ieee754::Ieee754;
    /// let x: f32 = 1.0;
    /// assert_eq!(x.next(), 1.000000119209);
    /// ```
    fn next(self) -> Self;

    /// Return the previous value before `self`.
    ///
    /// Calling this on NaN or negative infinity will yield nonsense.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use ieee754::Ieee754;
    /// let x: f32 = 1.0;
    /// assert_eq!(x.prev(), 0.99999995);
    /// ```
    fn prev(self) -> Self;

    /// Return the unit-in-the-last-place ulp of `self`. That is,
    /// `x.abs().next() - x.abs()`, but handling overflow properly.
    ///
    /// Returns `None` if `self` is not finite.
    ///
    /// # Examples
    ///
    /// Single precision:
    ///
    /// ```rust
    /// use std::f32;
    /// use ieee754::Ieee754;
    ///
    /// assert_eq!(0_f32.ulp(), Some(1.4e-45));
    ///
    /// assert_eq!(1_f32.ulp(), Some(1.1920928955078125e-07));
    /// assert_eq!((-1_f32).ulp(), Some(1.1920928955078125e-07));
    ///
    /// // 2^23
    /// assert_eq!(8_388_608_f32.ulp(), Some(1.0));
    /// // 2^24 - 1, the largest f32 with ULP 1
    /// assert_eq!(16_777_215_f32.ulp(), Some(1.0));
    /// // 2^24
    /// assert_eq!(16_777_216_f32.ulp(), Some(2.0));
    ///
    /// // non-finite
    /// assert_eq!(f32::INFINITY.ulp(), None);
    /// assert_eq!(f32::NAN.ulp(), None);
    /// ```
    ///
    /// Double precision:
    ///
    /// ```rust
    /// use std::f64;
    /// use ieee754::Ieee754;
    ///
    /// assert_eq!(0_f64.ulp(), Some(4.9e-324));
    ///
    /// assert_eq!(1_f64.ulp(), Some(2.220446049250313e-16));
    /// assert_eq!((-1_f64).ulp(), Some(2.220446049250313e-16));
    ///
    /// // 2^52
    /// assert_eq!(4_503_599_627_370_496_f64.ulp(), Some(1.0));
    /// // 2^53 - 1, the largest f64 with ULP 1
    /// assert_eq!(9_007_199_254_740_991_f64.ulp(), Some(1.0));
    /// // 2^53
    /// assert_eq!(9_007_199_254_740_992_f64.ulp(), Some(2.0));
    ///
    /// // non-finite
    /// assert_eq!(f64::INFINITY.ulp(), None);
    /// assert_eq!(f64::NAN.ulp(), None);
    /// ```
    fn ulp(self) -> Option<Self>;

    /// View `self` as a collection of bits.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use ieee754::Ieee754;
    /// let x: f32 = 1.0;
    /// assert_eq!(x.bits(), 0x3f80_0000);
    /// ```
    fn bits(self) -> Self::Bits;

    /// View a collections of bits as a floating point number.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use ieee754::Ieee754;
    /// let float: f32 = Ieee754::from_bits(0xbf80_0000);
    /// assert_eq!(float, -1.0);
    /// ```
    fn from_bits(x: Self::Bits) -> Self;

    /// Get the bias of the stored exponent.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use ieee754::Ieee754;
    ///
    /// assert_eq!(f32::exponent_bias(), 127);
    /// assert_eq!(f64::exponent_bias(), 1023);
    /// ```
    fn exponent_bias() -> Self::Exponent;

    /// Break `self` into the three constituent parts of an IEEE754 float.
    ///
    /// The exponent returned is the raw bits, use `exponent_bias` to
    /// compute the offset required or use `decompose` to obtain this
    /// in precomputed form.
    ///
    /// # Examples
    ///
    /// Single precision:
    ///
    /// ```rust
    /// use ieee754::Ieee754;
    ///
    /// assert_eq!(1_f32.decompose_raw(), (false, 127, 0));
    /// assert_eq!(1234.567_f32.decompose_raw(), (false, 137, 0x1a5225));
    ///
    /// assert_eq!((-0.525_f32).decompose_raw(), (true, 126, 0x66666));
    ///
    /// assert_eq!(std::f32::INFINITY.decompose_raw(), (false, 255, 0));
    ///
    /// let (sign, expn, signif) = std::f32::NAN.decompose_raw();
    /// assert_eq!((sign, expn), (false, 255));
    /// assert!(signif != 0);
    /// ```
    ///
    /// Double precision:
    ///
    /// ```rust
    /// use ieee754::Ieee754;
    ///
    /// assert_eq!(1_f64.decompose_raw(), (false, 1023, 0));
    /// assert_eq!(1234.567_f64.decompose_raw(), (false, 1033, 0x34a449ba5e354));
    ///
    /// assert_eq!((-0.525_f64).decompose_raw(), (true, 1022, 0xcccc_cccc_cccd));
    ///
    /// assert_eq!(std::f64::INFINITY.decompose_raw(), (false, 2047, 0));
    ///
    /// let (sign, expn, signif) = std::f64::NAN.decompose_raw();
    /// assert_eq!((sign, expn), (false, 2047));
    /// assert!(signif != 0);
    /// ```
    fn decompose_raw(self) -> (bool, Self::RawExponent, Self::Significand);

    /// Create a `Self` out of the three constituent parts of an IEEE754 float.
    ///
    /// This returns (-1)<sup><code>sign</code></sup> ×
    /// 1.<code>signif</code> × 2<sup><code>expn</code> - bias</sup>, where
    ///
    /// - `sign` is treated as if `true` == `1` (meaning `true` is
    ///   negative),
    /// - 1.<code>signif</code> refers to placing the bits of `signif`
    ///   as the fractional part of a number between 1 and 2, and
    /// - bias is the exponent bias for this float (see [`exponent_bias`]).
    ///
    /// The exponent should be the raw bits: use `exponent_bias` to
    /// compute the offset required, or use `recompose` to feed in an
    /// unbiased exponent.
    ///
    /// # Examples
    ///
    /// Single precision:
    ///
    /// ```rust
    /// use ieee754::Ieee754;
    ///
    /// assert_eq!(f32::recompose_raw(false, 127, 0), 1.0);
    /// assert_eq!(f32::recompose_raw(false, 137, 0x1a5225), 1234.567);
    /// assert_eq!(f32::recompose_raw(true, 126, 0x66666), -0.525);
    ///
    /// assert_eq!(f32::recompose_raw(false, 255, 0), std::f32::INFINITY);
    ///
    /// assert!(f32::recompose_raw(false, 255, 1).is_nan());
    /// ```
    ///
    /// Double precision:
    ///
    /// ```rust
    /// use ieee754::Ieee754;
    ///
    /// assert_eq!(f64::recompose_raw(false, 1023, 0), 1.0);
    /// assert_eq!(f64::recompose_raw(false, 1033, 0x34a449ba5e354), 1234.567);
    /// assert_eq!(f64::recompose_raw(true, 1022, 0xcccc_cccc_cccd), -0.525);
    ///
    /// assert_eq!(f64::recompose_raw(false, 2047, 0), std::f64::INFINITY);
    ///
    /// assert!(f64::recompose_raw(false, 2047, 1).is_nan());
    /// ```
    fn recompose_raw(sign: bool, expn: Self::RawExponent, signif: Self::Significand) -> Self;

    /// Break `self` into the three constituent parts of an IEEE754 float.
    ///
    /// The exponent returned is the true exponent, after accounting
    /// for the bias it is stored with. The significand does not
    /// include the implicit highest bit (if it exists), e.g. the
    /// 24-bit for single precision.
    ///
    /// # Examples
    ///
    /// Single precision:
    ///
    /// ```rust
    /// use ieee754::Ieee754;
    ///
    /// assert_eq!(1_f32.decompose(), (false, 0, 0));
    /// assert_eq!(1234.567_f32.decompose(), (false, 10, 0x1a5225));
    ///
    /// assert_eq!((-0.525_f32).decompose(), (true, -1, 0x66666));
    ///
    /// assert_eq!(std::f32::INFINITY.decompose(), (false, 128, 0));
    ///
    /// let (sign, expn, signif) = std::f32::NAN.decompose();
    /// assert_eq!((sign, expn), (false, 128));
    /// assert!(signif != 0);
    /// ```
    ///
    /// Double precision:
    ///
    /// ```rust
    /// use ieee754::Ieee754;
    ///
    /// assert_eq!(1_f64.decompose(), (false, 0, 0));
    /// assert_eq!(1234.567_f64.decompose(), (false, 10, 0x34a449ba5e354));
    ///
    /// assert_eq!((-0.525_f64).decompose(), (true, -1, 0xcccc_cccc_cccd));
    ///
    /// assert_eq!(std::f64::INFINITY.decompose(), (false, 1024, 0));
    ///
    /// let (sign, expn, signif) = std::f64::NAN.decompose();
    /// assert_eq!((sign, expn), (false, 1024));
    /// assert!(signif != 0);
    /// ```
    fn decompose(self) -> (bool, Self::Exponent, Self::Significand);

    /// Create a `Self` out of the three constituent parts of an IEEE754 float.
    ///
    /// This returns (-1)<sup><code>sign</code></sup> ×
    /// 1.<code>signif</code> × 2<sup><code>expn</code></sup>, where
    ///
    /// - `sign` is treated as if `true` == `1` (meaning `true` is
    ///   negative), and
    /// - 1.<code>signif</code> refers to placing the bits of `signif`
    ///   as the fractional part of a number between 1 and 2.
    ///
    /// The exponent should be the true exponent, not accounting for any
    /// bias. The significand should not include the implicit highest
    /// bit (if it exists), e.g. the 24-th bit for single precision.
    ///
    /// # Examples
    ///
    /// Single precision:
    ///
    /// ```rust
    /// use ieee754::Ieee754;
    ///
    /// // normal numbers
    /// assert_eq!(f32::recompose(false, 0, 0), 1.0);
    /// assert_eq!(f32::recompose(false, 10, 0x1a5225), 1234.567);
    /// assert_eq!(f32::recompose(true, -1, 0x66666), -0.525);
    ///
    /// // infinity
    /// assert_eq!(f32::recompose(false, 128, 0), std::f32::INFINITY);
    ///
    /// // NaN
    /// assert!(f32::recompose(false, 128, 1).is_nan());
    /// ```
    ///
    /// Double precision:
    ///
    /// ```rust
    /// use ieee754::Ieee754;
    ///
    /// // normal numbers
    /// assert_eq!(f64::recompose(false, 0, 0), 1.0);
    /// assert_eq!(f64::recompose(false, 10, 0x34a449ba5e354), 1234.567);
    /// assert_eq!(f64::recompose(true, -1, 0xcccc_cccc_cccd), -0.525);
    ///
    /// // infinity
    /// assert_eq!(f64::recompose(false, 1024, 0), std::f64::INFINITY);
    ///
    /// // NaN
    /// assert!(f64::recompose(false, 1024, 1).is_nan());
    /// ```
    fn recompose(sign: bool, expn: Self::Exponent, signif: Self::Significand) -> Self;

    /// Compare `x` and `y` using the IEEE-754 `totalOrder` predicate
    /// (Section 5.10).
    ///
    /// This orders NaNs before or after all non-NaN floats, depending
    /// on the sign bit. Using -qNaN to represent a quiet NaN with
    /// negative sign bit and similarly for a signalling NaN (sNaN),
    /// the order is:
    ///
    /// ```txt
    /// -qNaN < -sNaN < -∞ < -12.34 < -0.0 < +0.0 < +12.34 < +∞ < +sNaN < +qNaN
    /// ```
    ///
    /// (NaNs are ordered according to their payload.)
    ///
    /// # Examples
    ///
    /// Sorting:
    ///
    /// ```rust
    /// use std::f32;
    ///
    /// use ieee754::Ieee754;
    ///
    /// let mut data = vec![0.0, f32::NEG_INFINITY, -1.0, f32::INFINITY,
    ///                     f32::NAN, -0.0, 12.34e5, -f32::NAN];
    /// data.sort_by(|a, b| a.total_cmp(b));
    ///
    /// assert_eq!(format!("{:.0?}", data),
    ///            "[NaN, -inf, -1, -0, 0, 1234000, inf, NaN]");
    /// ```
    ///
    /// Single precision:
    ///
    /// ```rust
    /// use std::cmp::Ordering;
    /// use std::f32;
    ///
    /// use ieee754::Ieee754;
    ///
    /// // normal comparison
    /// assert_eq!(0_f32.total_cmp(&0_f32), Ordering::Equal);
    /// assert_eq!(0_f32.total_cmp(&1_f32), Ordering::Less);
    /// assert_eq!(1e10_f32.total_cmp(&f32::NEG_INFINITY), Ordering::Greater);
    ///
    /// // signed zero
    /// assert_eq!(0_f32.total_cmp(&-0_f32), Ordering::Greater);
    ///
    /// // NaNs
    /// assert_eq!(f32::NAN.total_cmp(&0_f32), Ordering::Greater);
    /// assert_eq!(f32::NAN.total_cmp(&f32::INFINITY), Ordering::Greater);
    /// assert_eq!((-f32::NAN).total_cmp(&f32::NEG_INFINITY), Ordering::Less);
    /// ```
    ///
    /// Double precision:
    ///
    /// ```rust
    /// use std::cmp::Ordering;
    /// use std::f64;
    ///
    /// use ieee754::Ieee754;
    ///
    /// // normal comparison
    /// assert_eq!(0_f64.total_cmp(&0_f64), Ordering::Equal);
    /// assert_eq!(0_f64.total_cmp(&1_f64), Ordering::Less);
    /// assert_eq!(1e10_f64.total_cmp(&f64::NEG_INFINITY), Ordering::Greater);
    ///
    /// // signed zero
    /// assert_eq!(0_f64.total_cmp(&-0_f64), Ordering::Greater);
    ///
    /// // NaNs
    /// assert_eq!(f64::NAN.total_cmp(&0_f64), Ordering::Greater);
    /// assert_eq!(f64::NAN.total_cmp(&f64::INFINITY), Ordering::Greater);
    /// assert_eq!((-f64::NAN).total_cmp(&f64::NEG_INFINITY), Ordering::Less);
    /// ```
    fn total_cmp(&self, other: &Self) -> Ordering;

    /// Return the absolute value of `x`.
    ///
    /// This provides a no_std/core-only version of the built-in `abs` in
    /// `std`, until
    /// [#50145](https://github.com/rust-lang/rust/issues/50145) is
    /// addressed.
    ///
    /// # Examples
    ///
    /// Single precision:
    ///
    /// ```rust
    /// #![no_std]
    /// # extern crate std; // this makes this "test" a lie, unfortunately
    /// # extern crate ieee754;
    /// use core::f32;
    ///
    /// use ieee754::Ieee754;
    ///
    /// # fn main() {
    /// assert_eq!((0_f32).abs(), 0.0);
    ///
    /// assert_eq!((12.34_f32).abs(), 12.34);
    /// assert_eq!((-12.34_f32).abs(), 12.34);
    ///
    /// assert_eq!(f32::INFINITY.abs(), f32::INFINITY);
    /// assert_eq!(f32::NEG_INFINITY.abs(), f32::INFINITY);
    /// assert!(f32::NAN.abs().is_nan());
    /// # }
    /// ```
    ///
    /// Double precision:
    ///
    /// ```rust
    /// #![no_std]
    /// # extern crate std; // this makes this "test" a lie, unfortunately
    /// # extern crate ieee754;
    /// use core::f64;
    ///
    /// use ieee754::Ieee754;
    ///
    /// # fn main() {
    /// assert_eq!((0_f64).abs(), 0.0);
    ///
    /// assert_eq!((12.34_f64).abs(), 12.34);
    /// assert_eq!((-12.34_f64).abs(), 12.34);
    ///
    /// assert_eq!(f64::INFINITY.abs(), f64::INFINITY);
    /// assert_eq!(f64::NEG_INFINITY.abs(), f64::INFINITY);
    /// assert!(f64::NAN.abs().is_nan());
    /// # }
    /// ```
    fn abs(self) -> Self;

    /// Return a float with the magnitude of `self` but the sign of
    /// `sign`.
    ///
    /// If `sign` is NaN, this still uses its sign bit, and does not
    /// (necessarily) return NaN.
    ///
    /// # Examples
    ///
    /// Single precision:
    ///
    /// ```rust
    /// use std::f32;
    ///
    /// use ieee754::Ieee754;
    ///
    /// // normal numbers
    /// assert_eq!(1_f32.copy_sign(1.0), 1.0);
    /// assert_eq!(2_f32.copy_sign(-1.0), -2.0);
    /// assert_eq!((-3_f32).copy_sign(1.0), 3.0);
    /// assert_eq!((-4_f32).copy_sign(-1.0), -4.0);
    ///
    /// // infinities
    /// assert_eq!(5_f32.copy_sign(f32::NEG_INFINITY), -5.0);
    /// assert_eq!(f32::NEG_INFINITY.copy_sign(1.0), f32::INFINITY);
    ///
    /// // signs of zeros matter
    /// assert_eq!((-6_f32).copy_sign(0.0), 6.0);
    /// assert_eq!(7_f32.copy_sign(-0.0), -7.0);
    ///
    /// // NaNs only propagate on the self argument
    /// assert!(f32::NAN.copy_sign(1.0).is_nan());
    /// assert_eq!(8_f32.copy_sign(-f32::NAN), -8.0);
    /// ```
    ///
    /// Double precision:
    ///
    /// ```rust
    /// use std::f64;
    ///
    /// use ieee754::Ieee754;
    ///
    /// // normal numbers
    /// assert_eq!(1_f64.copy_sign(1.0), 1.0);
    /// assert_eq!(2_f64.copy_sign(-1.0), -2.0);
    /// assert_eq!((-3_f64).copy_sign(1.0), 3.0);
    /// assert_eq!((-4_f64).copy_sign(-1.0), -4.0);
    ///
    /// // infinities
    /// assert_eq!(5_f64.copy_sign(f64::NEG_INFINITY), -5.0);
    /// assert_eq!(f64::NEG_INFINITY.copy_sign(1.0), f64::INFINITY);
    ///
    /// // signs of zeros matter
    /// assert_eq!((-6_f64).copy_sign(0.0), 6.0);
    /// assert_eq!(7_f64.copy_sign(-0.0), -7.0);
    ///
    /// // NaNs only propagate on the self argument
    /// assert!(f64::NAN.copy_sign(1.0).is_nan());
    /// assert_eq!(8_f64.copy_sign(-f64::NAN), -8.0);
    /// ```
    fn copy_sign(self, sign: Self) -> Self;

    /// Return the sign of `x`.
    ///
    /// This provides a no_std/core-only function similar to the
    /// built-in `signum` in `std` (until
    /// [#50145](https://github.com/rust-lang/rust/issues/50145) is
    /// addressed). This `sign` function differs at two values; it
    /// matches the mathematical definitions when `self == 0.0` :
    ///
    /// | `x` | `x.signum()` (`std`) | `x.sign()` (`ieee754`) |
    /// |--:|--:|--:|
    /// |< 0.0|−1.0|−1.0|
    /// |−0.0|−1.0|**−0.0**|
    /// |+0.0|+1.0|**+0.0**|
    /// |> 0.0|+1.0|+1.0|
    /// |NaN|NaN|NaN|
    ///
    /// # Examples
    ///
    /// Single precision:
    ///
    /// ```rust
    /// use std::f32;
    /// use std::cmp::Ordering;
    ///
    /// use ieee754::Ieee754;
    ///
    /// // zeros
    /// assert_eq!(0_f32.sign().total_cmp(&0.0), Ordering::Equal);
    /// assert_eq!((-0_f32).sign().total_cmp(&-0.0), Ordering::Equal);
    ///
    /// // normal numbers
    /// assert_eq!((12.34_f32).sign(), 1.0);
    /// assert_eq!((-12.34_f32).sign(), -1.0);
    ///
    /// // extremes
    /// assert_eq!(f32::INFINITY.sign(), 1.0);
    /// assert_eq!(f32::NEG_INFINITY.sign(), -1.0);
    /// assert!(f32::NAN.sign().is_nan());
    /// ```
    ///
    /// Double precision:
    ///
    /// ```rust
    /// use std::f64;
    /// use std::cmp::Ordering;
    ///
    /// use ieee754::Ieee754;
    ///
    /// // zeros
    /// assert_eq!(0_f64.sign().total_cmp(&0.0), Ordering::Equal);
    /// assert_eq!((-0_f64).sign().total_cmp(&-0.0), Ordering::Equal);
    ///
    /// // normal numbers
    /// assert_eq!((12.34_f64).sign(), 1.0);
    /// assert_eq!((-12.34_f64).sign(), -1.0);
    ///
    /// // extremes
    /// assert_eq!(f64::INFINITY.sign(), 1.0);
    /// assert_eq!(f64::NEG_INFINITY.sign(), -1.0);
    /// assert!(f64::NAN.sign().is_nan());
    /// ```
    fn sign(self) -> Self;

    /// Compute the (generalized) **signed** relative error of `self`
    /// as an approximation to `exact`.
    ///
    /// This computes the signed value: positive indicates `self` in
    /// the opposite direction to 0 from `exact`; negative indicates
    /// `self` is in the same direction as 0 from `exact`. Use
    /// `x.rel_error(exact).abs()` to get the non-signed relative
    /// error.
    ///
    /// The "generalized" refers to `exact` being 0 or ±∞ the handling
    /// of which is designed to indicate a "failure" (infinite error),
    /// if `self` doesn't precisely equal `exact`. This behaviour is
    /// designed for checking output of algorithms on floats when it
    /// is often desirable to match 0.0 and ±∞ perfectly.
    ///
    /// The values of this function are:
    ///
    /// |`exact`|`x`|`x.rel_error(exact)`|
    /// |--:|--:|--:|
    /// |NaN|any value|NaN|
    /// |any value|NaN|NaN|
    /// |0|equal to `exact`|0|
    /// |0|not equal to `exact`|signum(`x`) × ∞|
    /// |±∞|equal to `exact`|0|
    /// |±∞|not equal to `exact`|-∞|
    /// |any other value|any value|`(x - exact) / exact`|
    ///
    /// The sign of a zero-valued argument has no effect on the result
    /// of this function.
    ///
    /// # Examples
    ///
    /// Single precision:
    ///
    /// ```rust
    /// use std::f32;
    ///
    /// use ieee754::Ieee754;
    ///
    /// assert_eq!(4_f32.rel_error(4.0), 0.0);
    /// assert_eq!(3_f32.rel_error(4.0), -0.25);
    /// assert_eq!(5_f32.rel_error(4.0), 0.25);
    ///
    /// // zero
    /// assert_eq!(0_f32.rel_error(0.0), 0.0);
    /// assert_eq!(1_f32.rel_error(0.0), f32::INFINITY);
    /// assert_eq!((-1_f32).rel_error(0.0), f32::NEG_INFINITY);
    ///
    /// // infinities
    /// assert_eq!(f32::INFINITY.rel_error(f32::INFINITY), 0.0);
    /// assert_eq!(0_f32.rel_error(f32::INFINITY), f32::NEG_INFINITY);
    ///
    /// assert_eq!(f32::NEG_INFINITY.rel_error(f32::NEG_INFINITY), 0.0);
    /// assert_eq!(0_f32.rel_error(f32::NEG_INFINITY), f32::NEG_INFINITY);
    ///
    /// // NaNs
    /// assert!(f32::NAN.rel_error(4.0).is_nan());
    /// assert!(4_f32.rel_error(f32::NAN).is_nan());
    /// ```
    ///
    /// Double precision:
    ///
    /// ```rust
    /// use std::f64;
    /// use ieee754::Ieee754;
    ///
    /// assert_eq!(4_f64.rel_error(4.0), 0.0);
    /// assert_eq!(3_f64.rel_error(4.0), -0.25);
    /// assert_eq!(5_f64.rel_error(4.0), 0.25);
    ///
    /// // zero
    /// assert_eq!(0_f64.rel_error(0.0), 0.0);
    /// assert_eq!(1_f64.rel_error(0.0), f64::INFINITY);
    /// assert_eq!((-1_f64).rel_error(0.0), f64::NEG_INFINITY);
    ///
    /// // infinities
    /// assert_eq!(f64::INFINITY.rel_error(f64::INFINITY), 0.0);
    /// assert_eq!(0_f64.rel_error(f64::INFINITY), f64::NEG_INFINITY);
    ///
    /// assert_eq!(f64::NEG_INFINITY.rel_error(f64::NEG_INFINITY), 0.0);
    /// assert_eq!(0_f64.rel_error(f64::NEG_INFINITY), f64::NEG_INFINITY);
    ///
    /// // NaNs
    /// assert!(f64::NAN.rel_error(4.0).is_nan());
    /// assert!(4_f64.rel_error(f64::NAN).is_nan());
    /// ```
    fn rel_error(self, exact: Self) -> Self;
}
