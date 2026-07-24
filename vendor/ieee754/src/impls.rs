use core::mem;
use core::cmp::Ordering;
use {Iter, Ieee754};

macro_rules! mask{
    ($bits: expr; $current: expr => $($other: expr),*) => {
        ($bits >> (0 $(+ $other)*)) & ((1 << $current) - 1)
    }
}
macro_rules! unmask {
    ($x: expr => $($other: expr),*) => {
        $x << (0 $(+ $other)*)
    }
}

macro_rules! mk_impl {
    ($f: ident, $bits: ty, $signed_bits: ty,
     $expn: ty, $expn_raw: ty, $signif: ty,
     $expn_n: expr, $signif_n: expr) => {
        impl Ieee754 for $f {
            type Bits = $bits;
            type Exponent = $expn;
            type RawExponent = $expn_raw;
            type Significand = $signif;
            #[inline]
            fn upto(self, lim: Self) -> Iter<Self> {
                assert!(self <= lim);
                // map -0.0 to 0.0, i.e. ensure that any zero is
                // stored in a canonical manner. This is necessary to
                // use bit-hacks for the comparison in next.
                #[inline(always)]
                fn canon(x: $f) -> $f { if x == 0.0 { 0.0 } else { x } }

                ::iter::new_iter(canon(self), canon(lim))
            }
            #[inline]
            fn ulp(self) -> Option<Self> {
                let (_sign, expn, _signif) = self.decompose_raw();

                const MAX_EXPN: $expn_raw = ((1u64 << $expn_n) - 1) as $expn_raw;
                match expn {
                    MAX_EXPN => None,
                    0 => Some($f::recompose_raw(false, 0, 1)),
                    _ => {
                        let ulp_expn = expn.saturating_sub($signif_n);
                        if ulp_expn == 0 {
                            Some($f::recompose_raw(false, 0, 1 << (expn - 1)))
                        } else {
                            Some($f::recompose_raw(false, ulp_expn, 0))
                        }
                    }
                }
            }

            #[inline]
            fn next(self) -> Self {
                let abs_mask = (!(0 as Self::Bits)) >> 1;
                let mut bits = self.bits();
                if self == 0.0 {
                    bits = 1;
                } else if self < 0.0 {
                    bits -= 1;
                    if bits == !abs_mask {
                        // normalise -0.0 to +0.0
                        bits = 0
                    }
                } else {
                    bits += 1
                }
                Ieee754::from_bits(bits)
            }
            #[inline]
            fn prev(self) -> Self {
                let abs_mask = (!(0 as Self::Bits)) >> 1;
                let mut bits = self.bits();
                if self < 0.0 {
                     bits += 1;
                } else if bits & abs_mask == 0 {
                     bits = 1 | !abs_mask;
                } else {
                     bits -= 1;
                }
                Ieee754::from_bits(bits)
            }

            #[inline]
            fn exponent_bias() -> Self::Exponent {
                (1 << ($expn_n - 1)) - 1
            }

            #[inline]
            fn bits(self) -> Self::Bits {
                unsafe {mem::transmute(self)}
            }
            #[inline]
            fn from_bits(bits: Self::Bits) -> Self {
                unsafe {mem::transmute(bits)}
            }
            #[inline]
            fn decompose_raw(self) -> (bool, Self::RawExponent, Self::Significand) {
                let bits = self.bits();

                (mask!(bits; 1 => $expn_n, $signif_n) != 0,
                 mask!(bits; $expn_n => $signif_n) as Self::RawExponent,
                 mask!(bits; $signif_n => ) as Self::Significand)

            }
            #[inline]
            fn recompose_raw(sign: bool, expn: Self::RawExponent, signif: Self::Significand) -> Self {
                Ieee754::from_bits(
                    unmask!(sign as Self::Bits => $expn_n, $signif_n) |
                    unmask!(expn as Self::Bits => $signif_n) |
                    unmask!(signif as Self::Bits => ))
            }

            #[inline]
            fn decompose(self) -> (bool, Self::Exponent, Self::Significand) {
                let (sign, expn, signif) = self.decompose_raw();
                (sign, expn as Self::Exponent - Self::exponent_bias(),
                 signif)
            }
            #[inline]
            fn recompose(sign: bool, expn: Self::Exponent, signif: Self::Significand) -> Self {
                Self::recompose_raw(sign,
                                    (expn + Self::exponent_bias()) as Self::RawExponent,
                                    signif)
            }

            #[inline]
            fn total_cmp(&self, other: &Self) -> Ordering {
                #[inline]
                fn cmp_key(x: $f) -> $signed_bits {
                    let bits = x.bits();
                    let sign_bit = bits & (1 << ($expn_n + $signif_n));
                    let mask = ((sign_bit as $signed_bits) >> ($expn_n + $signif_n)) as $bits >> 1;
                    (bits ^ mask) as $signed_bits
                }
                cmp_key(*self).cmp(&cmp_key(*other))
            }

            #[inline]
            fn abs(self) -> Self {
                let (_, e, s) = self.decompose_raw();
                Self::recompose_raw(false, e, s)
            }

            #[inline]
            fn copy_sign(self, sign: Self) -> Self {
                let (s, _, _) = sign.decompose_raw();
                let (_, e, m) = self.decompose_raw();
                Self::recompose_raw(s, e, m)
            }

            #[inline]
            fn sign(self) -> Self {
                if self == 0.0 || self != self { // is NaN? (.is_nan() added to core in 1.27)
                    self
                } else {
                    1.0.copy_sign(self)
                }
            }

            #[inline]
            fn rel_error(self, exact: Self) -> Self {
                use core::$f;
                if exact == 0.0 {
                    if self == 0.0 {
                        0.0
                    } else {
                        self.sign() * $f::INFINITY
                    }
                } else if exact.abs() == $f::INFINITY {
                    if self == exact {
                        0.0
                    } else if self != self { // is NaN? (.is_nan() added to core in 1.27)
                        self
                    } else {
                        $f::NEG_INFINITY
                    }
                } else {
                    // NaNs propagate
                    (self - exact) / exact
                }
            }
        }

        #[cfg(test)]
        mod $f {
            use std::prelude::v1::*;
            use std::$f;

            use {Ieee754};
            #[test]
            fn upto() {
                assert_eq!((0.0 as $f).upto(0.0).collect::<Vec<_>>(),
                           &[0.0]);
                assert_eq!($f::recompose(false, 1, 1).upto($f::recompose(false, 1, 10)).count(),
                           10);

                assert_eq!($f::recompose(true, -$f::exponent_bias(), 10)
                           .upto($f::recompose(false, -$f::exponent_bias(), 10)).count(),
                           21);
            }
            #[test]
            fn upto_rev() {
                assert_eq!(0.0_f32.upto(0.0_f32).rev().collect::<Vec<_>>(),
                           &[0.0]);

                assert_eq!($f::recompose(false, 1, 1)
                           .upto($f::recompose(false, 1, 10)).rev().count(),
                           10);
                assert_eq!($f::recompose(true, -$f::exponent_bias(), 10)
                           .upto($f::recompose(false, -$f::exponent_bias(), 10)).rev().count(),
                           21);
            }

            #[test]
            fn upto_infinities() {
                use std::$f as f;
                assert_eq!(f::MAX.upto(f::INFINITY).collect::<Vec<_>>(),
                           &[f::MAX, f::INFINITY]);
                assert_eq!(f::NEG_INFINITY.upto(f::MIN).collect::<Vec<_>>(),
                           &[f::NEG_INFINITY, f::MIN]);
            }
            #[test]
            fn upto_infinities_rev() {
                use std::$f as f;
                assert_eq!(f::MAX.upto(f::INFINITY).rev().collect::<Vec<_>>(),
                           &[f::INFINITY, f::MAX]);
                assert_eq!(f::NEG_INFINITY.upto(f::MIN).rev().collect::<Vec<_>>(),
                           &[f::MIN, f::NEG_INFINITY]);
            }

            #[test]
            fn upto_size_hint() {
                let mut iter =
                    $f::recompose(true, -$f::exponent_bias(), 10)
                    .upto($f::recompose(false, -$f::exponent_bias(), 10));

                assert_eq!(iter.size_hint(), (21, Some(21)));
                for i in (0..21).rev() {
                    assert!(iter.next().is_some());
                    assert_eq!(iter.size_hint(), (i, Some(i)));
                }
                assert_eq!(iter.next(), None);
                assert_eq!(iter.size_hint(), (0, Some(0)))
            }

            #[test]
            fn upto_size_hint_rev() {
                let mut iter =
                    $f::recompose(true, -$f::exponent_bias(), 10)
                    .upto($f::recompose(false, -$f::exponent_bias(), 10))
                    .rev();

                assert_eq!(iter.size_hint(), (21, Some(21)));
                for i in (0..21).rev() {
                    assert!(iter.next().is_some());
                    assert_eq!(iter.size_hint(), (i, Some(i)));
                }
                assert_eq!(iter.next(), None);
                assert_eq!(iter.size_hint(), (0, Some(0)))
            }

            #[test]
            fn next_prev_order() {
                let cases = [0.0 as $f, -0.0, 1.0, 1.0001, 1e30, -1.0, -1.0001, -1e30];
                for &x in &cases {
                    assert!(x.next() > x);
                    assert!(x.prev() < x);
                }
            }

            #[test]
            fn ulp_smoke() {
                let smallest_subnormal = $f::recompose_raw(false, 0, 1);
                let smallest_normal = $f::recompose_raw(false, 1, 0);
                assert_eq!((0.0 as $f).ulp(), Some(smallest_subnormal));
                assert_eq!(smallest_subnormal.ulp(), Some(smallest_subnormal));
                assert_eq!($f::recompose_raw(true, 0, 9436).ulp(),
                           Some(smallest_subnormal));
                assert_eq!(smallest_normal.ulp(), Some(smallest_subnormal));

                assert_eq!((1.0 as $f).ulp(),
                           Some($f::recompose(false, -$signif_n, 0)));

                assert_eq!((-123.456e30 as $f).ulp(),
                           Some($f::recompose(false, 106 - $signif_n, 0)));

                assert_eq!($f::INFINITY.ulp(), None);
                assert_eq!($f::NEG_INFINITY.ulp(), None);
                assert_eq!($f::NAN.ulp(), None);
            }

            #[test]
            fn ulp_aggressive() {
                fn check_ulp(x: $f, ulp: $f) {
                    println!("  {:e} {:e}", x, ulp);
                    assert_eq!(x.ulp(), Some(ulp));
                    // with signed-magnitude we need to be moving away
                    let same_sign_ulp = if x < 0.0 { -ulp } else { ulp };

                    assert_ne!(x + same_sign_ulp, x, "adding ulp should be different");

                    if ulp / 2.0 > 0.0 {
                        // floats break ties like this by rounding to
                        // even (in the default mode), so adding half
                        // a ulp may be a new value depending on the
                        // significand.
                        if x.decompose().2 & 1 == 0 {
                            assert_eq!(x + same_sign_ulp / 2.0, x);
                        } else {
                            assert_eq!(x + same_sign_ulp / 2.0, x + same_sign_ulp);
                        }
                    }
                    // no ties to worry about
                    assert_eq!(x + same_sign_ulp / 4.0, x);
                }

                let smallest_subnormal = $f::recompose_raw(false, 0, 1);
                let mut ulp = smallest_subnormal;

                check_ulp(0.0, ulp);

                let mut pow2 = smallest_subnormal;
                for i in 0..200 {
                    println!("{}", i);
                    check_ulp(pow2, ulp);
                    check_ulp(-pow2, ulp);

                    let (_, e, _) = pow2.decompose_raw();
                    if e > 0 {
                        for &signif in &[1,
                                         // random numbers
                                         9436, 1577069,
                                         // last two for this exponent
                                         (1 << $signif_n) - 2, (1 << $signif_n) - 1] {
                            check_ulp($f::recompose_raw(false, e, signif), ulp);
                            check_ulp($f::recompose_raw(true, e, signif), ulp);
                        }
                    }

                    pow2 *= 2.0;
                    if i >= $signif_n {
                        ulp *= 2.0;
                    }
                }
            }

            #[test]
            fn abs() {
                let this_abs = <$f as Ieee754>::abs;
                assert!(this_abs($f::NAN).is_nan());

                let cases = [0.0 as $f, -1.0, 1.0001,
                             // denormals
                             $f::recompose_raw(false, 0, 123), $f::recompose(true, 0, 123),
                             $f::NEG_INFINITY, $f::INFINITY];
                for x in &cases {
                    assert_eq!(this_abs(*x), x.abs());
                }
            }

            #[test]
            fn total_cmp() {
                let nan_exp = $f::NAN.decompose_raw().1;
                let q = 1 << ($signif_n - 1);

                let qnan0 = $f::recompose_raw(false, nan_exp, q);
                let qnan1 = $f::recompose_raw(false, nan_exp, q | 1);
                let qnanlarge = $f::recompose_raw(false, nan_exp, q | (q - 1));

                let snan1 = $f::recompose_raw(false, nan_exp, 1);
                let snan2 = $f::recompose_raw(false, nan_exp, 2);
                let snanlarge = $f::recompose_raw(false, nan_exp, q - 1);

                let subnormal = $f::recompose_raw(false, 0, 1);

                // it's a total order, so we can literally write
                // options in order, and compare them all, using their
                // indices as ground-truth. NB. the snans seem to
                // get canonicalized to qnan on some versions of i686
                // Linux (using `cross` on Travis CI), so we can't
                // include them.
                let include_snan = cfg!(not(target_arch = "x86"));
                // -qNaN
                let mut cases = vec![-qnanlarge, -qnan1, -qnan0];
                // -sNaN
                if include_snan {
                    cases.extend_from_slice(&[-snanlarge, -snan2, -snan1]);
                }
                // Numbers (note -0, +0)
                cases.extend_from_slice(&[
                    $f::NEG_INFINITY,
                    -1e15, -1.001, -1.0, -0.999, -1e-15, -subnormal,
                    -0.0, 0.0,
                    subnormal, 1e-15, 0.999, 1.0, 1.001, 1e15,
                    $f::INFINITY
                ]);
                // +sNaN
                if include_snan {
                    cases.extend_from_slice(&[snan1, snan2, snanlarge]);
                }
                // +qNaN
                cases.extend_from_slice(&[qnan0, qnan1, qnanlarge]);

                for (ix, &x) in cases.iter().enumerate() {
                    for (iy, &y) in cases.iter().enumerate() {
                        let computed = x.total_cmp(&y);
                        let expected = ix.cmp(&iy);
                        assert_eq!(
                            computed, expected,
                            "{:e} ({}, {:?}) cmp {:e} ({}, {:?}), got: {:?}, expected: {:?}",
                            x, ix, x.decompose(),
                            y, iy, y.decompose(),
                            computed, expected);
                    }
                }
            }

            #[test]
            fn copy_sign() {
                let positives = [$f::NAN, $f::INFINITY, 1e10, 1.0, 0.0];

                for &pos in &positives {
                    assert_eq!((1 as $f).copy_sign(pos), 1.0);
                    assert_eq!((1 as $f).copy_sign(-pos), -1.0);
                    assert_eq!((-1 as $f).copy_sign(pos), 1.0);
                    assert_eq!((-1 as $f).copy_sign(-pos), -1.0);

                    assert_eq!($f::INFINITY.copy_sign(pos), $f::INFINITY);
                    assert_eq!($f::INFINITY.copy_sign(-pos), $f::NEG_INFINITY);
                    assert_eq!($f::NEG_INFINITY.copy_sign(pos), $f::INFINITY);
                    assert_eq!($f::NEG_INFINITY.copy_sign(-pos), $f::NEG_INFINITY);

                    assert!($f::NAN.copy_sign(pos).is_nan());
                    assert!($f::NAN.copy_sign(-pos).is_nan());
                }
            }

            #[test]
            fn sign() {
                assert!($f::NAN.sign().is_nan());

                assert_eq!($f::NEG_INFINITY.sign(), -1.0);
                assert_eq!((-1e10 as $f).sign(), -1.0);
                assert_eq!((-1.0 as $f).sign(), -1.0);
                assert_eq!((-0.0 as $f).sign().bits(), (-0.0 as $f).bits());
                assert_eq!((0.0 as $f).sign().bits(), (0.0 as $f).bits());
                assert_eq!((1.0 as $f).sign(), 1.0);
                assert_eq!((1e10 as $f).sign(), 1.0);
                assert_eq!($f::INFINITY.sign(), 1.0);
            }

            #[test]
            fn rel_error() {
                let zer: $f = 0.0;
                let one: $f = 1.0;
                let two: $f = 2.0;

                assert_eq!(zer.rel_error(one), -1.0);
                assert_eq!(one.rel_error(one), 0.0);
                assert_eq!((-one).rel_error(one), -2.0);
                assert_eq!(two.rel_error(one), 1.0);
                assert_eq!((-two).rel_error(one), -3.0);

                assert_eq!(zer.rel_error(-one), -1.0);
                assert_eq!(one.rel_error(-one), -2.0);
                assert_eq!((-one).rel_error(-one), 0.0);
                assert_eq!(two.rel_error(-one), -3.0);
                assert_eq!((-two).rel_error(-one), 1.0);

                assert_eq!(zer.rel_error(two), -1.0);
                assert_eq!(one.rel_error(two), -0.5);
                assert_eq!((-one).rel_error(two), -1.5);
                assert_eq!(two.rel_error(two), 0.0);
                assert_eq!((-two).rel_error(two), -2.0);

                assert_eq!(zer.rel_error(-two), -1.0);
                assert_eq!(one.rel_error(-two), -1.5);
                assert_eq!((-one).rel_error(-two), -0.5);
                assert_eq!(two.rel_error(-two), -2.0);
                assert_eq!((-two).rel_error(-two), 0.0);
            }
            #[test]
            fn rel_error_edge_cases() {
                let nan = $f::NAN;
                let inf = $f::INFINITY;
                let zer: $f = 0.0;
                let one: $f = 1.0;

                assert!(nan.rel_error(nan).is_nan());
                assert!(zer.rel_error(nan).is_nan());
                assert!(nan.rel_error(zer).is_nan());

                assert_eq!(zer.rel_error(zer), 0.0);
                assert_eq!(zer.rel_error(-zer), 0.0);
                assert_eq!((-zer).rel_error(zer), 0.0);
                assert_eq!((-zer).rel_error(-zer), 0.0);
                assert_eq!(one.rel_error(zer), inf);
                assert_eq!((-one).rel_error(zer), -inf);
                assert_eq!(inf.rel_error(zer), inf);
                assert_eq!((-inf).rel_error(zer), -inf);


                assert_eq!(inf.rel_error(inf), 0.0);
                assert_eq!(inf.rel_error(-inf), -inf);
                assert_eq!((-inf).rel_error(inf), -inf);
                assert_eq!((-inf).rel_error(-inf), 0.0);
                assert_eq!(zer.rel_error(inf), -inf);
                assert_eq!(zer.rel_error(-inf), -inf);
            }
        }
    }
}

mk_impl!(f32, u32, i32, i16, u8, u32, 8, 23);
mk_impl!(f64, u64, i64, i16, u16, u64, 11, 52);
