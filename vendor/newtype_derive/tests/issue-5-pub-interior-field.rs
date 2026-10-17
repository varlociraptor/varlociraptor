/*
Copyright ⓒ 2015 rust-custom-derive contributors.

Licensed under the MIT license (see LICENSE or <http://opensource.org
/licenses/MIT>) or the Apache License, Version 2.0 (see LICENSE of
<http://www.apache.org/licenses/LICENSE-2.0>), at your option. All
files in the project carrying such notice may not be copied, modified,
or distributed except according to those terms.
*/
#![recursion_limit = "128"]
#![cfg_attr(feature = "std-unstable", feature(zero_one))]
#![allow(deprecated)]
#[macro_use] extern crate custom_derive;
#[macro_use] extern crate newtype_derive;

use std::fmt::{self, Binary, Debug, Display, LowerExp, LowerHex, Octal, Pointer,
    UpperExp, UpperHex};

macro_rules! impl_fmt {
    (impl $tr:ident for $name:ident: $msg:expr) => {
        impl $tr for $name {
            fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
                write!(fmt, $msg)
            }
        }
    };
}

custom_derive! {
    #[derive(Copy, Clone, Eq, PartialEq, Debug,
        NewtypeAdd, NewtypeAdd(&self), NewtypeAdd(i32), NewtypeAdd(&self, i32),
        NewtypeBitAnd, NewtypeBitAnd(&self),
        NewtypeBitOr, NewtypeBitOr(&self),
        NewtypeBitXor, NewtypeBitXor(&self),
        NewtypeDiv, NewtypeDiv(&self),
        NewtypeMul, NewtypeMul(&self),
        NewtypeRem, NewtypeRem(&self),
        NewtypeSub, NewtypeSub(&self),

        NewtypeShl(), NewtypeShl(&self), NewtypeShl(usize), NewtypeShl(&self, usize),
        NewtypeShr(), NewtypeShr(&self), NewtypeShr(usize), NewtypeShr(&self, usize),

        NewtypeNeg, NewtypeNeg(&self),
        NewtypeNot, NewtypeNot(&self),

        NewtypeFrom
        )]
    pub struct Dummy1(pub i32);
}

custom_derive! {
    #[derive(Clone, Eq, PartialEq, Debug,
        NewtypeFrom,
        NewtypeDeref, NewtypeDerefMut,
        NewtypeIndex(usize), NewtypeIndexMut(usize)
        )]
    pub struct Dummy2(Vec<i32>);
}

struct Dummy3Inner;

impl_fmt!(impl Binary for Dummy3Inner: "binary");
impl_fmt!(impl Debug for Dummy3Inner: "debug");
impl_fmt!(impl Display for Dummy3Inner: "display");
impl_fmt!(impl LowerExp for Dummy3Inner: "lowerexp");
impl_fmt!(impl LowerHex for Dummy3Inner: "lowerhex");
impl_fmt!(impl Octal for Dummy3Inner: "octal");
impl_fmt!(impl Pointer for Dummy3Inner: "pointer");
impl_fmt!(impl UpperExp for Dummy3Inner: "upperexp");
impl_fmt!(impl UpperHex for Dummy3Inner: "upperhex");

custom_derive! {
    #[derive(
        NewtypeBinary,
        NewtypeDebug,
        NewtypeDisplay,
        NewtypeLowerExp,
        NewtypeLowerHex,
        NewtypeOctal,
        NewtypePointer,
        NewtypeUpperExp,
        NewtypeUpperHex
    )]
    struct Dummy3(Dummy3Inner);
}

#[test]
fn test_pub_interior_fields() {
    let _ = Dummy1(0);
    let _ = Dummy2(vec![0]);
    let _ = Dummy3(Dummy3Inner);
}

#[cfg(feature = "std-unstable")]
mod std_unstable {
    custom_derive! {
        #[derive(Eq, PartialEq, Debug, NewtypeOne, NewtypeZero)]
        struct Dummy4(pub i32);
    }

    #[test]
    fn test_pub_interior_fields_std_unstable() {
        let _ = Dummy4(0);
    }
}
