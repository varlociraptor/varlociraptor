/*
Copyright ⓒ 2015 rust-custom-derive contributors.

Licensed under the MIT license (see LICENSE or <http://opensource.org
/licenses/MIT>) or the Apache License, Version 2.0 (see LICENSE of
<http://www.apache.org/licenses/LICENSE-2.0>), at your option. All
files in the project carrying such notice may not be copied, modified,
or distributed except according to those terms.
*/
#![cfg(feature="std-unstable")]
#![allow(deprecated)]
#![cfg_attr(feature="std-unstable", feature(zero_one))]
#![cfg_attr(feature="std-unstable", feature(iter_arith_traits))]

#[macro_use] extern crate custom_derive;
#[macro_use] extern crate newtype_derive;

use std::num::{One, Zero};

custom_derive! {
    #[derive(
        Clone, Eq, PartialEq, Debug,
        NewtypeOne, NewtypeZero,
        NewtypeSum, NewtypeSum(&Self),
        NewtypeProduct, NewtypeProduct(&Self),
    )]
    struct Dummy(i32);
}

#[test]
fn test_one_zero() {
    assert_eq!(Dummy::zero(), Dummy(0));
    assert_eq!(Dummy::one(), Dummy(1));
}

#[test]
fn test_sum_product() {
    let dummies = &[Dummy(2), Dummy(3)];
    assert_eq!(dummies.into_iter().sum::<Dummy>(), Dummy(5));
    assert_eq!(dummies.into_iter().cloned().sum::<Dummy>(), Dummy(5));
    assert_eq!(dummies.into_iter().product::<Dummy>(), Dummy(6));
    assert_eq!(dummies.into_iter().cloned().product::<Dummy>(), Dummy(6));
}
