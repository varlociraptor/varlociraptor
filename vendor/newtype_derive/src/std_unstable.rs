/*
Copyright ⓒ 2015 rust-custom-derive contributors.

Licensed under the MIT license (see LICENSE or <http://opensource.org
/licenses/MIT>) or the Apache License, Version 2.0 (see LICENSE of
<http://www.apache.org/licenses/LICENSE-2.0>), at your option. All
files in the project carrying such notice may not be copied, modified,
or distributed except according to those terms.
*/
#![cfg(feature = "std-unstable")]

#[macro_export]
macro_rules! NewtypeOne {
    (() $(pub)* struct $name:ident(pub $_t0:ty);) => {
        impl ::std::num::One for $name {
            fn one() -> Self {
                $name(::std::num::One::one())
            }
        }
    };

    (() $(pub)* struct $name:ident($_t0:ty);) => {
        impl ::std::num::One for $name {
            fn one() -> Self {
                $name(::std::num::One::one())
            }
        }
    };
}

#[macro_export]
macro_rules! NewtypeProduct {
    ($arg:tt $(pub)* struct $name:ident(pub $t0:ty);) => {
        NewtypeProduct! { $arg struct $name($t0); }
    };

    (() $(pub)* struct $name:ident($t0:ty);) => {
        impl ::std::iter::Product<$name> for $name {
            fn product<I>(iter: I) -> Self
            where I: Iterator<Item=$name> {
                $name(iter.map(|e| e.0).product::<$t0>())
            }
        }
    };

    ((&Self) $(pub)* struct $name:ident($t0:ty);) => {
        impl<'a> ::std::iter::Product<&'a $name> for $name {
            fn product<I>(iter: I) -> Self
            where I: Iterator<Item=&'a $name> {
                $name(iter.map(|e| &e.0).product::<$t0>())
            }
        }
    };
}

#[macro_export]
macro_rules! NewtypeSum {
    ($arg:tt $(pub)* struct $name:ident(pub $t0:ty);) => {
        NewtypeSum! { $arg struct $name($t0); }
    };

    (() $(pub)* struct $name:ident($t0:ty);) => {
        impl ::std::iter::Sum<$name> for $name {
            fn sum<I>(iter: I) -> Self
            where I: Iterator<Item=$name> {
                $name(iter.map(|e| e.0).sum::<$t0>())
            }
        }
    };

    ((&Self) $(pub)* struct $name:ident($t0:ty);) => {
        impl<'a> ::std::iter::Sum<&'a $name> for $name {
            fn sum<I>(iter: I) -> Self
            where I: Iterator<Item=&'a $name> {
                $name(iter.map(|e| &e.0).sum::<$t0>())
            }
        }
    };
}

#[macro_export]
macro_rules! NewtypeZero {
    (() $(pub)* struct $name:ident(pub $_t0:ty);) => {
        impl ::std::num::Zero for $name {
            fn zero() -> Self {
                $name(::std::num::Zero::zero())
            }
        }
    };

    (() $(pub)* struct $name:ident($_t0:ty);) => {
        impl ::std::num::Zero for $name {
            fn zero() -> Self {
                $name(::std::num::Zero::zero())
            }
        }
    };
}
