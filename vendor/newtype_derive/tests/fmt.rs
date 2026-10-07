/*
Copyright ⓒ 2015 rust-custom-derive contributors.

Licensed under the MIT license (see LICENSE or <http://opensource.org
/licenses/MIT>) or the Apache License, Version 2.0 (see LICENSE of
<http://www.apache.org/licenses/LICENSE-2.0>), at your option. All
files in the project carrying such notice may not be copied, modified,
or distributed except according to those terms.
*/
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

struct Dummy;

impl_fmt!(impl Binary for Dummy: "binary");
impl_fmt!(impl Debug for Dummy: "debug");
impl_fmt!(impl Display for Dummy: "display");
impl_fmt!(impl LowerExp for Dummy: "lowerexp");
impl_fmt!(impl LowerHex for Dummy: "lowerhex");
impl_fmt!(impl Octal for Dummy: "octal");
impl_fmt!(impl Pointer for Dummy: "pointer");
impl_fmt!(impl UpperExp for Dummy: "upperexp");
impl_fmt!(impl UpperHex for Dummy: "upperhex");

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
    struct Wrapper(Dummy);
}

#[test]
fn test_fmt() {
    let a = Wrapper(Dummy);

    assert_eq!(&*format!("{:b}", a), "binary");
    assert_eq!(&*format!("{:?}", a), "debug");
    assert_eq!(&*format!("{}", a), "display");
    assert_eq!(&*format!("{:e}", a), "lowerexp");
    assert_eq!(&*format!("{:x}", a), "lowerhex");
    assert_eq!(&*format!("{:o}", a), "octal");
    assert_eq!(&*format!("{:p}", a), "pointer");
    assert_eq!(&*format!("{:E}", a), "upperexp");
    assert_eq!(&*format!("{:X}", a), "upperhex");
}
