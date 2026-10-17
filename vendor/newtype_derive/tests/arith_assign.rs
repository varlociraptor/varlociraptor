/*
Copyright ⓒ 2015 rust-custom-derive contributors.

Licensed under the MIT license (see LICENSE or <http://opensource.org
/licenses/MIT>) or the Apache License, Version 2.0 (see LICENSE of
<http://www.apache.org/licenses/LICENSE-2.0>), at your option. All
files in the project carrying such notice may not be copied, modified,
or distributed except according to those terms.
*/
#![cfg(op_assign)]
#![recursion_limit = "128"]
#[macro_use] extern crate custom_derive;
#[macro_use] extern crate newtype_derive;

custom_derive! {
    #[derive(Copy, Clone, Eq, PartialEq, Debug,
        NewtypeAddAssign, NewtypeAddAssign(&Self), NewtypeAddAssign(i32),
        NewtypeBitAndAssign, NewtypeBitAndAssign(&Self),
        NewtypeBitOrAssign, NewtypeBitOrAssign(&Self),
        NewtypeBitXorAssign, NewtypeBitXorAssign(&Self),
        NewtypeDivAssign, NewtypeDivAssign(&Self),
        NewtypeMulAssign, NewtypeMulAssign(&Self),
        NewtypeRemAssign, NewtypeRemAssign(&Self),
        NewtypeSubAssign, NewtypeSubAssign(&Self),
        NewtypeShlAssign, NewtypeShlAssign(&Self), NewtypeShlAssign(i32),
        NewtypeShrAssign, NewtypeShrAssign(&Self), NewtypeShrAssign(i32),
        NewtypeFrom
        )]
    pub struct Dummy(i32);
}

macro_rules! oa {
    (@as_stmt $s:stmt) => { $s };

    ($var:ident $op:tt $rhs:expr) => {
        {
            let mut $var = $var;
            oa!(@as_stmt $var $op $rhs);
            $var
        }
    };
}

#[test]
fn test_arith_assign() {
    let a = Dummy::from(4);
    let b = Dummy::from(7);

    assert_eq!(oa!(a += b), Dummy::from(4 + 7));
    assert_eq!(oa!(a += &b), Dummy::from(4 + 7));

    assert_eq!(oa!(a += b), Dummy::from(4 + 7));
    assert_eq!(oa!(a += &b), Dummy::from(4 + 7));
    assert_eq!(oa!(a += 7), Dummy::from(4 + 7));
    assert_eq!(oa!(a += 7), Dummy::from(4 + 7));
    assert_eq!(oa!(a &= b), Dummy::from(4 & 7));
    assert_eq!(oa!(a &= &b), Dummy::from(4 & 7));
    assert_eq!(oa!(a |= b), Dummy::from(4 | 7));
    assert_eq!(oa!(a |= &b), Dummy::from(4 | 7));
    assert_eq!(oa!(a ^= b), Dummy::from(4 ^ 7));
    assert_eq!(oa!(a ^= &b), Dummy::from(4 ^ 7));
    assert_eq!(oa!(a /= b), Dummy::from(4 / 7));
    assert_eq!(oa!(a /= &b), Dummy::from(4 / 7));
    assert_eq!(oa!(a *= b), Dummy::from(4 * 7));
    assert_eq!(oa!(a *= &b), Dummy::from(4 * 7));
    assert_eq!(oa!(a %= b), Dummy::from(4 % 7));
    assert_eq!(oa!(a %= &b), Dummy::from(4 % 7));
    assert_eq!(oa!(a -= b), Dummy::from(4 - 7));
    assert_eq!(oa!(a -= &b), Dummy::from(4 - 7));

    assert_eq!(oa!(a <<= b), Dummy::from(4 << 7));
    assert_eq!(oa!(a <<= &b), Dummy::from(4 << 7));
    assert_eq!(oa!(a <<= 7), Dummy::from(4 << 7));

    assert_eq!(oa!(a >>= b), Dummy::from(4 >> 7));
    assert_eq!(oa!(a >>= &b), Dummy::from(4 >> 7));
    assert_eq!(oa!(a >>= 7), Dummy::from(4 >> 7));
}
