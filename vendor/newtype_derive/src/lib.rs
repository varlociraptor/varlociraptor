/*
Copyright ⓒ 2015 rust-custom-derive contributors.

Licensed under the MIT license (see LICENSE or <http://opensource.org
/licenses/MIT>) or the Apache License, Version 2.0 (see LICENSE of
<http://www.apache.org/licenses/LICENSE-2.0>), at your option. All
files in the project carrying such notice may not be copied, modified,
or distributed except according to those terms.
*/
/*!
This crate provides several macros for deriving implementations of various traits for "newtype" wrappers (*i.e.* tuple structs with a single element).  That is, given a tuple struct with exactly one field (*e.g.* `struct Buckets(i32)`), these macros will derive "obvious" implementations of traits such as `Add`, `Neg`, `Index`, `Deref`, `From`, etc.

All of these macros are designed to be used with the [`custom_derive`](https://crates.io/crates/custom_derive) crate, though they can be used independent of it.

# Example

Create a simple integer wrapper with some arithmetic operators:

```rust
#[macro_use] extern crate custom_derive;
#[macro_use] extern crate newtype_derive;

custom_derive! {
    #[derive(NewtypeFrom, NewtypeAdd, NewtypeMul(i32))]
    pub struct Happy(i32);
}

# fn main() {
// Let's add some happy little ints.
let a = Happy::from(6);
let b = Happy::from(7);
let c = (a + b) * 3;
let d: i32 = c.into();
assert_eq!(d, 39);
# }
```

Create a "deref-transparent" wrapper around a type:

```rust
#[macro_use] extern crate custom_derive;
#[macro_use] extern crate newtype_derive;

custom_derive! {
    #[derive(NewtypeFrom,
        NewtypeDeref, NewtypeDerefMut,
        NewtypeIndex(usize), NewtypeIndexMut(usize)
        )]
    pub struct I32Array(Vec<i32>);
}

# fn main() {
let mut arr = I32Array::from(vec![1, 2, 3]);
arr.push(4);
arr[2] = 5;
assert_eq!(&**arr, &[1, 2, 5, 4]);
assert_eq!(arr.len(), 4);
# }
```

# Overview

This crate provides macros to derive implementations of the following traits for newtype structs:

- Binary Arithmetic Operators: Add, BitAnd, BitOr, BitXor, Div, Mul, Rem, Sub, Shl, Shr, plus the corresponding *Assign traits.
- Unary Arithmetic Operators: Neg, Not.
- Other Operators: Deref, DerefMut, Index, IndexMut.
- Formatting: Binary, Debug, Display, LowerExp, LowerHex, Octal, Pointer, UpperExp, UpperHex.
- Miscellaneous: From.
- Unstable: One, Product, Sum, Zero (requires the `std-unstable` feature).

All of these macros are named `Newtype$Trait`.

None of these macros currently support generic newtype structs.

## Binary Arithmetic Operators

Each of the binary arithmetic operators accept several deriving forms.  To use `Add` on a struct `T` as an example:

- `NewtypeAdd`: `impl Add<T, Output=T> for T`
- `NewtypeAdd(&self)`: `impl<'a> Add<&'a T, Output=T> for &'a T`
- `NewtypeAdd(U)`: `impl Add<U, Output=T> for T`
- `NewtypeAdd(&self, U)`: `impl<'a> Add<U, Output=T> for &'a T`
- `NewtypeAdd(*)`: All four combinations of `T` and `&T`

The `*Assign` variants accept zero or one argument only.  For example:

- `NewtypeAddAssign`: `impl AddAssign<T> for T`
- `NewtypeAddAssign(&Self)`: `impl<'a> Add<&'a T> for &'a T`
- `NewtypeAddAssign(U)`: `impl Add<U> for T`
- `NewtypeAddAssign(*)`: Implements for `T` and `&T`.

In all cases, the implementation unwraps the newtype (where necessary), forwards to the wrapped value's implementation, then re-wraps the result in the newtype.

## Unary Arithmetic Operators

Each of the binary arithmetic operators accept several deriving forms.  To use `Neg` on a struct `T` as an example:

- `NewtypeNeg`: `impl Neg<Output=T> for T`
- `NewtypeNeg(&self)`: `impl<'a> Neg<Output=T> for &'a T`
- `NewtypeNeg(*)`: both of the above

In all cases, the implementation unwraps the newtype, forwards to the wrapped value's implementation, then re-wraps the result in the newtype.

## Other Operators

`NewtypeDeref` and `NewtypeDerefMut` only support the argument-less form, and implements the corresponding trait such that the newtype structure derefs to a pointer to the wrapped value.

`NewtypeIndex` and `NewtypeIndexMut` must be used as `NewtypeIndex(usize)`, where the argument is the type to use for indexing.  The call is forwarded to the wrapped value's implementation.

## Formatting

The deriving macros for the formatting traits in [`std::fmt`][] forward to the wrapped value's implementation.

[`std::fmt`]: http://doc.rust-lang.org/std/fmt/index.html

## Miscellaneous

`NewtypeFrom` implements `std::convert::From` twice: once for converting from the wrapped type to the newtype, and once for converting from the newtype to the wrapped type.

`NewtypeProduct` and `NewtypeSum` optionally support specifying `&Self` as an argument to generate an implementation that accepts an iterator of borrowed pointers (*e.g.* `NewtypeSum(&Self)`).

## Using Without `custom_derive!`

Although designed to be used with `custom_derive!`, all of the macros in this crate can be used without it.  The following:

```rust
# #[macro_use] extern crate custom_derive;
# #[macro_use] extern crate newtype_derive;
custom_derive! {
    #[derive(Copy, Clone, Debug, NewtypeFrom, NewtypeAdd, NewtypeAdd(f32))]
    pub struct Meters(f32);
}
# fn main() {}
```

Can also be written as:

```rust
# #[macro_use] extern crate newtype_derive;
#[derive(Copy, Clone, Debug)]
pub struct Meters(f32);

NewtypeFrom! { () pub struct Meters(f32); }
NewtypeAdd! { () pub struct Meters(f32); }
NewtypeAdd! { (f32) pub struct Meters(f32); }
# fn main() {}
```
*/
/*
# `Newtype$binop` Template

Given `/\/\/\s*(ntbop\s+([A-Za-z0-9]+),\s*([a-z_]+))\n(^#\[.+?\]$\n)*^macro_rules!.*?\{$\n(^ +.*?$\n)*^\}$/`,

```
// \1
#[macro_export]
macro_rules! Newtype\2 {
    ((*) $($tts:tt)*) => {
        Newtype\2! { () $($tts)* }
        Newtype\2! { (&self) $($tts)* }
        Newtype\2! { (&Self) $($tts)* }
        Newtype\2! { (&self, Self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::\2)::\3, kind: simple, item: $($tts)* }
    };
    ((&self) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::\2)::\3, kind: simple_ref, item: $($tts)* }
    };
    ((&self, $($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::\2)::\3, kind: ref_rhs_rewrap($($rhs)*), item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::\2)::\3, kind: rhs_rewrap($($rhs)*), item: $($tts)* }
    };
}
```

# `Newtype$binopass` Template

Given `/\/\/\s*(ntbopass\s+([A-Za-z0-9]+),\s*([a-z_]+))\n(^#\[.+?\]$\n)*^macro_rules!.*?\{$\n(^ +.*?$\n)*^\}$/`,

```
// \1
#[macro_export]
#[cfg(op_assign)]
macro_rules! Newtype\2 {
    ((*) $($tts:tt)*) => {
        Newtype\2! { () $($tts)* }
        Newtype\2! { (&self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::\2)::\3, kind: simple, item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::\2)::\3, kind: rhs($($rhs)*), item: $($tts)* }
    };
}
```

# `Newtype$unop` Template

Given `/\/\/\s*(ntuop\s+([A-Za-z0-9]+),\s*([a-z_]+))\n(^#\[.+?\]$\n)*^macro_rules!.*?\{$\n(^ +.*?$\n)*^\}$/`,

```
// \1
#[macro_export]
macro_rules! Newtype\2 {
    ((*) $($tts:tt)*) => {
        Newtype\2! { () $($tts)* }
        Newtype\2! { (&self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_un_op! { trait: (::std::ops::\2)::\3, kind: simple, item: $($tts)* }
    };
    ((&self) $($tts:tt)*) => {
        newtype_wrap_un_op! { trait: (::std::ops::\2)::\3, kind: simple_ref, item: $($tts)* }
    };
}
```
*/
#![cfg_attr(not(feature = "std"), no_std)]

mod std_unstable;

#[doc(hidden)]
#[macro_export]
macro_rules! newtype_as_item {
    ($i:item) => {$i};
}

#[doc(hidden)]
#[macro_export]
macro_rules! newtype_wrap_bin_op {
    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: simple,
        item: $(pub)* struct $name:ident(pub $t:ty);
    ) => {
        newtype_as_item! {
            impl $($tr)*<$name> for $name {
                type Output = $name;
                fn $meth(self, rhs: Self) -> $name {
                    $name((self.0).$meth(rhs.0))
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: simple,
        item: $(pub)* struct $name:ident($t:ty);
    ) => {
        newtype_as_item! {
            impl $($tr)*<$name> for $name {
                type Output = $name;
                fn $meth(self, rhs: Self) -> $name {
                    $name((self.0).$meth(rhs.0))
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: simple_ref,
        item: $(pub)* struct $name:ident(pub $t:ty);
    ) => {
        newtype_as_item! {
            impl<'a> $($tr)*<&'a $name> for &'a $name {
                type Output = $name;
                fn $meth(self, rhs: Self) -> $name {
                    $name((self.0).$meth(rhs.0))
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: simple_ref,
        item: $(pub)* struct $name:ident($t:ty);
    ) => {
        newtype_as_item! {
            impl<'a> $($tr)*<&'a $name> for &'a $name {
                type Output = $name;
                fn $meth(self, rhs: Self) -> $name {
                    $name((self.0).$meth(rhs.0))
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: rhs_rewrap(&Self),
        item: $(pub)* struct $name:ident(pub $t:ty);
    ) => {
        newtype_as_item! {
            impl<'a> $($tr)*<&'a $name> for $name {
                type Output = $name;
                fn $meth(self, rhs: &'a $name) -> $name {
                    $name((self.0).$meth(&rhs.0))
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: rhs_rewrap(&Self),
        item: $(pub)* struct $name:ident($t:ty);
    ) => {
        newtype_as_item! {
            impl<'a> $($tr)*<&'a $name> for $name {
                type Output = $name;
                fn $meth(self, rhs: &'a $name) -> $name {
                    $name((self.0).$meth(&rhs.0))
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: rhs_rewrap($rhs:ty),
        item: $(pub)* struct $name:ident(pub $t:ty);
    ) => {
        newtype_as_item! {
            impl $($tr)*<$rhs> for $name {
                type Output = $name;
                fn $meth(self, rhs: $rhs) -> $name {
                    $name((self.0).$meth(rhs))
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: rhs_rewrap($rhs:ty),
        item: $(pub)* struct $name:ident($t:ty);
    ) => {
        newtype_as_item! {
            impl $($tr)*<$rhs> for $name {
                type Output = $name;
                fn $meth(self, rhs: $rhs) -> $name {
                    $name((self.0).$meth(rhs))
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: ref_rhs_rewrap(Self),
        item: $(pub)* struct $name:ident(pub $t:ty);
    ) => {
        newtype_as_item! {
            impl<'a> $($tr)*<$name> for &'a $name {
                type Output = $name;
                fn $meth(self, rhs: $name) -> $name {
                    $name((self.0).$meth(rhs.0))
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: ref_rhs_rewrap(Self),
        item: $(pub)* struct $name:ident($t:ty);
    ) => {
        newtype_as_item! {
            impl<'a> $($tr)*<$name> for &'a $name {
                type Output = $name;
                fn $meth(self, rhs: $name) -> $name {
                    $name((self.0).$meth(rhs.0))
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: ref_rhs_rewrap($rhs:ty),
        item: $(pub)* struct $name:ident(pub $t:ty);
    ) => {
        newtype_as_item! {
            impl<'a> $($tr)*<$rhs> for &'a $name {
                type Output = $name;
                fn $meth(self, rhs: $rhs) -> $name {
                    $name((self.0).$meth(rhs))
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: ref_rhs_rewrap($rhs:ty),
        item: $(pub)* struct $name:ident($t:ty);
    ) => {
        newtype_as_item! {
            impl<'a> $($tr)*<$rhs> for &'a $name {
                type Output = $name;
                fn $meth(self, rhs: $rhs) -> $name {
                    $name((self.0).$meth(rhs))
                }
            }
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! newtype_wrap_bin_op_assign {
    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: simple,
        item: $(pub)* struct $name:ident(pub $t:ty);
    ) => {
        newtype_as_item! {
            impl $($tr)*<$name> for $name {
                fn $meth(&mut self, rhs: Self) {
                    (self.0).$meth(rhs.0)
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: simple,
        item: $(pub)* struct $name:ident($t:ty);
    ) => {
        newtype_as_item! {
            impl $($tr)*<$name> for $name {
                fn $meth(&mut self, rhs: Self) {
                    (self.0).$meth(rhs.0)
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: rhs(&Self),
        item: $(pub)* struct $name:ident(pub $t:ty);
    ) => {
        newtype_as_item! {
            impl<'a> $($tr)*<&'a $name> for $name {
                fn $meth(&mut self, rhs: &'a $name) {
                    (self.0).$meth(rhs.0)
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: rhs(&Self),
        item: $(pub)* struct $name:ident($t:ty);
    ) => {
        newtype_as_item! {
            impl<'a> $($tr)*<&'a $name> for $name {
                fn $meth(&mut self, rhs: &'a $name) {
                    (self.0).$meth(rhs.0)
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: rhs($rhs:ty),
        item: $(pub)* struct $name:ident(pub $t:ty);
    ) => {
        newtype_as_item! {
            impl $($tr)*<$rhs> for $name {
                fn $meth(&mut self, rhs: $rhs) {
                    (self.0).$meth(rhs)
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: rhs($rhs:ty),
        item: $(pub)* struct $name:ident($t:ty);
    ) => {
        newtype_as_item! {
            impl $($tr)*<$rhs> for $name {
                fn $meth(&mut self, rhs: $rhs) {
                    (self.0).$meth(rhs)
                }
            }
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! newtype_wrap_un_op {
    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: simple,
        item: $(pub)* struct $name:ident(pub $t:ty);
    ) => {
        newtype_as_item! {
            impl $($tr)* for $name {
                type Output = $name;
                fn $meth(self) -> $name {
                    $name((self.0).$meth())
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: simple,
        item: $(pub)* struct $name:ident($t:ty);
    ) => {
        newtype_as_item! {
            impl $($tr)* for $name {
                type Output = $name;
                fn $meth(self) -> $name {
                    $name((self.0).$meth())
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: simple_ref,
        item: $(pub)* struct $name:ident(pub $t:ty);
    ) => {
        newtype_as_item! {
            impl<'a> $($tr)* for &'a $name {
                type Output = $name;
                fn $meth(self) -> $name {
                    $name((self.0).$meth())
                }
            }
        }
    };

    (
        trait: ($($tr:tt)*)::$meth:ident,
        kind: simple_ref,
        item: $(pub)* struct $name:ident($t:ty);
    ) => {
        newtype_as_item! {
            impl<'a> $($tr)* for &'a $name {
                type Output = $name;
                fn $meth(self) -> $name {
                    $name((self.0).$meth())
                }
            }
        }
    };
}

// ntbop Add,      add
#[macro_export]
macro_rules! NewtypeAdd {
    ((*) $($tts:tt)*) => {
        NewtypeAdd! { () $($tts)* }
        NewtypeAdd! { (&self) $($tts)* }
        NewtypeAdd! { (&Self) $($tts)* }
        NewtypeAdd! { (&self, Self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Add)::add, kind: simple, item: $($tts)* }
    };
    ((&self) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Add)::add, kind: simple_ref, item: $($tts)* }
    };
    ((&self, $($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Add)::add, kind: ref_rhs_rewrap($($rhs)*), item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Add)::add, kind: rhs_rewrap($($rhs)*), item: $($tts)* }
    };
}

// ntbopass AddAssign, add_assign
#[macro_export]
#[cfg(op_assign)]
macro_rules! NewtypeAddAssign {
    ((*) $($tts:tt)*) => {
        NewtypeAddAssign! { () $($tts)* }
        NewtypeAddAssign! { (&self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::AddAssign)::add_assign, kind: simple, item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::AddAssign)::add_assign, kind: rhs($($rhs)*), item: $($tts)* }
    };
}

// ntbop BitAnd,   bitand
#[macro_export]
macro_rules! NewtypeBitAnd {
    ((*) $($tts:tt)*) => {
        NewtypeBitAnd! { () $($tts)* }
        NewtypeBitAnd! { (&self) $($tts)* }
        NewtypeBitAnd! { (&Self) $($tts)* }
        NewtypeBitAnd! { (&self, Self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::BitAnd)::bitand, kind: simple, item: $($tts)* }
    };
    ((&self) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::BitAnd)::bitand, kind: simple_ref, item: $($tts)* }
    };
    ((&self, $($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::BitAnd)::bitand, kind: ref_rhs_rewrap($($rhs)*), item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::BitAnd)::bitand, kind: rhs_rewrap($($rhs)*), item: $($tts)* }
    };
}

// ntbopass BitAndAssign, bitand_assign
#[macro_export]
#[cfg(op_assign)]
macro_rules! NewtypeBitAndAssign {
    ((*) $($tts:tt)*) => {
        NewtypeBitAndAssign! { () $($tts)* }
        NewtypeBitAndAssign! { (&self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::BitAndAssign)::bitand_assign, kind: simple, item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::BitAndAssign)::bitand_assign, kind: rhs($($rhs)*), item: $($tts)* }
    };
}

// ntbop BitOr,    bitor
#[macro_export]
macro_rules! NewtypeBitOr {
    ((*) $($tts:tt)*) => {
        NewtypeBitOr! { () $($tts)* }
        NewtypeBitOr! { (&self) $($tts)* }
        NewtypeBitOr! { (&Self) $($tts)* }
        NewtypeBitOr! { (&self, Self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::BitOr)::bitor, kind: simple, item: $($tts)* }
    };
    ((&self) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::BitOr)::bitor, kind: simple_ref, item: $($tts)* }
    };
    ((&self, $($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::BitOr)::bitor, kind: ref_rhs_rewrap($($rhs)*), item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::BitOr)::bitor, kind: rhs_rewrap($($rhs)*), item: $($tts)* }
    };
}

// ntbopass BitOrAssign, bitor_assign
#[macro_export]
#[cfg(op_assign)]
macro_rules! NewtypeBitOrAssign {
    ((*) $($tts:tt)*) => {
        NewtypeBitOrAssign! { () $($tts)* }
        NewtypeBitOrAssign! { (&self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::BitOrAssign)::bitor_assign, kind: simple, item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::BitOrAssign)::bitor_assign, kind: rhs($($rhs)*), item: $($tts)* }
    };
}

// ntbop BitXor,   bitxor
#[macro_export]
macro_rules! NewtypeBitXor {
    ((*) $($tts:tt)*) => {
        NewtypeBitXor! { () $($tts)* }
        NewtypeBitXor! { (&self) $($tts)* }
        NewtypeBitXor! { (&Self) $($tts)* }
        NewtypeBitXor! { (&self, Self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::BitXor)::bitxor, kind: simple, item: $($tts)* }
    };
    ((&self) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::BitXor)::bitxor, kind: simple_ref, item: $($tts)* }
    };
    ((&self, $($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::BitXor)::bitxor, kind: ref_rhs_rewrap($($rhs)*), item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::BitXor)::bitxor, kind: rhs_rewrap($($rhs)*), item: $($tts)* }
    };
}

// ntbopass BitXorAssign, bitxor_assign
#[macro_export]
#[cfg(op_assign)]
macro_rules! NewtypeBitXorAssign {
    ((*) $($tts:tt)*) => {
        NewtypeBitXorAssign! { () $($tts)* }
        NewtypeBitXorAssign! { (&self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::BitXorAssign)::bitxor_assign, kind: simple, item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::BitXorAssign)::bitxor_assign, kind: rhs($($rhs)*), item: $($tts)* }
    };
}

// ntbop Div,      div
#[macro_export]
macro_rules! NewtypeDiv {
    ((*) $($tts:tt)*) => {
        NewtypeDiv! { () $($tts)* }
        NewtypeDiv! { (&self) $($tts)* }
        NewtypeDiv! { (&Self) $($tts)* }
        NewtypeDiv! { (&self, Self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Div)::div, kind: simple, item: $($tts)* }
    };
    ((&self) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Div)::div, kind: simple_ref, item: $($tts)* }
    };
    ((&self, $($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Div)::div, kind: ref_rhs_rewrap($($rhs)*), item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Div)::div, kind: rhs_rewrap($($rhs)*), item: $($tts)* }
    };
}

// ntbopass DivAssign, div_assign
#[macro_export]
#[cfg(op_assign)]
macro_rules! NewtypeDivAssign {
    ((*) $($tts:tt)*) => {
        NewtypeDivAssign! { () $($tts)* }
        NewtypeDivAssign! { (&self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::DivAssign)::div_assign, kind: simple, item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::DivAssign)::div_assign, kind: rhs($($rhs)*), item: $($tts)* }
    };
}

// ntbop Mul,      mul
#[macro_export]
macro_rules! NewtypeMul {
    ((*) $($tts:tt)*) => {
        NewtypeMul! { () $($tts)* }
        NewtypeMul! { (&self) $($tts)* }
        NewtypeMul! { (&Self) $($tts)* }
        NewtypeMul! { (&self, Self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Mul)::mul, kind: simple, item: $($tts)* }
    };
    ((&self) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Mul)::mul, kind: simple_ref, item: $($tts)* }
    };
    ((&self, $($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Mul)::mul, kind: ref_rhs_rewrap($($rhs)*), item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Mul)::mul, kind: rhs_rewrap($($rhs)*), item: $($tts)* }
    };
}

// ntbopass MulAssign, mul_assign
#[macro_export]
#[cfg(op_assign)]
macro_rules! NewtypeMulAssign {
    ((*) $($tts:tt)*) => {
        NewtypeMulAssign! { () $($tts)* }
        NewtypeMulAssign! { (&self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::MulAssign)::mul_assign, kind: simple, item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::MulAssign)::mul_assign, kind: rhs($($rhs)*), item: $($tts)* }
    };
}

// ntbop Rem,      rem
#[macro_export]
macro_rules! NewtypeRem {
    ((*) $($tts:tt)*) => {
        NewtypeRem! { () $($tts)* }
        NewtypeRem! { (&self) $($tts)* }
        NewtypeRem! { (&Self) $($tts)* }
        NewtypeRem! { (&self, Self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Rem)::rem, kind: simple, item: $($tts)* }
    };
    ((&self) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Rem)::rem, kind: simple_ref, item: $($tts)* }
    };
    ((&self, $($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Rem)::rem, kind: ref_rhs_rewrap($($rhs)*), item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Rem)::rem, kind: rhs_rewrap($($rhs)*), item: $($tts)* }
    };
}

// ntbopass RemAssign, rem_assign
#[macro_export]
#[cfg(op_assign)]
macro_rules! NewtypeRemAssign {
    ((*) $($tts:tt)*) => {
        NewtypeRemAssign! { () $($tts)* }
        NewtypeRemAssign! { (&self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::RemAssign)::rem_assign, kind: simple, item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::RemAssign)::rem_assign, kind: rhs($($rhs)*), item: $($tts)* }
    };
}

// ntbop Sub,      sub
#[macro_export]
macro_rules! NewtypeSub {
    ((*) $($tts:tt)*) => {
        NewtypeSub! { () $($tts)* }
        NewtypeSub! { (&self) $($tts)* }
        NewtypeSub! { (&Self) $($tts)* }
        NewtypeSub! { (&self, Self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Sub)::sub, kind: simple, item: $($tts)* }
    };
    ((&self) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Sub)::sub, kind: simple_ref, item: $($tts)* }
    };
    ((&self, $($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Sub)::sub, kind: ref_rhs_rewrap($($rhs)*), item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Sub)::sub, kind: rhs_rewrap($($rhs)*), item: $($tts)* }
    };
}

// ntbopass SubAssign, sub_assign
#[macro_export]
#[cfg(op_assign)]
macro_rules! NewtypeSubAssign {
    ((*) $($tts:tt)*) => {
        NewtypeSubAssign! { () $($tts)* }
        NewtypeSubAssign! { (&self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::SubAssign)::sub_assign, kind: simple, item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::SubAssign)::sub_assign, kind: rhs($($rhs)*), item: $($tts)* }
    };
}

// ntbop Shl,      shl
#[macro_export]
macro_rules! NewtypeShl {
    ((*) $($tts:tt)*) => {
        NewtypeShl! { () $($tts)* }
        NewtypeShl! { (&self) $($tts)* }
        NewtypeShl! { (&Self) $($tts)* }
        NewtypeShl! { (&self, Self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Shl)::shl, kind: simple, item: $($tts)* }
    };
    ((&self) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Shl)::shl, kind: simple_ref, item: $($tts)* }
    };
    ((&self, $($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Shl)::shl, kind: ref_rhs_rewrap($($rhs)*), item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Shl)::shl, kind: rhs_rewrap($($rhs)*), item: $($tts)* }
    };
}

// ntbopass ShlAssign, shl_assign
#[macro_export]
#[cfg(op_assign)]
macro_rules! NewtypeShlAssign {
    ((*) $($tts:tt)*) => {
        NewtypeShlAssign! { () $($tts)* }
        NewtypeShlAssign! { (&self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::ShlAssign)::shl_assign, kind: simple, item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::ShlAssign)::shl_assign, kind: rhs($($rhs)*), item: $($tts)* }
    };
}

// ntbop Shr,      shr
#[macro_export]
macro_rules! NewtypeShr {
    ((*) $($tts:tt)*) => {
        NewtypeShr! { () $($tts)* }
        NewtypeShr! { (&self) $($tts)* }
        NewtypeShr! { (&Self) $($tts)* }
        NewtypeShr! { (&self, Self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Shr)::shr, kind: simple, item: $($tts)* }
    };
    ((&self) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Shr)::shr, kind: simple_ref, item: $($tts)* }
    };
    ((&self, $($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Shr)::shr, kind: ref_rhs_rewrap($($rhs)*), item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op! { trait: (::std::ops::Shr)::shr, kind: rhs_rewrap($($rhs)*), item: $($tts)* }
    };
}

// ntbopass ShrAssign, shr_assign
#[macro_export]
#[cfg(op_assign)]
macro_rules! NewtypeShrAssign {
    ((*) $($tts:tt)*) => {
        NewtypeShrAssign! { () $($tts)* }
        NewtypeShrAssign! { (&self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::ShrAssign)::shr_assign, kind: simple, item: $($tts)* }
    };
    (($($rhs:tt)*) $($tts:tt)*) => {
        newtype_wrap_bin_op_assign! { trait: (::std::ops::ShrAssign)::shr_assign, kind: rhs($($rhs)*), item: $($tts)* }
    };
}

// ntuop Neg,      neg
#[macro_export]
macro_rules! NewtypeNeg {
    ((*) $($tts:tt)*) => {
        NewtypeNeg! { () $($tts)* }
        NewtypeNeg! { (&self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_un_op! { trait: (::std::ops::Neg)::neg, kind: simple, item: $($tts)* }
    };
    ((&self) $($tts:tt)*) => {
        newtype_wrap_un_op! { trait: (::std::ops::Neg)::neg, kind: simple_ref, item: $($tts)* }
    };
}

// ntuop Not,      not
#[macro_export]
macro_rules! NewtypeNot {
    ((*) $($tts:tt)*) => {
        NewtypeNot! { () $($tts)* }
        NewtypeNot! { (&self) $($tts)* }
    };
    (() $($tts:tt)*) => {
        newtype_wrap_un_op! { trait: (::std::ops::Not)::not, kind: simple, item: $($tts)* }
    };
    ((&self) $($tts:tt)*) => {
        newtype_wrap_un_op! { trait: (::std::ops::Not)::not, kind: simple_ref, item: $($tts)* }
    };
}

#[macro_export]
macro_rules! NewtypeDeref {
    (() $(pub)* struct $name:ident(pub $t0:ty);) => {
        impl ::std::ops::Deref for $name {
            type Target = $t0;
            fn deref(&self) -> &Self::Target {
                &self.0
            }
        }
    };

    (() $(pub)* struct $name:ident($t0:ty);) => {
        impl ::std::ops::Deref for $name {
            type Target = $t0;
            fn deref(&self) -> &Self::Target {
                &self.0
            }
        }
    };
}

#[macro_export]
macro_rules! NewtypeDerefMut {
    (() $(pub)* struct $name:ident(pub $t0:ty);) => {
        impl ::std::ops::DerefMut for $name {
            fn deref_mut(&mut self) -> &mut Self::Target {
                &mut self.0
            }
        }
    };

    (() $(pub)* struct $name:ident($t0:ty);) => {
        impl ::std::ops::DerefMut for $name {
            fn deref_mut(&mut self) -> &mut Self::Target {
                &mut self.0
            }
        }
    };
}

#[macro_export]
macro_rules! NewtypeIndex {
    (($index_ty:ty) $(pub)* struct $name:ident(pub $t0:ty);) => {
        impl ::std::ops::Index<$index_ty> for $name {
            type Output = <$t0 as ::std::ops::Index<$index_ty>>::Output;
            fn index(&self, index: $index_ty) -> &Self::Output {
                (&self.0).index(index)
            }
        }
    };

    (($index_ty:ty) $(pub)* struct $name:ident($t0:ty);) => {
        impl ::std::ops::Index<$index_ty> for $name {
            type Output = <$t0 as ::std::ops::Index<$index_ty>>::Output;
            fn index(&self, index: $index_ty) -> &Self::Output {
                (&self.0).index(index)
            }
        }
    };
}

#[macro_export]
macro_rules! NewtypeIndexMut {
    (($index_ty:ty) $(pub)* struct $name:ident(pub $t0:ty);) => {
        impl ::std::ops::IndexMut<$index_ty> for $name {
            fn index_mut(&mut self, index: $index_ty) -> &mut Self::Output {
                (&mut self.0).index_mut(index)
            }
        }
    };

    (($index_ty:ty) $(pub)* struct $name:ident($t0:ty);) => {
        impl ::std::ops::IndexMut<$index_ty> for $name {
            fn index_mut(&mut self, index: $index_ty) -> &mut Self::Output {
                (&mut self.0).index_mut(index)
            }
        }
    };
}

#[macro_export]
macro_rules! NewtypeFrom {
    (() $(pub)* struct $name:ident(pub $t0:ty);) => {
        impl ::std::convert::From<$t0> for $name {
            fn from(v: $t0) -> Self {
                $name(v)
            }
        }
        impl ::std::convert::From<$name> for $t0 {
            fn from(v: $name) -> Self {
                v.0
            }
        }
    };

    (() $(pub)* struct $name:ident($t0:ty);) => {
        impl ::std::convert::From<$t0> for $name {
            fn from(v: $t0) -> Self {
                $name(v)
            }
        }
        impl ::std::convert::From<$name> for $t0 {
            fn from(v: $name) -> Self {
                v.0
            }
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! newtype_fmt {
    ($fmt_trait:ident, $name:ident) => {
        impl ::std::fmt::$fmt_trait for $name {
            fn fmt(&self, fmt: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
                ::std::fmt::$fmt_trait::fmt(&self.0, fmt)
            }
        }
    };
}

#[macro_export]
macro_rules! NewtypeBinary {
    (() $(pub)* struct $name:ident $_field:tt;) => {
        newtype_fmt! { Binary, $name }
    };
}

#[macro_export]
macro_rules! NewtypeDebug {
    (() $(pub)* struct $name:ident $_field:tt;) => {
        newtype_fmt! { Debug, $name }
    };
}

#[macro_export]
macro_rules! NewtypeDisplay {
    (() $(pub)* struct $name:ident $_field:tt;) => {
        newtype_fmt! { Display, $name }
    };
}

#[macro_export]
macro_rules! NewtypeLowerExp {
    (() $(pub)* struct $name:ident $_field:tt;) => {
        newtype_fmt! { LowerExp, $name }
    };
}

#[macro_export]
macro_rules! NewtypeLowerHex {
    (() $(pub)* struct $name:ident $_field:tt;) => {
        newtype_fmt! { LowerHex, $name }
    };
}

#[macro_export]
macro_rules! NewtypeOctal {
    (() $(pub)* struct $name:ident $_field:tt;) => {
        newtype_fmt! { Octal, $name }
    };
}

#[macro_export]
macro_rules! NewtypePointer {
    (() $(pub)* struct $name:ident $_field:tt;) => {
        newtype_fmt! { Pointer, $name }
    };
}

#[macro_export]
macro_rules! NewtypeUpperExp {
    (() $(pub)* struct $name:ident $_field:tt;) => {
        newtype_fmt! { UpperExp, $name }
    };
}

#[macro_export]
macro_rules! NewtypeUpperHex {
    (() $(pub)* struct $name:ident $_field:tt;) => {
        newtype_fmt! { UpperHex, $name }
    };
}
