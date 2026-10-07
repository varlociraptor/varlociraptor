#![no_std]
//! Macros for easy operator overloading.
//!
//! The primary macro to learn is `impl_op!(<op> <closure>);`
//! where `<op>` is an operator and `<closure>` is a closure with the same signature as the trait function associated with `<op>`.
//! The macro you'll actually want to use most of the time, however, is [`impl_op_ex!`](macro.impl_op_ex.html). It works the same way as `impl_op!` but with some extra magic behind the scenes.
//!
//! To use, import the macros with `use auto_ops::*;`. Remember that you can only overload operators between one or more types defined in the current crate (respecting Rust orphan rules).
//! # Examples
//! All of the following examples are run with the following struct defined:
//!
//! ```
//! #[derive(Clone, Debug, PartialEq)]
//! struct DonkeyKong {
//!     pub bananas: i32,
//! }
//! impl DonkeyKong {
//!     pub fn new(bananas: i32) -> DonkeyKong {
//!         DonkeyKong { bananas: bananas }
//!     }
//! }
//! ```
//! ## Binary operators
//!
//! ```
//! // impl_op!(op |a: LHS, b: RHS| -> OUT {...});
//! // impl_op!(op |a: LHS, b: &RHS| -> OUT {...});
//! // impl_op!(op |a: &LHS, b: RHS| -> OUT {...});
//! // impl_op!(op |a: &LHS, b: &RHS| -> OUT {...});
//! // where
//! // OP  : +, -, *, /, %, &, |, ^, <<, >>
//! // a, b: variable names
//!
//! use auto_ops::impl_op;
//! # #[derive(Clone, Debug, PartialEq)]
//! # struct DonkeyKong {
//! #     pub bananas: i32,
//! # }
//! # impl DonkeyKong {
//! #     pub fn new(bananas: i32) -> DonkeyKong {
//! #         DonkeyKong { bananas: bananas }
//! #     }
//! # }
//!
//! impl_op!(- |a: DonkeyKong, b: i32| -> DonkeyKong { DonkeyKong::new(a.bananas - b) });
//! impl_op!(+ |a: &DonkeyKong, b: &DonkeyKong| -> i32 { a.bananas + b.bananas });
//!
//! fn main() {
//!     let dk = DonkeyKong::new(3) - 1;
//!     assert_eq!(DonkeyKong::new(2), dk);
//!     let total_bananas = &dk + &DonkeyKong::new(4);
//!     assert_eq!(6, total_bananas);
//! }
//! ```
//! ## Assignment operators
//! ```
//! // impl_op!(OP |a: &mut LHS, b: RHS| {...});
//! // impl_op!(op |a: &mut LHS, b: &RHS| {...})
//! // where
//! // op  : +=, -=, *=, /=, %=, &=, |=, ^=, <<=, >>=
//! // a, b: variable names
//!
//! use auto_ops::impl_op;
//! # #[derive(Clone, Debug, PartialEq)]
//! # struct DonkeyKong {
//! #     pub bananas: i32,
//! # }
//! # impl DonkeyKong {
//! #     pub fn new(bananas: i32) -> DonkeyKong {
//! #         DonkeyKong { bananas: bananas }
//! #     }
//! # }
//!
//! impl_op!(+= |a: &mut DonkeyKong, b: DonkeyKong| { a.bananas += b.bananas });
//! impl_op!(+= |a: &mut DonkeyKong, b: &DonkeyKong| { a.bananas += b.bananas });
//!
//! fn main() {
//!     let mut dk = DonkeyKong::new(3);
//!     dk += DonkeyKong::new(1);
//!     assert_eq!(DonkeyKong::new(4), dk);
//!     dk += &DonkeyKong::new(1);
//!     assert_eq!(DonkeyKong::new(5), dk);
//! }
//! ```
//! ## Unary operators
//! ```
//! // impl_op!(OP |a: LHS| -> OUT {...});
//! // impl_op!(op |a: &LHS| -> OUT {...})
//! // where
//! // OP: !, -
//! // a : variable name
//!
//! use auto_ops::impl_op;
//! # #[derive(Clone, Debug, PartialEq)]
//! # struct DonkeyKong {
//! #     pub bananas: i32,
//! # }
//! # impl DonkeyKong {
//! #     pub fn new(bananas: i32) -> DonkeyKong {
//! #         DonkeyKong { bananas: bananas }
//! #     }
//! # }
//!
//! impl_op!(- |a: DonkeyKong| -> DonkeyKong { DonkeyKong::new(-a.bananas) });
//! impl_op!(- |a: &DonkeyKong| -> DonkeyKong { DonkeyKong::new(-a.bananas) });
//!
//! fn main() {
//!     let dk = -DonkeyKong::new(3);
//!     assert_eq!(DonkeyKong::new(-3), dk);
//!     assert_eq!(DonkeyKong::new(3), -&dk);
//! }
//! ```
//! # Limitations
//! * The output type of any operation must be an owned type (i.e. `impl_op!(+ |a: DonkeyKong b: i32| -> &DonkeyKong {...})` is invalid).
//! * Types that have an unqualified lifetime or associated type are invalid
//!
//! ```compile_fail
//! // impl_op!(+ |a: SomeType<'a>, b: SomeType<'a>| -> SomeType<'a> {...}) // INVALID
//! // impl_op!(+ |a: SomeType<T>, b: SomeType<T>| -> SomeType<T> {...})    // INVALID
//! impl_op!(+ |a: SomeType<i32>, b: SomeType<i32>| -> SomeType<i32> {...}) // VALID
//! ```
mod assignment;
mod binary;
mod unary;

/// Overloads an operator using the given closure as its body.
///
/// See the [module level documentation](index.html) for more information.
#[macro_export]
macro_rules! impl_op {
    ($op:tt |$lhs_i:ident : &mut $lhs:ty, $rhs_i:ident : &$rhs:ty| $body:block) => {
        $crate::_parse_assignment_op!($op, $lhs, &$rhs, lhs, rhs, {
            |$lhs_i: &mut $lhs, $rhs_i: &$rhs| -> () { $body }(lhs, rhs);
        });
    };
    ($op:tt |$lhs_i:ident : &mut $lhs:ty, $rhs_i:ident : $rhs:ty| $body:block) => {
        $crate::_parse_assignment_op!($op, $lhs, $rhs, lhs, rhs, {
            |$lhs_i: &mut $lhs, $rhs_i: $rhs| -> () { $body }(lhs, rhs);
        });
    };
    ($op:tt |$lhs_i:ident : &$lhs:ty| -> $out:ty $body:block) => {
        $crate::_parse_unary_op!($op, &$lhs, $out, lhs, {
            |$lhs_i: &$lhs| -> $out { $body }(lhs)
        });
    };
    ($op:tt |$lhs_i:ident : &$lhs:ty, $rhs_i:ident : &$rhs:ty| -> $out:ty $body:block) => {
        $crate::_parse_binary_op!($op, &$lhs, &$rhs, $out, lhs, rhs, {
            |$lhs_i: &$lhs, $rhs_i: &$rhs| -> $out { $body }(lhs, rhs)
        });
    };
    ($op:tt |$lhs_i:ident : &$lhs:ty, $rhs_i:ident : $rhs:ty| -> $out:ty $body:block) => {
        $crate::_parse_binary_op!($op, &$lhs, $rhs, $out, lhs, rhs, {
            |$lhs_i: &$lhs, $rhs_i: $rhs| -> $out { $body }(lhs, rhs)
        });
    };
    ($op:tt |$lhs_i:ident : $lhs:ty| -> $out:ty $body:block) => {
        $crate::_parse_unary_op!($op, $lhs, $out, lhs, {
            |$lhs_i: $lhs| -> $out { $body }(lhs)
        });
    };
    ($op:tt |$lhs_i:ident : $lhs:ty, $rhs_i:ident : &$rhs:ty| -> $out:ty $body:block) => {
        $crate::_parse_binary_op!($op, $lhs, &$rhs, $out, lhs, rhs, {
            |$lhs_i: $lhs, $rhs_i: &$rhs| -> $out { $body }(lhs, rhs)
        });
    };
    ($op:tt |$lhs_i:ident : $lhs:ty, $rhs_i:ident : $rhs:ty| -> $out:ty $body:block) => {
        $crate::_parse_binary_op!($op, $lhs, $rhs, $out, lhs, rhs, {
            |$lhs_i: $lhs, $rhs_i: $rhs| -> $out { $body }(lhs, rhs)
        });
    };
}

/// Overloads an operator using the given closure as its body. Generates overloads for both owned and borrowed variants where possible.
///
/// Used with the same syntax as `impl_op!` (see the [module level documentation](index.html) for more information).
///
/// Expands any borrowed inputs into both owned and borrowed variants.
///
/// `impl_op_ex!(op |a: &LHS, b: &RHS| -> OUT {...});`
/// gets expanded to
///
/// ```compile_fail
/// impl_op!(op |a: LHS, b: RHS| -> OUT {...});
/// impl_op!(op |a: LHS, b: &RHS| -> OUT {...});
/// impl_op!(op |a: &LHS, b: RHS| -> OUT {...});
/// impl_op!(op |a: &LHS, b: &RHS| -> OUT {...});
/// ```
///
/// and `impl_op_ex!(op |a: &LHS, b: RHS| -> OUT {...});`
/// gets expanded to
///
/// ```compile_fail
/// impl_op!(op |a: LHS, b: RHS| -> OUT {...});
/// impl_op!(op |a: &LHS, b: RHS| -> OUT {...});
/// ```
/// # Examples
/// ```
/// use auto_ops::impl_op_ex;
/// # #[derive(Clone, Debug, PartialEq)]
/// # struct DonkeyKong {
/// #     pub bananas: i32,
/// # }
/// # impl DonkeyKong {
/// #     pub fn new(bananas: i32) -> DonkeyKong {
/// #         DonkeyKong { bananas: bananas }
/// #     }
/// #  }
///
/// impl_op_ex!(+ |a: &DonkeyKong, b: &DonkeyKong| -> i32 { a.bananas + b.bananas });
///
/// fn main() {
///     let total_bananas = &DonkeyKong::new(2) + &DonkeyKong::new(4);
///     assert_eq!(6, total_bananas);
///     let total_bananas = &DonkeyKong::new(2) + DonkeyKong::new(4);
///     assert_eq!(6, total_bananas);
///     let total_bananas = DonkeyKong::new(2) + &DonkeyKong::new(4);
///     assert_eq!(6, total_bananas);
///     let total_bananas = DonkeyKong::new(2) + DonkeyKong::new(4);
///     assert_eq!(6, total_bananas);
/// }
#[macro_export]
macro_rules! impl_op_ex {
    ($op:tt |$lhs_i:ident : &mut $lhs:ty, $rhs_i:ident : &$rhs:ty| $body:block) => (
        $crate::_parse_assignment_op!($op, $lhs, &$rhs, lhs, rhs, {|$lhs_i : &mut $lhs, $rhs_i : &$rhs| -> () {$body} (lhs, rhs);});
        $crate::_parse_assignment_op!($op, $lhs, $rhs, lhs, rhs, {|$lhs_i : &mut $lhs, $rhs_i : &$rhs| -> () {$body} (lhs, &rhs);});
    );
    ($op:tt |$lhs_i:ident : &mut $lhs:ty, $rhs_i:ident : $rhs:ty| $body:block) => (
        $crate::_parse_assignment_op!($op, $lhs, $rhs, lhs, rhs, {|$lhs_i : &mut $lhs, $rhs_i : $rhs| -> () {$body} (lhs, rhs);});
    );
    ($op:tt |$lhs_i:ident : &$lhs:ty| -> $out:ty $body:block) => (
        $crate::_parse_unary_op!($op, &$lhs, $out, lhs, {|$lhs_i : &$lhs| -> $out {$body} (lhs)});
        $crate::_parse_unary_op!($op, $lhs, $out, lhs, {|$lhs_i : &$lhs| -> $out {$body} (&lhs)});
    );
    ($op:tt |$lhs_i:ident : &$lhs:ty, $rhs_i:ident : &$rhs:ty| -> $out:ty $body:block) => (
        $crate::impl_op!($op |$lhs_i : &$lhs, $rhs_i : &$rhs| -> $out $body);
        $crate::_parse_binary_op!($op, &$lhs, $rhs, $out, lhs, rhs, {|$lhs_i : &$lhs, $rhs_i : &$rhs| -> $out {$body} (lhs, &rhs)});
        $crate::_parse_binary_op!($op, $lhs, &$rhs, $out, lhs, rhs, {|$lhs_i : &$lhs, $rhs_i : &$rhs| -> $out {$body} (&lhs, rhs)});
        $crate::_parse_binary_op!($op, $lhs, $rhs, $out, lhs, rhs, {|$lhs_i : &$lhs, $rhs_i : &$rhs| -> $out {$body} (&lhs, &rhs)});
    );
    ($op:tt |$lhs_i:ident : &$lhs:ty, $rhs_i:ident : $rhs:ty| -> $out:ty $body:block) => (
        $crate::impl_op!($op |$lhs_i : &$lhs, $rhs_i : $rhs| -> $out $body);
        $crate::_parse_binary_op!($op, $lhs, $rhs, $out, lhs, rhs, {|$lhs_i : &$lhs, $rhs_i : $rhs| -> $out {$body} (&lhs, rhs)});
    );
    ($op:tt |$lhs_i:ident : $lhs:ty|  -> $out:ty $body:block) => (
        $crate::_parse_unary_op!($op, $lhs, $out, lhs, {|$lhs_i : $lhs| -> $out {$body} (lhs)});
    );
    ($op:tt |$lhs_i:ident : $lhs:ty, $rhs_i:ident : &$rhs:ty| -> $out:ty $body:block) => (
        $crate::impl_op!($op |$lhs_i : $lhs, $rhs_i : &$rhs| -> $out $body);
        $crate::_parse_binary_op!($op, $lhs, $rhs, $out, lhs, rhs, {|$lhs_i : $lhs, $rhs_i : &$rhs| -> $out {$body} (lhs, &rhs)});
    );
    ($op:tt |$lhs_i:ident : $lhs:ty, $rhs_i:ident : $rhs:ty| -> $out:ty $body:block) => (
        $crate::impl_op!($op |$lhs_i : $lhs, $rhs_i : $rhs| -> $out $body);
    );
}

/// Overloads a binary operator commutatively using the given closure as its body.
///
/// Used with the same syntax as `impl_op!` (see the [module level documentation](index.html) for more information).
/// Can only be used with binary operators, and the operation must be between two different types.
///
/// An operator is commutative if A <op> B == B <op> A. Common commutative operators are `+` and `*`.
///
/// ```compile_fail
/// impl_op_commutative!(op |a: LHS, b: RHS| -> OUT {...});
/// // where LHS != RHS
/// ```
///
/// gets expanded to
///
/// ```compile_fail
/// impl_op!(op |a: LHS, b: RHS| -> OUT {...});
/// impl_op!(op |a: RHS, b: LHS| -> OUT {...});
/// ```
/// Make sure that LHS != RHS, and that the operator you are trying to overload is a commutative one.
/// See the examples for what happens when you try `impl_op_commutative!` on the `-` operator (which isn't usually commutative).
/// # Examples
/// ```
/// use auto_ops::impl_op_commutative;
/// # #[derive(Clone, Debug, PartialEq)]
/// # struct DonkeyKong {
/// #     pub bananas: i32,
/// # }
/// # impl DonkeyKong {
/// #     pub fn new(bananas: i32) -> DonkeyKong {
/// #         DonkeyKong { bananas: bananas }
/// #     }
/// #  }
///
/// impl_op_commutative!(+ |a: DonkeyKong, b: i32| -> i32 { a.bananas + b });
/// // Don't do this unless you know what you are doing:
/// impl_op_commutative!(- |a: DonkeyKong, b: i32| -> i32 { a.bananas - b });
///
/// fn main() {
///     let total_bananas = DonkeyKong::new(5) + 1;
///     assert_eq!(6, total_bananas);
///     let total_bananas = 1 + DonkeyKong::new(5);
///     assert_eq!(6, total_bananas);
///     let total_bananas = DonkeyKong::new(5) - 1;
///     assert_eq!(4, total_bananas);
///     let total_bananas = 1 - DonkeyKong::new(5);
///     assert_eq!(4, total_bananas);
///     // notice that in this case (5 - 1 == 4) and (1 - 5 == 1): that is the definition of a
///     // commutative operator, but probably not what you want for the '-' operator
/// }
#[macro_export]
macro_rules! impl_op_commutative {
    ($op:tt |$lhs_i:ident : &$lhs:ty, $rhs_i:ident : &$rhs:ty| -> $out:ty $body:block) => (
        $crate::impl_op!($op |$lhs_i : &$lhs, $rhs_i : &$rhs| -> $out $body);
        $crate::_parse_binary_op!($op, &$rhs, &$lhs, $out, lhs, rhs, {|$lhs_i : &$lhs, $rhs_i : &$rhs| -> $out {$body} (rhs, lhs)});
    );
    ($op:tt |$lhs_i:ident : &$lhs:ty, $rhs_i:ident : $rhs:ty| -> $out:ty $body:block) => (
        $crate::impl_op!($op |$lhs_i : &$lhs, $rhs_i : $rhs| -> $out $body);
        $crate::_parse_binary_op!($op, $rhs, &$lhs, $out, lhs, rhs, {|$lhs_i : &$lhs, $rhs_i : $rhs| -> $out {$body} (rhs, lhs)});
    );
    ($op:tt |$lhs_i:ident : $lhs:ty, $rhs_i:ident : &$rhs:ty| -> $out:ty $body:block) => (
        $crate::impl_op!($op |$lhs_i : $lhs, $rhs_i : &$rhs| -> $out $body);
        $crate::_parse_binary_op!($op, &$rhs, $lhs, $out, lhs, rhs, {|$lhs_i : $lhs, $rhs_i : &$rhs| -> $out {$body} (rhs, lhs)});
    );
    ($op:tt |$lhs_i:ident : $lhs:ty, $rhs_i:ident : $rhs:ty| -> $out:ty $body:block) => (
        $crate::impl_op!($op |$lhs_i : $lhs, $rhs_i : $rhs| -> $out $body);
        $crate::_parse_binary_op!($op, $rhs, $lhs, $out, lhs, rhs, {|$lhs_i : $lhs, $rhs_i : $rhs| -> $out {$body} (rhs, lhs)});
    );
}

/// Overloads a binary operator commutatively using the given closure as its body. Generates overloads for both owned and borrowed variants where possible.
///
/// See [`impl_op_commutative!`](macro.impl_op_commutative.html) for usage.
///
/// Expands borrowed inputs to both borrowed and owned variants.
///
/// ```compile_fail
/// impl_op_ex_commutative!(op |a: &LHS, b: &RHS| -> OUT {...});
/// // where LHS != RHS
/// ```
///
/// gets expanded to
///
/// ```compile_fail
/// impl_op!(op |a: &LHS, b: &RHS| -> OUT {...});
/// impl_op!(op |a: &LHS, b: RHS| -> OUT {...});
/// impl_op!(op |a: LHS, b: &RHS| -> OUT {...});
/// impl_op!(op |a: LHS, b: RHS| -> OUT {...});
///
/// impl_op!(op |a: &RHS, b: &LHS| -> OUT {...});
/// impl_op!(op |a: &RHS, b: LHS| -> OUT {...});
/// impl_op!(op |a: RHS, b: &LHS| -> OUT {...});
/// impl_op!(op |a: RHS, b: LHS| -> OUT {...});
/// ```
/// # Examples
/// ```
/// use auto_ops::impl_op_ex_commutative;
/// # #[derive(Clone, Debug, PartialEq)]
/// # struct DonkeyKong {
/// #     pub bananas: i32,
/// # }
/// # impl DonkeyKong {
/// #     pub fn new(bananas: i32) -> DonkeyKong {
/// #         DonkeyKong { bananas: bananas }
/// #     }
/// #  }
/// # #[derive(Clone, Debug, PartialEq)]
/// # struct DiddyKong {
/// #     pub bananas: i32,
/// # }
/// # impl DiddyKong {
/// #     pub fn new(bananas: i32) -> DiddyKong {
/// #         DiddyKong { bananas: bananas }
/// #     }
/// #  }
///
/// impl_op_ex_commutative!(+ |a: &DonkeyKong, b: &DiddyKong| -> i32 { a.bananas + b.bananas });
/// impl_op_ex_commutative!(+ |a: &DonkeyKong, b: i32| -> i32 { a.bananas + b });
///
/// fn main() {
///     let total_bananas = &DonkeyKong::new(5) + &DiddyKong::new(1);
///     assert_eq!(6, total_bananas);
///     let total_bananas = &DonkeyKong::new(5) + DiddyKong::new(1);
///     assert_eq!(6, total_bananas);
///     let total_bananas = DonkeyKong::new(5) + &DiddyKong::new(1);
///     assert_eq!(6, total_bananas);
///     let total_bananas = DonkeyKong::new(5) + DiddyKong::new(1);
///     assert_eq!(6, total_bananas);
///
///     let total_bananas = &DiddyKong::new(1) + &DonkeyKong::new(5);
///     assert_eq!(6, total_bananas);
///     let total_bananas = &DiddyKong::new(1) + DonkeyKong::new(5);
///     assert_eq!(6, total_bananas);
///     let total_bananas = DiddyKong::new(1) + &DonkeyKong::new(5);
///     assert_eq!(6, total_bananas);
///     let total_bananas = DiddyKong::new(1) + DonkeyKong::new(5);
///     assert_eq!(6, total_bananas);
///
///     let total_bananas = &DonkeyKong::new(5) + 1;
///     assert_eq!(6, total_bananas);
///     let total_bananas = DonkeyKong::new(5) + 1;
///     assert_eq!(6, total_bananas);
///
///     let total_bananas = 1 + &DonkeyKong::new(5);
///     assert_eq!(6, total_bananas);
///     let total_bananas = 1 + DonkeyKong::new(5);
///     assert_eq!(6, total_bananas);
/// }
#[macro_export]
macro_rules! impl_op_ex_commutative {
    ($op:tt |$lhs_i:ident : &$lhs:ty, $rhs_i:ident : &$rhs:ty| -> $out:ty $body:block) => (
        $crate::impl_op_ex!($op |$lhs_i : &$lhs, $rhs_i : &$rhs| -> $out $body);

        $crate::_parse_binary_op!($op, &$rhs, &$lhs, $out, lhs, rhs, {|$lhs_i : &$lhs, $rhs_i : &$rhs| -> $out {$body} (rhs, lhs)});
        $crate::_parse_binary_op!($op, &$rhs, $lhs, $out, lhs, rhs, {|$lhs_i : &$lhs, $rhs_i : &$rhs| -> $out {$body} (&rhs, lhs)});
        $crate::_parse_binary_op!($op, $rhs, &$lhs, $out, lhs, rhs, {|$lhs_i : &$lhs, $rhs_i : &$rhs| -> $out {$body} (rhs, &lhs)});
        $crate::_parse_binary_op!($op, $rhs, $lhs, $out, lhs, rhs, {|$lhs_i : &$lhs, $rhs_i : &$rhs| -> $out {$body} (&rhs, &lhs)});
    );
    ($op:tt |$lhs_i:ident : &$lhs:ty, $rhs_i:ident : $rhs:ty| -> $out:ty $body:block) => (
        $crate::impl_op_ex!($op |$lhs_i : &$lhs, $rhs_i : $rhs| -> $out $body);

        $crate::_parse_binary_op!($op, $rhs, &$lhs, $out, lhs, rhs, {|$lhs_i : &$lhs, $rhs_i : $rhs| -> $out {$body} (rhs, lhs)});
        $crate::_parse_binary_op!($op, $rhs, $lhs, $out, lhs, rhs, {|$lhs_i : &$lhs, $rhs_i : $rhs| -> $out {$body} (&rhs, lhs)});
    );
    ($op:tt |$lhs_i:ident : $lhs:ty, $rhs_i:ident : &$rhs:ty| -> $out:ty $body:block) => (
        $crate::impl_op_ex!($op |$lhs_i : $lhs, $rhs_i : &$rhs| -> $out $body);

        $crate::_parse_binary_op!($op, &$rhs, $lhs, $out, lhs, rhs, {|$lhs_i : $lhs, $rhs_i : &$rhs| -> $out {$body} (rhs, lhs)});
        $crate::_parse_binary_op!($op, $rhs, $lhs, $out, lhs, rhs, {|$lhs_i : $lhs, $rhs_i : &$rhs| -> $out {$body} (rhs, &lhs)});
    );
    ($op:tt |$lhs_i:ident : $lhs:ty, $rhs_i:ident : $rhs:ty| -> $out:ty $body:block) => (
        $crate::impl_op_commutative!($op |$lhs_i : $lhs, $rhs_i : $rhs| -> $out $body);
    );
}
