//! Low-level manipulations of IEEE754 floating-point numbers.
//!
//! # Installation
//!
//! Add this to your Cargo.toml
//!
//! ```toml
//! [dependencies]
//! ieee754 = "0.2"
//! ```
//!
//! # Examples
//!
//! ```rust
//! use ieee754::Ieee754;
//!
//! // there are 840 single-precision floats between 1.0 and 1.0001
//! // (inclusive).
//! assert_eq!(1_f32.upto(1.0001).count(), 840);
//! ```

#![no_std]
#[cfg(test)] #[macro_use] extern crate std;

mod iter;
mod impls;
mod traits;

pub use traits::{Bits, Ieee754};
pub use iter::Iter;

#[inline]
#[doc(hidden)]
pub fn abs<F: Ieee754>(x: F) -> F {
    x.abs()
}
