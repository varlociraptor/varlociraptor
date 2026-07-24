# auto_ops [![Build Status]][travis] [![Latest Version]][crates.io]

[Build Status]: https://api.travis-ci.org/carbotaniuman/auto_ops.svg?branch=master
[travis]: https://travis-ci.org/carbotaniuman/auto_ops
[Latest Version]: https://img.shields.io/crates/v/auto_ops.svg
[crates.io]: https://crates.io/crates/auto_ops

Macros for easy operator overloading.

[Documentation](https://docs.rs/auto_ops/)

This library is forked from the original [impl_ops](https://github.com/brianwp3000/impl_ops) by brianwp3000.

This library makes writing multiple `impl std::ops::<op>` blocks much faster, especially when you want operators defined for both owned and borrowed variants of the inputs.

To use, import the macros with `use auto_ops::*;`. Remember that you can only overload operators between one or more types defined in the current crate (respecting Rust orphan rules).
# Examples
```rust
use auto_ops::*;

#[derive(Clone, Debug, PartialEq)]
struct DonkeyKong {
    pub bananas: i32,
}
impl DonkeyKong {
    pub fn new(bananas: i32) -> DonkeyKong {
        DonkeyKong { bananas: bananas }
    }
}

impl_op_ex!(+ |a: &DonkeyKong, b: &DonkeyKong| -> DonkeyKong { DonkeyKong::new(a.bananas + b.bananas) });
impl_op_ex!(+= |a: &mut DonkeyKong, b: &DonkeyKong| { a.bananas += b.bananas });

fn main() {
    assert_eq!(DonkeyKong::new(5), DonkeyKong::new(4) + DonkeyKong::new(1));
    assert_eq!(DonkeyKong::new(5), DonkeyKong::new(4) + &DonkeyKong::new(1));
    assert_eq!(DonkeyKong::new(5), &DonkeyKong::new(4) + DonkeyKong::new(1));
    assert_eq!(DonkeyKong::new(5), &DonkeyKong::new(4) + &DonkeyKong::new(1));

    let mut dk = DonkeyKong::new(4);
    dk += DonkeyKong::new(1);
    dk += &DonkeyKong::new(1);
    assert_eq!(DonkeyKong::new(6), dk);
}
```

# Roadmap
With Rust lifetime inference changes, implementations for generic (over types and lifetimes) impls are being worked on.
