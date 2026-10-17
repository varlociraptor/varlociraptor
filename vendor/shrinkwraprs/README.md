# shrinkwraprs [![pipeline status](https://gitlab.com/williamyaoh/shrinkwraprs/badges/master/pipeline.svg)](https://gitlab.com/williamyaoh/shrinkwraprs/commits/master)

[![latest version](https://img.shields.io/crates/v/shrinkwraprs.svg)](https://crates.io/crates/shrinkwraprs)
[![api documentation](https://docs.rs/shrinkwraprs/badge.svg)](https://docs.rs/shrinkwraprs)
[![license](https://img.shields.io/badge/license-BSD--3-ff69b4.svg)](https://gitlab.com/williamyaoh/shrinkwraprs/blob/master/LICENSE)

Making wrapper types allows us to give more compile-time
guarantees about our code being correct:

```rust
// Now we can't mix up widths and heights; the compiler will yell at us!
struct Width(u64);
struct Height(u64);
```

But... they're kind of a pain to work with. If you ever need to get at
that wrapped `u64`, you need to constantly pattern-match back and forth
to wrap and unwrap the values.

`shrinkwraprs` aims to alleviate this pain by allowing you to derive
implementations of various conversion traits by deriving
`Shrinkwrap`.

## Functionality implemented

Currently, using `#[derive(Shrinkwrap)]` will derive the following traits
for all structs:

* `AsRef<InnerType>`
* `Borrow<InnerType>`
* `Deref<Target=InnerType>`

Additionally, using `#[shrinkwrap(mutable)]` will also
derive the following traits:

* `AsMut<InnerType>`
* `BorrowMut<InnerType>`
* `DerefMut<Target=InnerType>`

Finally, one more option is `#[shrinkwrap(transformers)]`, which will derive
some useful inherent functions for transforming the wrapped data:

* `fn transform<F>(&mut self, mut f: F) -> &mut Self where F: FnMut(&mut InnerType)`
* `fn siphon<F, T>(self, mut f: F) -> T where F: FnMut(InnerType) -> T`

...where `transform` makes it easy to chain updates on the inner value, and
`siphon` allows you to easily move out the inner value to produce a value
of a different type.

`transform` will have the same visibility as the inner field, which ensures that
`transform` doesn't leak the possibility of changing the inner value
(potentially in invariant-violating ways). `siphon` has the same visibility as
the struct itself, since it *doesn't* provide a direct way for callers to break
your data.

## Cool, how do I use it?

First, add `shrinkwraprs` as a dependency in your `Cargo.toml`:

```toml
[dependencies]

shrinkwraprs = "0.3.0"
```

Then, just slap a `#[derive(Shrinkwrap)]` on any structs you want
convenient-ified:

```rust
#[macro_use] extern crate shrinkwraprs;

#[derive(Shrinkwrap)]
struct Email(String);

fn main() {
    let email = Email("chiya+snacks@natsumeya.jp".into());

    let is_discriminated_email =
        email.contains("+");  // Woohoo, we can use the email like a string!

    /* ... */
}
```

If you have multiple fields, but there's only one field you want to be able
to deref/borrow as, mark it with `#[shrinkwrap(main_field)]`:

```rust
#[derive(Shrinkwrap)]
struct Email {
    spamminess: f64,
    #[shrinkwrap(main_field)] addr: String
}

#[derive(Shrinkwrap)]
struct CodeSpan(u32, u32, #[shrinkwrap(main_field)] Token);
```

If you also want to be able to modify the wrapped value directly,
add the attribute `#[shrinkwrap(mutable)]` as well:

```rust
#[derive(Shrinkwrap)]
#[shrinkwrap(mutable)]
struct InputBuffer {
    buffer: String
}

...
let mut input_buffer = /* ... */;
input_buffer.push_str("some values");
...
```
