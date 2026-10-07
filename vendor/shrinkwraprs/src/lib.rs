//! # shrinkwraprs
//!
//! Making wrapper types allows us to give more compile-time
//! guarantees about our code being correct:
//!
//! ```ignore
//! // Now we can't mix up widths and heights; the compiler will yell at us!
//! struct Width(u64);
//! struct Height(u64);
//! ```
//!
//! But... they're kind of a pain to work with. If you ever need to get at
//! that wrapped `u64`, you need to constantly pattern-match back and forth
//! to wrap and unwrap the values.
//!
//! `shrinkwraprs` aims to alleviate this pain by allowing you to derive
//! implementations of various conversion traits by deriving
//! `Shrinkwrap`.
//!
//! ## Functionality implemented
//!
//! Currently, using `#[derive(Shrinkwrap)]` will derive the following traits
//! for all structs:
//!
//! * `AsRef<InnerType>`
//! * `Borrow<InnerType>`
//! * `Deref<Target=InnerType>`
//!
//! Additionally, using `#[shrinkwrap(mutable)]` will also
//! derive the following traits:
//!
//! * `AsMut<InnerType>`
//! * `BorrowMut<InnerType>`
//! * `DerefMut<Target=InnerType>`
//!
//! Finally, one more option is `#[shrinkwrap(transformers)]`, which will derive
//! some useful inherent functions for transforming the wrapped data:
//!
//! * `fn transform<F>(&mut self, mut f: F) -> &mut Self where F: FnMut(&mut InnerType)`
//! * `fn siphon<F, T>(self, mut f: F) -> T where F: FnMut(InnerType) -> T`
//!
//! ...where `transform` makes it easy to chain updates on the inner value, and
//! `siphon` allows you to easily move out the inner value to produce a value
//! of a different type.
//!
//! `transform` will have the same visibility as the inner field, which ensures that
//! `transform` doesn't leak the possibility of changing the inner value
//! (potentially in invariant-violating ways). `siphon` has the same visibility as
//! the struct itself, since it *doesn't* provide a direct way for callers to break
//! your data.
//!
//! ## Cool, how do I use it?
//!
//! ```ignore
//! #[macro_use] extern crate shrinkwraprs;
//!
//! #[derive(Shrinkwrap)]
//! struct Email(String);
//!
//! fn main() {
//!     let email = Email("chiya+snacks@natsumeya.jp".into());
//!
//!     let is_discriminated_email =
//!         email.contains("+");  // Woohoo, we can use the email like a string!
//!
//!     /* ... */
//! }
//! ```
//!
//! If you have multiple fields, but there's only one field you want to be able
//! to deref/borrow as, mark it with `#[shrinkwrap(main_field)]`:
//!
//! ```ignore
//! #[derive(Shrinkwrap)]
//! struct Email {
//!     spamminess: f64,
//!     #[shrinkwrap(main_field)] addr: String
//! }
//!
//! #[derive(Shrinkwrap)]
//! struct CodeSpan(u32, u32, #[shrinkwrap(main_field)] Token);
//! ```
//!
//! If you also want to be able to modify the wrapped value directly,
//! add the attribute `#[shrinkwrap(mutable)]` as well:
//!
//! ```ignore
//! #[derive(Shrinkwrap)]
//! #[shrinkwrap(mutable)]
//! struct InputBuffer {
//!     buffer: String
//! }
//!
//! ...
//! let mut input_buffer = /* ... */;
//! input_buffer.push_str("some values");
//! ...
//! ```

// Additionally, perhaps subsume some functionality from
// [`from_variants`](https://crates.io/crates/from_variants)?

#![cfg_attr(feature = "strict", deny(warnings))]
#![recursion_limit="128"]

extern crate proc_macro;
extern crate proc_macro2;
extern crate syn;
#[macro_use] extern crate quote;
extern crate itertools;
#[macro_use] extern crate bitflags;

use proc_macro2::{TokenStream, Span};
use quote::ToTokens;
use syn::spanned::Spanned;

mod ast;
mod visibility;

#[proc_macro_derive(Shrinkwrap, attributes(shrinkwrap))]
pub fn shrinkwrap(tokens: proc_macro::TokenStream) -> proc_macro::TokenStream {
  use ast::{ShrinkwrapFlags, validate_derive_input};
  use visibility::field_visibility;
  use visibility::FieldVisibility::*;

  let input: syn::DeriveInput = syn::parse(tokens)
    .unwrap();
  let validate_result = validate_derive_input(input);

  let mut tokens = TokenStream::new();

  match validate_result {
    Err(error) => error.to_tokens(&mut tokens),
    Ok((details, input)) => {
      impl_immut_borrows(&details, &input)
        .to_tokens(&mut tokens);

      if details.flags.contains(ShrinkwrapFlags::SW_TRANSFORMERS) {
        impl_transformers(&details, &input)
          .to_tokens(&mut tokens);
      }

      if details.flags.contains(ShrinkwrapFlags::SW_MUT) {
        // Make sure that the inner field isn't less visible than the outer struct.
        if !details.flags.contains(ast::ShrinkwrapFlags::SW_IGNORE_VIS) {
          match field_visibility(&details.visibility, &input.inner_field.vis) {
            Restricted => {
              let outer_vis = show_visibility(&details.visibility);
              let inner_vis = show_visibility(&input.inner_field.vis);

              let msg = format!(
"Encountered an error while mutably Shrinkwrapping this field.
This field is less visible than its containing struct; the containing
struct

    `{}'

has visibility

    `{}'

while this field has visibility

    `{}'

Implementing mutable shrinkwraps could allow outside code to modify the
inner value, even when visibility modifiers say that they can't.

Some ways to solve this problem:

1) Change the visibility of the inner field to be the same as its
   containing struct: `{}'
2) Turn off mutable shrinkwraps if you don't need them
3) Override this check and implement mutable shrinkwraps anyways by using
   #[shrinkwrap(unsafe_ignore_visibility)] on your struct",
                &details.ident,
                outer_vis,
                inner_vis,
                outer_vis
              );

              quote_spanned!(input.inner_field.span()=> compile_error!(#msg);)
                .to_tokens(&mut tokens);
            },
            CantDetermine => {
              let outer_vis = show_visibility(&details.visibility);
              let inner_vis = show_visibility(&input.inner_field.vis);

              let msg = format!(
"Encountered an error while mutably Shrinkwrapping this field.
I can't determine whether the inner field is at least as visible as its
containing struct; the containing struct

    `{}'

has visibility

    `{}'

while this field has visibility

    `{}'

If the inner field is less visible than its containing struct, implementing
mutable shrinkwraps could allow outside code to modify the inner value,
even when visibility modifiers say that they can't.

Some ways to solve this problem:

1) Change the visibility of the inner field to be the same as its
   containing struct: `{}'
2) Turn off mutable shrinkwraps if you don't need them
3) Override this check and implement mutable shrinkwraps anyways by using
   #[shrinkwrap(unsafe_ignore_visibility)] on your struct",
                &details.ident,
                outer_vis,
                inner_vis,
                outer_vis
              );

              quote_spanned!(input.inner_field.span()=> compile_error!(#msg);)
                .to_tokens(&mut tokens);
            },
            _ => ()
          }
        }

        impl_mut_borrows(&details, &input)
          .to_tokens(&mut tokens);
      }
    }
  }

  tokens.into()
}

// When generating our code, we need to be careful not to leak things into the
// surrounding code. For example, we don't use imports unless they're inside a
// scope, because otherwise we'd be inserting invisible imports whenever a user
// used #[derive(Shrinkwrap)].

fn impl_immut_borrows(details: &ast::StructDetails, input: &ast::Struct) -> TokenStream {
  let &ast::StructDetails { ref ident, ref generics, .. } = details;
  let &ast::Struct { ref inner_field, ref inner_field_name, .. } = input;

  let inner_type = &inner_field.ty;

  let (impl_generics, ty_generics, where_clause) = generics.split_for_impl();
  let rust = syn::Ident::new(RUST, Span::call_site());

  quote! {
    impl #impl_generics ::#rust::ops::Deref for #ident #ty_generics #where_clause {
      type Target = #inner_type;
      fn deref(&self) -> &Self::Target {
        &self.#inner_field_name
      }
    }

    impl #impl_generics ::#rust::borrow::Borrow<#inner_type> for #ident #ty_generics #where_clause {
      fn borrow(&self) -> &#inner_type {
        &self.#inner_field_name
      }
    }

    impl #impl_generics ::#rust::convert::AsRef<#inner_type> for #ident #ty_generics #where_clause {
      fn as_ref(&self) -> &#inner_type {
        &self.#inner_field_name
      }
    }
  }
}

fn impl_mut_borrows(details: &ast::StructDetails, input: &ast::Struct) -> TokenStream {
  let &ast::StructDetails { ref ident, ref generics, .. } = details;
  let &ast::Struct { ref inner_field, ref inner_field_name, .. } = input;

  let inner_type = &inner_field.ty;

  let (impl_generics, ty_generics, where_clause) = generics.split_for_impl();
  let rust = syn::Ident::new(RUST, Span::call_site());

  quote! {
    impl #impl_generics ::#rust::ops::DerefMut for #ident #ty_generics #where_clause {
      fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.#inner_field_name
      }
    }

    impl #impl_generics ::#rust::borrow::BorrowMut<#inner_type> for #ident #ty_generics #where_clause {
      fn borrow_mut(&mut self) -> &mut #inner_type {
        &mut self.#inner_field_name
      }
    }

    impl #impl_generics ::#rust::convert::AsMut<#inner_type> for #ident #ty_generics #where_clause {
      fn as_mut(&mut self) -> &mut #inner_type {
        &mut self.#inner_field_name
      }
    }
  }
}

fn impl_transformers(details: &ast::StructDetails, input: &ast::Struct) -> TokenStream {
  let &ast::StructDetails { ref ident, ref generics, .. } = details;
  let &ast::Struct { ref inner_field, ref inner_field_name, .. } = input;

  let inner_type = &inner_field.ty;
  let inner_visibility = &inner_field.vis;

  let (impl_generics, ty_generics, where_clause) = generics.split_for_impl();

  // This is a *massive* hack to avoid variable capture, but I can't figure out
  // how to get `quote` to enforce hygiene or generate a gensym.
  let f = quote!( __SHRINKWRAP_F );
  let t = quote!( __SHRINKWRAP_T );

  quote! {
    #[allow(dead_code, non_camel_case_types)]
    impl #impl_generics #ident #ty_generics #where_clause {
      /// Create a new value by calling a function over the wrapped value,
      /// consuming the outer value in the process.
      pub fn siphon<#t, #f: FnMut(#inner_type) -> #t>(self, mut f: #f) -> #t {
        f(self.#inner_field_name)
      }

      /// Update the outer value by calling a function on the wrapped value, which
      /// might update the wrapped value in place. Returns a mutable ref to the
      /// original outer value for easy chaining.
      #inner_visibility fn transform<#f>(&mut self, mut f: #f) -> &mut Self
        where #f: FnMut(&mut #inner_type)
      {
        f(&mut self.#inner_field_name);
        self
      }
    }
  }
}

/// Turn the visibility AST into what you might see in source code.
fn show_visibility(vis: &syn::Visibility) -> String {
  match vis {
    syn::Visibility::Inherited => "<private>".into(),
    _ => {
      let mut tokens = TokenStream::new();
      vis.to_tokens(&mut tokens);

      format!("{}", tokens)
    }
  }
}

#[cfg(feature = "std")]
const RUST: &str = "std";
#[cfg(not(feature = "std"))]
const RUST: &str = "core";
