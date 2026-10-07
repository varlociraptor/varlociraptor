//! A procedural macro that allows you to derive `std::ops::Deref` and `std::ops::DerefMut` for
//! your structs. This macro can only be derived on structs **with atleast one field**. You can
//! specify which field you want to be deref'ed to with the `#[deref]` and allow mutable
//! dereferencing with `#[deref(mutable)]`.
//!
//! ## Deriving `std::ops::Deref`
//! ```ignore
//! use std::collections::HashMap;
//!
//! use derefable::Derefable;
//!
//! #[derive(Default, Derefable)]
//! struct Map {
//!     #[deref]
//!     inner: HashMap<&'static str, &'static str>
//! }
//!
//! fn main() {
//!     let map = Map::default();
//!
//!     assert!(map.is_empty());
//! }
//! ```
//!
//! ## Deriving `std::ops::DerefMut`
//! ```ignore
//! use std::collections::HashMap;
//!
//! use derefable::Derefable;
//!
//! #[derive(Default, Derefable)]
//! struct MutableMap {
//!     #[deref(mutable)]
//!     inner: HashMap<&'static str, &'static str>
//! }
//!
//! fn main() {
//!     let mut map = MutableMap::default();
//!
//!     map.insert("Hello", "World");
//!
//!     assert_eq!(map.get("Hello"), Some("World"));
//! }
//! ```

extern crate proc_macro;
extern crate proc_macro2;

use proc_macro::TokenStream;
use proc_macro2::{Ident, Span};
use quote::{quote, ToTokens};
use syn::*;

#[proc_macro_derive(Derefable, attributes(deref))]
pub fn derefable_derive(input: TokenStream) -> TokenStream {
    let ast = syn::parse(input).unwrap();

    impl_derefable(&ast)
}

fn impl_derefable(ast: &syn::DeriveInput) -> TokenStream {
    let name = &ast.ident;

    let data = match ast.data {
        Data::Struct(ref s) => s,
        _ => {
            // name.span()
            //     .unstable()
            //     .error("`#[derive(Derefable)]` is only available for structs")
            //     .emit();

            // return TokenStream::new()
            panic!("`#[derive(Derefable)]` is only available for structs")
        }
    };

    match data.fields {
        Fields::Unit => {
            // data.span()
            //     .unstable()
            //     .error("`#[derive(Derefable)]` requires a field to be able to deref")
            //     .emit();

            // return TokenStream::new()
            panic!("`#[derive(Derefable)]` requires a field to be able to deref")
        }
        _ => {}
    }

    let mut deref_field = None;
    let mut is_field_mutable = false;

    for (i, field) in data.fields.iter().enumerate() {
        for attribute in &field.attrs {
            if let Ok(meta) = attribute.parse_meta() {
                match meta {
                    Meta::Word(ref ident) if ident == "deref" => {
                        if deref_field.is_none() {
                            deref_field = Some((field.clone(), i as u32));
                        } else {
                            // name.span()
                            //     .unstable()
                            //     .error("Only one field in a struct can be `#[deref]`")
                            //     .emit();

                            // return TokenStream::new()
                            panic!("Only one field in a struct can be `#[deref]`")
                        }
                    }

                    Meta::List(MetaList {
                        ref ident,
                        ref nested,
                        ..
                    }) if ident == "deref" => {
                        is_field_mutable = nested.iter().any(|nested_item| match nested_item {
                            NestedMeta::Meta(m) => match m {
                                Meta::Word(ident) => ident == "mutable",
                                _ => false,
                            },
                            _ => false,
                        });

                        if deref_field.is_none() {
                            deref_field = Some((field.clone(), i as u32));
                        } else {
                            // name.span()
                            //     .unstable()
                            //     .error("Only one field in a struct can be `#[deref]`")
                            //     .emit();

                            // return TokenStream::new()
                            panic!("Only one field in a struct can be `#[deref]`")
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    if deref_field.is_none() {
        // name.span()
        //     .unstable()
        //     .error("`#[derive(Derefable)]` requires one field to be marked `#[deref]`")
        //     .emit();

        // return TokenStream::new()
        panic!("`#[derive(Derefable)]` requires one field to be marked `#[deref]`");
    }

    let (field, index) = deref_field.unwrap();

    let target = field.ty;
    let ident = field
        .ident
        .map(Ident::into_token_stream)
        .unwrap_or_else(|| {
            Index {
                index,
                span: Span::call_site(),
            }
            .into_token_stream()
        });

    let mut_gen = if is_field_mutable {
        quote! {
            impl std::ops::DerefMut for #name {
                fn deref_mut(&mut self) -> &mut Self::Target {
                    &mut self.#ident
                }
            }
        }
    } else {
        quote! {}
    };

    let gen = quote! {
        impl std::ops::Deref for #name {
            type Target = #target;

            fn deref(&self) -> &Self::Target {
                &self.#ident
            }
        }

        #mut_gen
    };

    gen.into()
}
