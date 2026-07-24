//! We want to make sure that the struct that our caller passes us is in the
//! right form. However, we don't want to clutter up our code generation
//! logic with lots of error handling. So instead, we take in our `DeriveInput`
//! and do all the error handling in one place, transforming it into an AST
//! specific to our crate if it's valid.

use syn;
use syn::spanned::Spanned;

use proc_macro2::{Span, TokenStream};

use itertools::Itertools;

type Fields = Vec<syn::Field>;

bitflags! {
  /// Controls which code and implementations we generate.
  pub struct ShrinkwrapFlags: u32 {
    const SW_MUT          = 0b00000001;
    const SW_IGNORE_VIS   = 0b00000010;
    const SW_TRANSFORMERS = 0b00000100;
  }
}

pub struct StructDetails {
  pub flags: ShrinkwrapFlags,
  pub ident: syn::Ident,
  pub generics: syn::Generics,
  pub visibility: syn::Visibility
}

/// Represents either a tuple or bracketed struct with at least one field.
pub struct Struct {
  pub inner_field: syn::Field,
  pub inner_field_name: proc_macro2::TokenStream
}

/// Check if the input stream matches our required data structures.
/// The TokenStream on error contains a compile error pointing to the right place.
pub fn validate_derive_input(input: syn::DeriveInput)
  -> Result<(StructDetails, Struct), TokenStream>
{
  // We *don't* want to use `panic` and `unwrap` here, even though they're
  // safe, because we want our compile errors to be attached to the right
  // lines of code.

  use syn::{DeriveInput, DataStruct, FieldsUnnamed, FieldsNamed};
  use syn::Data::{Struct, Enum, Union};
  use syn::Fields::{Named, Unnamed};

  let whole_span = input.span();
  let DeriveInput { attrs, vis, ident, generics, data, .. } = input;

  let flags = shrinkwrap_flags(&attrs);
  let details = StructDetails { flags, ident, visibility: vis, generics };

  let strct = match data {
    Struct(DataStruct { fields: Unnamed(FieldsUnnamed { unnamed: fields, .. }), .. }) => {
      let fields = fields.into_iter().collect_vec();
      validate_tuple(whole_span, fields)
    },
    Struct(DataStruct { fields: Named(FieldsNamed { named: fields, .. }), .. }) => {
      let fields = fields.into_iter().collect_vec();
      validate_nontuple(whole_span, fields)
    },
    Struct(..) =>
      Err(compile_error_at(whole_span, "Shrinkwrap needs a struct with at least one field!")),
    Enum(..) =>
      Err(compile_error_at(whole_span, "Shrinkwrap does not support enums!")),
    Union(..) =>
      Err(compile_error_at(whole_span, "Shrinkwrap does not support C-style unions!"))
  }?;

  Ok((details, strct))
}

/// Specifically for working with attributes like #[shrinkwrap(..)], where
/// a name is combined with a list of attributes. Get the list of attributes
/// matching the tag.
fn tagged_attrs(tag: &str, attrs: &[syn::Attribute]) -> Vec<syn::NestedMeta> {
  use syn::{Meta, MetaList};

  let mut result = vec![];

  for attr in attrs {
    let meta = attr.parse_meta();

    if let Ok(Meta::List(MetaList { path, nested, .. })) = meta {
      if path.is_ident(tag) {
        result.extend(nested);
      }
    }
  }

  result
}

fn shrinkwrap_flags(attrs: &[syn::Attribute]) -> ShrinkwrapFlags {
  use syn::{Meta, NestedMeta};

  let meta = tagged_attrs("shrinkwrap", attrs);
  let mut flags = ShrinkwrapFlags::empty();

  for attr in meta {
    if let NestedMeta::Meta(Meta::Path(path)) = attr {
      if path.is_ident("mutable") {
        flags |= ShrinkwrapFlags::SW_MUT;
      } else if path.is_ident("unsafe_ignore_visibility") {
        flags |= ShrinkwrapFlags::SW_IGNORE_VIS;
      } else if path.is_ident("transformers") {
        flags |= ShrinkwrapFlags::SW_TRANSFORMERS;
      }
    }
  }

  flags
}

fn is_marked(field: &syn::Field) -> bool {
  use syn::{Meta, NestedMeta};

  let meta = tagged_attrs("shrinkwrap", &field.attrs);

  meta.into_iter().any(|meta| {
    if let NestedMeta::Meta(Meta::Path(path)) = meta {
      path.is_ident("main_field")
    } else {
      false
    }
  })
}

/// Only a single field, out of all a struct's fields, can be marked as
/// the main field that we deref to. So let's find that field.
/// We also return the 0-based number of the marked field.
fn find_marked_field(whole_span: Span, fields: Fields)
  -> Result<((usize, syn::Field), Fields), TokenStream>
{
  let (marked, unmarked) = fields.into_iter()
    .enumerate()
    .partition::<Vec<_>, _>(|&(_, ref field)| is_marked(field));
  let marked_len = marked.len();
  let single: Option<(_,)> = marked.into_iter()
    .collect_tuple();

  match (single, unmarked.len()) {
    (Some((field,)), _) => {
      let unmarked = unmarked.into_iter()
        .map(|(_, field)| field)
        .collect_vec();

      Ok((field, unmarked))
    }
    (None, 1) => {
      let single: (_,) = unmarked.into_iter()
        .collect_tuple()
        .unwrap();

      Ok((single.0, vec![]))
    },
    _ => if marked_len == 0 {
      Err(compile_error_at(
        whole_span,
        "Shrinkwrap doesn't know which field you want this struct to convert
to. Did you forget to mark a field with #[shrinkwrap(main_field)]?"
      ))
    } else {
      Err(compile_error_at(
        whole_span,
        "Shrinkwrap sees too many marked fields in this struct. Did you
accidentally mark more than one field with #[shrinkwrap(main_field)]?"
      ))
    }
  }
}

fn validate_tuple(whole_span: Span, fields: Fields) -> Result<Struct, TokenStream> {
  if fields.len() == 0 {
    return Err(compile_error_at(
      whole_span,
      "Shrinkwrap requires tuple structs to have at least one field!"
    ));
  }

  let ((marked_index, marked_field), _) = find_marked_field(whole_span, fields)?;
  let index: syn::Index = marked_index.into();

  Ok(Struct {
    inner_field: marked_field,
    inner_field_name: quote!( #index )
  })
}

fn validate_nontuple(whole_span: Span, fields: Fields) -> Result<Struct, TokenStream> {
  if fields.len() == 0 {
    return Err(compile_error_at(
      whole_span,
      "Shrinkwrap requires structs to have at least one field!"
    ));
  }

  let ((_, marked_field), _) = find_marked_field(whole_span, fields)?;
  let ident = marked_field.ident
    .clone()
    .unwrap();

  Ok(Struct {
    inner_field: marked_field,
    inner_field_name: quote!( #ident )
  })
}

fn compile_error_at(at: Span, msg: &str) -> TokenStream {
  quote_spanned!(at=> compile_error!(#msg);)
}

#[cfg(test)]
mod tests {
  use syn;
  use itertools::Itertools;

  use super::*;

  #[test]
  fn test_field_attribute_found() {
    let input = r"
      struct Foo {
        field1: u32,
        #[shrinkwrap(main_field)]
        field2: u32
      }
    ";

    let strct: syn::DeriveInput = syn::parse_str(input)
      .unwrap();

    match strct.data {
      syn::Data::Struct(syn::DataStruct { fields, .. }) => {
        let marked = fields.into_iter()
          .filter(|field| is_marked(field));
        let field: (syn::Field,) = marked
          .collect_tuple()
          .unwrap();
        let ident = field.0.ident
          .unwrap();

        assert_eq!(&ident, "field2");
      },
      _ => panic!()
    }
  }

  #[test]
  fn test_field_attribute_not_found() {
    let input = r"
      struct Foo {
        field1: u32,
        field2: u32
      }
    ";

    let strct: syn::DeriveInput = syn::parse_str(input)
      .unwrap();

    match strct.data {
      syn::Data::Struct(syn::DataStruct { fields, .. }) => {
        let marked = fields.into_iter()
          .filter(|field| is_marked(field))
          .collect_vec();
        assert_eq!(marked.len(), 0);
      },
      _ => panic!()
    }
  }
}
