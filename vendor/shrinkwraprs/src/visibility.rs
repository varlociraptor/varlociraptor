//! We want to make sure that providing mutable traits doesn't accidentally
//! leak the internal implementation details of a shrinkwrapped type.
//!
//! To do that, we need to make sure that the inner field has the same
//! visibility as the shrinkwrapped struct itself. If it doesn't, we can
//! give the user an error and refuse to generate implementations.

use syn;

use itertools::Itertools;

// When checking for visibility containment, we can make use of the guarantee
// that the langauge provides us that any visibility path must be a parent
// module of the current one. This means, for instance, that we don't have
// to worry about the possibility of the visibility paths "diverging".

#[derive(PartialEq)]
#[cfg_attr(test, derive(Debug))]
pub enum PathComponent {
  /// Effectively, this means private.
  Inherited,
  Pub,
  Crate,
  InSelf,
  InSuper,
  Mod(String)
}

#[cfg_attr(test, derive(PartialEq, Debug))]
pub enum FieldVisibility {
  /// The inner field is *at least* as visible as its containing struct.
  Visible,
  /// The inner field is less visible than its containing struct.
  Restricted,
  /// We can't figure out how the visibilities relate, probably due to the
  /// paths starting at different points (e.g. one is self and the other
  /// is ::a::b::c)
  CantDetermine
}

/// Check what the relation between the given struct's visibility and the
/// field's visibility is.
pub fn field_visibility(struct_vis: &syn::Visibility, field_vis: &syn::Visibility) -> FieldVisibility {
  let struct_vis = to_path(struct_vis);
  let field_vis = to_path(field_vis);

  fn check_head(struct_vis: &[PathComponent], field_vis: &[PathComponent]) -> FieldVisibility {
    match (struct_vis.split_first(), field_vis.split_first()) {
      (_, None)
        | (Some((&PathComponent::Inherited, _)), _)
        => FieldVisibility::Visible,
      (None, _)
        | (_, Some((&PathComponent::Inherited, _)))
        => FieldVisibility::Restricted,
      (Some((sh, sr)), Some((fh, fr))) => if sh == fh {
        check_head(sr, fr)
      } else {
        FieldVisibility::CantDetermine
      }
    }
  }

  // If the field is marked `pub`, then we know it's definitely visible...
  if &field_vis == &vec![ PathComponent::Pub ] {
    return FieldVisibility::Visible;
  }

  // ...and if that's not the case, but the struct is marked `pub`, we know
  // the field is definitely restricted.
  if &struct_vis == &vec![ PathComponent::Pub ] {
    return FieldVisibility::Restricted;
  }

  check_head(&struct_vis, &field_vis)
}

fn to_path(path: &syn::Visibility) -> Vec<PathComponent> {
  use syn::Visibility::*;

  // Visibilities of pub(self) and pub(in self) are equivalent to
  // no visibility modifier at all, so we need to make sure using
  // that doesn't report false positives.
  //
  // Technically this still doesn't catch *really* weird cases, like
  // if someone wrote pub(in self::self), but it should be good
  // enough.
  if is_pub_self(path) {
    return vec![ PathComponent::Inherited ];
  }

  match path {
    &Public(..) => vec![ PathComponent::Pub ],
    &Crate(..) => vec![ PathComponent::Pub, PathComponent::Crate ],
    &Inherited => vec![ PathComponent::Inherited ],
    &Restricted(ref vis) => to_path_restricted(&vis.path)
  }
}

fn to_path_restricted(path: &syn::Path) -> Vec<PathComponent> {
  let segments = path.segments.iter()
    .map(|path_segment| &path_segment.ident)
    .collect_vec();

  match segments.split_first() {
    None => vec![],
    Some((ident, rest)) => {
      let mut result;

      if *ident == "self" {
        result = vec![ PathComponent::InSelf ];
      } else if *ident == "super" {
        result = vec![ PathComponent::InSuper ];
      } else if *ident == "crate" {
        result = vec![ PathComponent::Pub, PathComponent::Crate ];
      } else {
        // We add these components in non-self/super paths to allow us to
        // match them up with visibilities like `pub` and `pub(crate)`.
        result = vec![ PathComponent::Pub, PathComponent::Crate, PathComponent::Mod(ident.to_string()) ];
      }

      let rest = rest.iter()
        .map(|ident| PathComponent::Mod(ident.to_string()));

      result.extend(rest);

      result
    }
  }
}

fn is_pub_self(vis: &syn::Visibility) -> bool {
  if let syn::Visibility::Restricted(syn::VisRestricted { ref path, .. }) = vis {
    path.is_ident("self")
  } else {
    false
  }
}

#[cfg(test)]
mod path_convert_tests {
  use std::convert::From;

  use syn::{self, Visibility};

  use super::{PathComponent, to_path};
  use super::PathComponent::*;

  impl<'a> From<&'a str> for PathComponent {
    fn from(input: &'a str) -> Self {
      Mod(input.to_string())
    }
  }

  macro_rules! vis_test {
    ($test_name:ident => $input:expr; $($component:expr),+) => {
      #[test]
      fn $test_name() {
        let vis: Visibility = syn::parse_str($input)
          .expect("path input is structured incorrectly!");
        let vis = to_path(&vis);

        let expected = vec![ $($component.into()),+ ];

        assert_eq!(&vis, &expected);
      }
    }
  }

  vis_test!(vis_test1 => "pub"; Pub);
  vis_test!(vis_test2 => "pub(crate)"; Pub, Crate);
  vis_test!(vis_test3 => ""; Inherited);
  vis_test!(vis_test4 => "pub(self)"; Inherited);
  vis_test!(vis_test5 => "pub(super)"; InSuper);
  vis_test!(vis_test6 => "pub(in ::a::b::c)"; Pub, Crate, "a", "b", "c");
  vis_test!(vis_test7 => "pub(in ::super::b)"; InSuper, "b");
}

#[cfg(test)]
mod field_visibility_tests {
  use syn::{self, Visibility};

  use super::field_visibility;
  use super::FieldVisibility::*;

  macro_rules! field_vis_test {
    ($test_name:ident => $struct_vis: expr; $field_vis: expr; $vis: expr) => {
      #[test]
      fn $test_name() {
        let struct_vis: Visibility = syn::parse_str($struct_vis)
          .expect("failed to parse struct visibility");
        let field_vis: Visibility = syn::parse_str($field_vis)
          .expect("failed to parse field visibility");

        let vis = field_visibility(&struct_vis, &field_vis);

        assert_eq!(vis, $vis);
      }
    }
  }

  field_vis_test!(test_field_vis1 => "pub"; "pub"; Visible);
  field_vis_test!(test_field_vis2 => ""; ""; Visible);
  field_vis_test!(test_field_vis3 => "pub(in a::b::c)"; "pub(in a::b)"; Visible);
  field_vis_test!(test_field_vis4 => "pub(in a::b)"; "pub(in a::b::c)"; Restricted);
  field_vis_test!(test_field_vis5 => "pub"; "pub(crate)"; Restricted);
  field_vis_test!(test_field_vis6 => "pub(crate)"; "pub(in a::b::c)"; Restricted);
  field_vis_test!(test_field_vis7 => "pub"; ""; Restricted);
  field_vis_test!(test_field_vis8 => ""; "pub"; Visible);
  field_vis_test!(test_field_vis9 => "pub(in a::b::c)"; "pub(self)"; Restricted);
  field_vis_test!(test_field_vis10 => "pub(in a::b::c)"; "pub(super)"; CantDetermine);
  field_vis_test!(test_field_vis11 => "pub"; "pub(self)"; Restricted);
  field_vis_test!(test_field_vis12 => "pub(in a::b::c)"; "pub"; Visible);
  field_vis_test!(test_field_vis13 => "pub(self)"; "pub(self)"; Visible);
  field_vis_test!(test_field_vis14 => "pub(super)"; "pub(super)"; Visible);
  field_vis_test!(test_field_vis15 => "pub(crate)"; "pub(crate)"; Visible);
  field_vis_test!(test_field_vis16 => "pub(in a::b::c)"; "pub(in a::b::c)"; Visible);
}
