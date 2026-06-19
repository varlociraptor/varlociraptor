# Changelog

## [Unreleased]

## [v0.3.0] -- 2019-12-03

* Removed the `map_ref` and `map_mut` functions because they were
  essentially useless.
* Renamed `map` to `siphon` to better reflect how it consumes the
  value.
* Added a new function `transform` in place of `map_mut`, which
  is actually useful because it allows chained updates.
* Added a new attribute `#[shrinkwrap(transformers)]`; `siphon`
  and `transform` are not derived by default, and must be enabled
  using this flag.

## [v0.2.3] -- 2019-11-18

* Improved error messages on visibility errors
* Attached error messages to correct source code spans instead of
  the derive attribute
* Fixed visibility checks around `pub(self)` and `pub(in self)` to
  be equivalent to private visibility

## [v0.2.2] -- 2019-11-17

* Updated dependencies to latest major versions

## [v0.2.1] -- 2019-01-24

* Added the ability to generate `#[nostd]`-compatible code through
  feature flags -- thanks @dazabani!

  * The default is still to emit code that uses `std`.
  * Use the feature flag `std` to emit `std`-compatible code, or omit
    it to emit code that uses `core`.

## [v0.2.0] -- 2018-02-10

* Added visibility checking on mutable derives to help prevent deriving
  mutable traits when inner field is less mutable than surrounding struct
* Added `#[shrinkwrap(unsafe_ignore_visibility)]` flag to override this
  behavior when desired.
* Removed `#[derive(ShrinkwrapMut)]` proc macro; replaced with
  `#[shrinkwrap(mutable)]` attribute.

## [v0.1.1] -- 2018-02-07

* Added a changelog
* Implemented mapping methods `map()`, `map_ref()`, `map_mut()` for
  mapping functions over wrapped values (useful for function call chaining)
* Added support for structs with lifetimes and generic parameters

## [v0.1.0] -- 2018-02-06

* Split out derivation of mutable traits (`DerefMut`, `BorrowMut`, `AsMut`) into
  separate derive trait `ShrinkwrapMut`

## [v0.0.2] -- 2018-02-04

* Fixed typoes in documentation -- no functionality changes

## [v0.0.1] -- 2018-02-04

* Initial release -- implemented `#[derive(Shrinkwrap)]` to auto-derive
  `Deref`, `DerefMut`, `Borrow`, `BorrowMut`, `AsRef`, and `AsMut`.
