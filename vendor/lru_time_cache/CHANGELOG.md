# Changelog

All notable changes to this project will be documented in this file. See [standard-version](https://github.com/conventional-changelog/standard-version) for commit guidelines.

### [0.11.11](https://github.com/maidsafe/lru_time_cache/compare/v0.11.10...v0.11.11) (2021-06-09)

### [0.11.10](https://github.com/maidsafe/lru_time_cache/compare/v0.11.9...v0.11.10) (2021-03-03)

### [0.11.9](https://github.com/maidsafe/lru_time_cache/compare/v0.11.8...v0.11.9) (2021-03-01)

### [0.11.8](https://github.com/maidsafe/lru_time_cache/compare/v0.11.7...v0.11.8) (2021-02-24)

### [0.11.7](https://github.com/maidsafe/lru_time_cache/compare/v0.11.6...v0.11.7) (2021-02-10)

### [0.11.6](https://github.com/maidsafe/lru_time_cache/compare/v0.11.5...v0.11.6) (2021-02-03)

### [0.11.5](https://github.com/maidsafe/lru_time_cache/compare/v0.11.4...v0.11.5) (2021-01-20)

### [0.11.4](https://github.com/maidsafe/lru_time_cache/compare/v0.11.3...v0.11.4) (2021-01-18)

### [0.11.3](https://github.com/maidsafe/lru_time_cache/compare/v0.11.2...v0.11.3) (2020-11-23)

### [0.11.2](https://github.com/maidsafe/lru_time_cache/compare/v0.11.1...v0.11.2) (2020-10-09)

### [0.11.1](https://github.com/maidsafe/lru_time_cache/compare/v0.11.0...v0.11.1) (2020-09-17)

### [0.11.0](https://github.com/maidsafe/lru_time_cache/compare/v0.10.0...v0.11.0) (2020-09-01)

* update to reference renamed sn_fake_clock crate

### [0.10.0](https://github.com/maidsafe/lru_time_cache/compare/0.9.0...v0.10.0) (2020-03-20)

* Move iterators into a separate module
* Split different test cases
* Update `LruCache::peek_iter()` order - most recently used items will be
  produced first.
* Fix edge cases related to time atomicity
* Fix atomicity of insert on entry
* Make library work in Rust stable 1.41.
* Use `next` to get the first element in the cache

### [0.9.0](https://github.com/maidsafe/lru_time_cache/compare/0.8.1...0.9.0) (2019-02-20)

* API to get expired or pushed out items from the LRU
* Update `LruCache::iter()` order - most recently used items will be produced
  first.
* Update `rand` dependency

### [0.8.1](https://github.com/maidsafe/lru_time_cache/compare/0.8.0...0.8.1) (2019-01-05)
* Update to dual license (MIT/BSD)

### [0.8.0](https://github.com/maidsafe/lru_time_cache/compare/0.7.0...0.8.0) (2018-01-05)
* Use rust 1.22.1 stable / 2017-12-02 nightly
* rustfmt 0.9.0 and clippy-0.0.175

### [0.7.0](https://github.com/maidsafe/lru_time_cache/compare/0.6.0...0.7.0) (2017-07-25)
* Use rust 1.19 stable / 2017-07-20 nightly
* rustfmt 0.9.0 and clippy-0.0.144
* Replace -Zno-trans with cargo check
* Make appveyor script using fixed version of stable
* Use cargo_install from QA

### [0.6.0](https://github.com/maidsafe/lru_time_cache/compare/0.5.0...0.6.0) (2017-04-12)
* Add support for using fake clock.
* CI, README, rustfmt and clippy cleanups.

### [0.5.0](https://github.com/maidsafe/lru_time_cache/compare/0.4.0...0.5.0) (2016-08-03)
* Add `iter` and remove obsolete `retrieve_all` methods.

### [0.4.0](https://github.com/maidsafe/lru_time_cache/compare/0.3.1...0.4.0) (2020-09-17)
* Add `clear`, `peek` and `peek_iter` methods.

### [0.3.1](https://github.com/maidsafe/lru_time_cache/compare/0.3.0...0.3.1) (2016-04-26)
* Fix arithmetic operation overflows.

### [0.3.0](https://github.com/maidsafe/lru_time_cache/compare/0.2.7...0.3.0) (2016-04-20)
* Remove dependency on the time crate.
* Use std::time::Duration in the API

### [0.2.7](https://github.com/maidsafe/lru_time_cache/compare/0.2.6...0.2.7) (2016-03-04)
* Updated dependencies.

### [0.2.6](https://github.com/maidsafe/lru_time_cache/compare/0.2.5...0.2.6) (2016-01-21)
* Allow non-Clone Value types.

### [0.2.5](https://github.com/maidsafe/lru_time_cache/compare/0.2.4...0.2.5) (2015-12-11)
* Update time to live when accessing elements.

### [0.2.4](https://github.com/maidsafe/lru_time_cache/compare/0.2.3...0.2.4) (2015-11-13)
* Update deprecated item, replaced by `std::thread::sleep`.

### [0.2.3](https://github.com/maidsafe/lru_time_cache/compare/0.2.2...0.2.3) (2015-11-13)
* Remove wildcard dependencies.

### [0.2.2] (2015-09-14)
* Removes expired values before accessing elements. Removed deprecated check method.

### [0.2.1] (2015-09-12)
* Provides a getter to fetch all key value pairs in order.
* Removed `add` function (deprecated in favor of the `insert` function from v0.1.6)

### [0.1.7 - 0.2.0] (2015-07-06)
* [#21] (https://github.com/maidsafe/lru_time_cache/issues/21) Enforced lint checks
* Rename `check` to `contains_key`

### [0.1.6] API additions (2015-05-25)
* Implement the `entry` function
* Implement the `insert` function as a replacement for `add` (with same semantics as Rust's standard `Map::insert` functions)
* Implement the `get_mut`

### [0.0.0 - 0.1.5] First implementation (2015-05-02)
* Implement add_key_value
* Test add_key_value (time and size based tests)
* Implement check
* Test check (time and size based tests)
* Implement get(key)
* Test get (time and size based tests)
* API version 0.8.0
* Implement delete_key
* Test delete_key (time and size based tests)
* API version 0.1.0
