# Changelog

## [2.2.0](https://github.com/rust-bio/hts-sys/compare/v2.1.4...v2.2.0) (2024-12-02)


### Features

* add hfile.h to wrapper.h so hopen appears in bindings ([#17](https://github.com/rust-bio/hts-sys/issues/17)) ([e9eaef7](https://github.com/rust-bio/hts-sys/commit/e9eaef7cba36d7a18537ed64d4656f33a7fea018))

## [2.1.4](https://github.com/rust-bio/hts-sys/compare/v2.1.3...v2.1.4) (2024-05-10)


### Bug Fixes

* add missing sam_mods to bindings ([#13](https://github.com/rust-bio/hts-sys/issues/13)) ([a4fc270](https://github.com/rust-bio/hts-sys/commit/a4fc270cc525a596bc54b1388e26db2f8cb9fe4c))

## [2.1.3](https://github.com/rust-bio/hts-sys/compare/v2.1.2...v2.1.3) (2024-05-06)


### Bug Fixes

* build on aarch64-apple-darwin ([#10](https://github.com/rust-bio/hts-sys/issues/10)) ([534c463](https://github.com/rust-bio/hts-sys/commit/534c4638a375d7f45d0de57564f6f42443d1ef97))

## [2.1.2](https://github.com/rust-bio/hts-sys/compare/v2.1.1...v2.1.2) (2024-03-15)


### Miscellaneous Chores

* release 2.1.2 ([32f185f](https://github.com/rust-bio/hts-sys/commit/32f185f5f963c8cd1574dfbb1858775be4dfc6eb))

## [2.1.1](https://github.com/rust-bio/hts-sys/compare/v2.1.0...v2.1.1) (2023-06-14)


### Bug Fixes

* remove the `hts-sys/` subpath from `--manifest-path` specification ([#5](https://github.com/rust-bio/hts-sys/issues/5)) ([53c2e23](https://github.com/rust-bio/hts-sys/commit/53c2e23c72da9dfbb68ebaa1717f4a1636709cb8))

## [2.1.0](https://github.com/rust-bio/hts-sys/compare/v2.0.3...v2.1.0) (2023-06-13)


### Features

* bindings for basemod API ([6e6c585](https://github.com/rust-bio/hts-sys/commit/6e6c5852c9281195429290f771dd44183679921b))

## [2.0.3](https://github.com/rust-bio/rust-htslib/compare/hts-sys-v2.0.2...hts-sys-v2.0.3) (2022-10-13)


### Performance Improvements

* update htslib and corresponding bindings to 1.16 ([#366](https://github.com/rust-bio/rust-htslib/issues/366)) ([f597ce0](https://github.com/rust-bio/rust-htslib/commit/f597ce0451e3f3c393166a7291486bbc2bde4c39))

### [2.0.2](https://www.github.com/rust-bio/rust-htslib/compare/hts-sys-v2.0.1...hts-sys-v2.0.2) (2021-07-17)


### Bug Fixes

* Configuration when cross-compiling. Even when cross-compiling, build.rs runs on the build host. Hence within build.rs `#[cfg(target_os)]` always reflects the host, not the target. Use $CARGO_CFG_TARGET_OS instead to query target properties. ([#329](https://www.github.com/rust-bio/rust-htslib/issues/329)) ([d5198e6](https://www.github.com/rust-bio/rust-htslib/commit/d5198e6c777fdbbfdd9c73a820f1be983a458ce2))

### [2.0.1](https://www.github.com/rust-bio/rust-htslib/compare/hts-sys-v2.0.0...hts-sys-v2.0.1) (2021-07-06)


### Bug Fixes

* dummy release ([74d1565](https://www.github.com/rust-bio/rust-htslib/commit/74d1565329fc862f1172c0925c7b66ceb8bcf988))
* dummy release ([af2f84e](https://www.github.com/rust-bio/rust-htslib/commit/af2f84eb0411507f8866b3cc05e9a6ba9d81d172))

## [2.0.0](https://www.github.com/rust-bio/rust-htslib/compare/hts-sys-v1.0.0...hts-sys-v2.0.0) (2021-07-06)


### ⚠ BREAKING CHANGES

* dummy major version bump to move away from previous versions that were following htslib versions.

### Bug Fixes

* dummy major version bump to move away from previous versions that were following htslib versions. ([aaa70a8](https://www.github.com/rust-bio/rust-htslib/commit/aaa70a85ef9a908d3b101f23879189e84a15d23f))

## [1.0.0](https://www.github.com/rust-bio/rust-htslib/compare/hts-sys-v0.1.0...hts-sys-v1.0.0) (2021-07-06)


### ⚠ BREAKING CHANGES

* bump to new major version (for technical reasons).
* dummy breaking change to increase hts-sys major version.

### Bug Fixes

* bump to new major version (for technical reasons). ([9c6db30](https://www.github.com/rust-bio/rust-htslib/commit/9c6db3060818692070db1411d63e113dc7effd64))
* dummy breaking change to increase hts-sys major version. ([93415cb](https://www.github.com/rust-bio/rust-htslib/commit/93415cbb82e4f11d257a2b2cedba2664f86a034d))
* dummy changes ([3af5ede](https://www.github.com/rust-bio/rust-htslib/commit/3af5ede13a6b44ce5d1e7f0eb90836a692e711ec))
* update changelog ([deef08f](https://www.github.com/rust-bio/rust-htslib/commit/deef08feb0b5ba2d8abf98f2cc6d327236da8aef))

## [2.0.0](https://www.github.com/rust-bio/rust-htslib/compare/hts-sys-v1.11.1-fix1...hts-sys-v2.0.0) -  (2021-07-06)


### Bug Fixes

* dummy release ([b97915f](https://www.github.com/rust-bio/rust-htslib/commit/b97915f2c70da4c914f2e69861bf78eec5979baf))
* trigger dummy release ([7c5a7de](https://www.github.com/rust-bio/rust-htslib/commit/7c5a7de33e2a92052126e5f44389d421974d1e02))


## 0.1.0 - 2021-07-06

### Bug Fixes

* dummy release
