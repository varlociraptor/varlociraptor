# Varlociraptor

[![Travis](https://img.shields.io/travis/varlociraptor/varlociraptor/master.svg?maxAge=2592000)](https://travis-ci.org/varlociraptor/varlociraptor)
[![Codecov](https://img.shields.io/codecov/c/github/varlociraptor/varlociraptor/master.svg)](https://codecov.io/gh/varlociraptor/varlociraptor)
[![API docs](https://img.shields.io/badge/API-documentation-blue.svg)](https://docs.rs/varlociraptor)
[![Crates.io](https://img.shields.io/crates/d/varlociraptor.svg)](https://crates.io/crates/varlociraptor)
[![Commitizen friendly](https://img.shields.io/badge/commitizen-friendly-brightgreen.svg)](http://commitizen.github.io/cz-cli/)

Varlociraptor implements a novel, unified fully uncertainty-aware approach to genomic variant calling in arbitrary scenarios. 

### Key features

* Calls SNVs and indels in all length ranges (from small to structural) with a unified statistical model.
* The statistical model entails all possible sources of uncertainty.
* Resulting variant calls can be filtered by false discovery rate. No parameter tuning necessary.
* Unbiased, maximum a posteriori allele frequency estimates are provided with each call.

### Calling modes

* Tumor-normal-calling, classifying variants as somatic in tumor, somatic in normal, germline, or absent.
* Generic, grammar based calling, allowing to classify arbitrary scenarios.

**For details, see the homepage: https://varlociraptor.github.io**
