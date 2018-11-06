# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.7.0] - 2018-11-06
### Changed
- overhaul of FDR machinery to make it only one tool that outputs a BCF filtered at a provided FDR threshold
- a number of performance optimizations leading to more than an order of magnitude speedup of ProSolo, a tool using this library, [with details described in its repo](https://github.com/ProSolo/prosolo/issues/2); most importantly these changes are:
  - [caching the CIGAR string of reads](https://github.com/PROSIC/libprosic/pull/34)
  - caching likelihood point estimates for reuse in different (two-dimensional) Events (see PRs [36](https://github.com/PROSIC/libprosic/pull/36) and [40](https://github.com/PROSIC/libprosic/pull/40))
  - [caching of `prob_rho()` function in the SingleCellBulkModel](https://github.com/PROSIC/libprosic/pull/46)
  - [minimize the number of `ln` operations performed, caching some values per program run, per pileup or per observation](https://github.com/PROSIC/libprosic/pull/48)
- [minor bugfix for `ContinuousAlleleFreqs` ranges depicting points](https://github.com/PROSIC/libprosic/pull/47)
- introduced `cargo fmt` use, including in continuous integration tests
- [insert size estimation from data and handling as `alignment_properties`, removing the necessity to provide them for each run via command-line arguments](https://github.com/PROSIC/libprosic/pull/41)
- [binarise MAPQ to 0 vs. maximum observed MAPQ at a site, to remove reference mapping bias](https://github.com/PROSIC/libprosic/pull/38) -- using the given instead of the binarised MAPQ could be selected via a command line option in the downstream tool that sets `use_mapq`
- overhaul of observation extraction mechanism and sampling probability calculation for indels, looking at full fragments (e.g. read pairs) jointly instead of looking at reads separately( see PRs [28](https://github.com/PROSIC/libprosic/pull/28), [29](https://github.com/PROSIC/libprosic/pull/29) and [30](https://github.com/PROSIC/libprosic/pull/30))
- [added functionality to calculate likelihoods and calls for candidate sites without a candidate alternative allele, i.e. homozygous reference candidates](https://github.com/PROSIC/libprosic/pull/24)
- dependency updates (`rust-htslib` to `0.22`, `rust-bio` to `0.23`)

## [0.6.0] - 2018-01-12
### Changed

- added functionality to filter by FDR thresholds based on sets of `Events` and improved code for summing up Event likelihoods (incl. fixing numerical overshoot issues)
- dependency updates (`rust-htslib` to `0.16`, `rust-bio` to `0.16` from `rustc-serialize` to `csv` and `serde`)
- some test fixes


## [0.5.0] - 2017-11-17
### Changed
- FDR control now works for sets of events.
- The pair HMM implementation has been finetuned and improved.
- Artifacts from weird mapper decisions are now better detected.


## [0.4.0] - 2017-07-06
### Changed
- Refactored and fixed the false discovery rate estimation code.
- Improved indel likelihood calculation (pairHMM based indel likelihoods; do not normalize indel likelihoods for ref and alt, only consider a small window around the indel)
- New SingleCellBulkModel for calling variants in single cell data against a bulk background sample.
- Update rust-bio dependency to 0.14.* or newer.
- Marginalization only over specified Events, not the full allele frequency ranges.
- Full-blown Cigar string parsing with CigarErrors and full single nucleotide variant (SNV) support, externalized via rust-htslib dependency.
- New Ranges with inclusive/exclusive end points.
- Added commitizen.


## [0.3.0] - 2017-05-04
### Changed
- Ensure fair sampling of reads supporting ALT and REF allele.


## [0.2.0] - 2017-03-28
### Changed
- calculate ALT and REF likelihood by realigning against pseudo-haplotypes

## [0.1.0] - 2016-11-02
### Added
- initial release
