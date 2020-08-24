# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [2.2.0] - 2020-08-24
- Allow parallelization via setting the number of threads.
- BCF output is now unsorted, and must be sorted afterwards with bcftools.
- Imprecise variants are skipped for now, until proper support is implemented.

## [2.1.0] - 2020-08-11
- Infer missing antisense breakends (sometimes callers only report one direction although the other is necessary as well).
- Support for single breakends.
- Support for arbitrary replacements.
- Fixed a bug with posterior odds filtration leading to underestimation of the odds.

## [2.0.1] - 2020-08-06
- Fixed allele frequency biases that occurred due to missed evidence when investigating breakends.
- Fixed pattern-too-long error that occurred in some corner cases.
- Fixed error occurring when investigation breakends without event tag. These are skipped for now. Special handling will be added later.

## [2.0.0] - 2020-07-09
- Add support for inversions, deletions and breakends.
- Slightly modified CLI options.

## [1.7.3] - 2020-06-17
- Polished TMB plots.
- Allow to omit insert size evidence. This is important for amplicon data, where indels do not impact insert size.

## [1.7.2] - 2020-06-03
- Fixed possible values in CLI.

## [1.7.1] - 2020-06-02
- Fixed compiler issue when building varlociraptor with certain rust versions.
- Use latest rust-bio release.

## [1.7.0] - 2020-05-27
- Improved TMB estimation plots, now offering three modes (hist, curve, and stratified), as well as being able to better see multiple scales.

## [1.6.4] - 2020-03-18
- Report negative SVLEN for deletions again. The htslib bug is actually fixed already. However, it is mandatory to rerun varlociraptor preprocess to avoid it downstream.

## [1.6.3] - 2020-03-17
- Cleanup of debug messages.

## [1.6.2] - 2020-03-17
- Fix odds filtering arguments.
- Improved error messages.
- Improved insert size estimation.
- Use htslib 1.10.
- Work around htslib bug that misinterprets negative SVLENs: SVLEN is now always positive.

## [1.5.0] - 2019-12-04
- Introduce SNV selectors in the event formula language.
- Work around a segmentation fault in Htslib when estimating tmb.

## [1.4.4] - 2019-12-03
- Modify observation format such that encoding cannot lead to corner cases where BCF interprets a value as a vector end marker.
- Fix test case generation

## [1.4.3] - 2019-12-02
- Make observation reading and writing more robust by controlling the length of INFO field vectors. This fixes a bug when reading leads to invalid strand information.

## [1.4.2] - 2019-12-01
- Avoid unsafe memory operations when writing and reading preprocessed observations.

## [1.4.1] - 2019-11-29
- Work around a segmentation fault caused by Htslib.
- Update command line usage instructions.

## [1.4.0] - 2019-11-28
- Separate calling and preprocessing of observations. This allows to easily reuse large parts of the computation when changing the scenario. Further, it allows to parallelize across samples.

## [1.3.0] - 2019-11-13
- Enable setting of the strandedness of the sequencing protocol at the command line (same or opposite).
- Allow contig-specific definition of the allele frequency universe in the grammar.

## [1.2.2] - 2019-10-04
- Explicitly annotate unit (PHRED, or linear) of probabilities in BCF output.

## [1.2.1] - 2019-09-27
- Fixed a bug in usage of htslib that could lead to a segmentation fault when estimating TMB.

## [1.2.0] - 2019-09-27
- Estimation of tumor mutational burden as a curve of expectations over minimum allele frequencies.
- Various small bug fixes.

## [1.1.0] - 2019-07-03
- Numerical improvements to strand bias model.

## [1.0.1] - 2019-05-15
- Fix numerical issue with indel variant calling.
- Fix indel window length issues with reads >100bp.

## [1.0.0] - 2019-05-09
- Rewrite of statistical model, using rust-bios new trait machinery.
- Generalization of the model, allowing arbitrary allele frequency events and sample numbers.
- Inclusion of strand bias in the model.
- Better VCF output.
- Various bug fixes.
- New test framework.
- CLI overhaul.
- Experimental tumor-normal CNV calling.
- Various improvements to calling accuracy and fixed corner cases.

## [0.7.0] - 2018-11-06
### Changed
- overhaul of FDR machinery to make it only one tool that outputs a BCF filtered at a provided FDR threshold
- a number of performance optimizations leading to more than an order of magnitude speedup of ProSolo, a tool using this library, [with details described in its repo](https://github.com/ProSolo/prosolo/issues/2); most importantly these changes are:
  - [caching the CIGAR string of reads](https://github.com/varlociraptor/varlociraptor/pull/34)
  - caching likelihood point estimates for reuse in different (two-dimensional) Events (see PRs [36](https://github.com/varlociraptor/varlociraptor/pull/36) and [40](https://github.com/varlociraptor/varlociraptor/pull/40))
  - [caching of `prob_rho()` function in the SingleCellBulkModel](https://github.com/varlociraptor/varlociraptor/pull/46)
  - [minimize the number of `ln` operations performed, caching some values per program run, per pileup or per observation](https://github.com/varlociraptor/varlociraptor/pull/48)
- remove use of flamegraphing crates in favor of [perf flamegraphing](https://gist.github.com/dlaehnemann/df31787c41bd50c0fe223df07cf6eb89)
- [minor bugfix for `ContinuousAlleleFreqs` ranges depicting points](https://github.com/varlociraptor/varlociraptor/pull/47)
- introduced `cargo fmt` use, including in continuous integration tests
- [insert size estimation from data and handling as `alignment_properties`, removing the necessity to provide them for each run via command-line arguments](https://github.com/varlociraptor/varlociraptor/pull/41)
- [binarise MAPQ to 0 vs. maximum observed MAPQ at a site, to remove reference mapping bias](https://github.com/varlociraptor/varlociraptor/pull/38) -- using the given instead of the binarised MAPQ could be selected via a command line option in the downstream tool that sets `use_mapq`
- overhaul of observation extraction mechanism and sampling probability calculation for indels, looking at full fragments (e.g. read pairs) jointly instead of looking at reads separately( see PRs [28](https://github.com/varlociraptor/varlociraptor/pull/28), [29](https://github.com/varlociraptor/varlociraptor/pull/29) and [30](https://github.com/varlociraptor/varlociraptor/pull/30))
- [added functionality to calculate likelihoods and calls for candidate sites without a candidate alternative allele, i.e. homozygous reference candidates](https://github.com/varlociraptor/varlociraptor/pull/24)
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
