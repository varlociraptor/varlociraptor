# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [4.8.0](https://www.github.com/varlociraptor/varlociraptor/compare/v4.7.0...v4.8.0) (2021-12-07)


### Features

* accurate consideration of multiallelic variants ([#239](https://www.github.com/varlociraptor/varlociraptor/issues/239)) ([d8e646e](https://www.github.com/varlociraptor/varlociraptor/commit/d8e646e567317bf02bf665c173ec6d58d3f06a0e))
* ability to display each processed record in the logging output of `varlociraptor call` via --log-mode each-record.

## [4.7.0](https://www.github.com/varlociraptor/varlociraptor/compare/v4.6.0...v4.7.0) (2021-11-17)


### Features

* do not anymore require genome size definition in scenarios ([#231](https://www.github.com/varlociraptor/varlociraptor/issues/231)) ([a707d02](https://www.github.com/varlociraptor/varlociraptor/commit/a707d02fb90e5f75182579f6c2bfc5ac06820e56))
* output sampled VAF densities (in the AFD format tag) ([#234](https://www.github.com/varlociraptor/varlociraptor/issues/234)) ([35ced8f](https://www.github.com/varlociraptor/varlociraptor/commit/35ced8f3066bc00367aafc8d2e76b8a0630d061d))


### Bug Fixes

* always keep replacements as they are, even if they are just simple insertions or deletions ([4a8e6ed](https://www.github.com/varlociraptor/varlociraptor/commit/4a8e6edda70e2e132853fcdbea100cc4d621428c))
* fixed determination of homopolymer indels from replacements ([#233](https://www.github.com/varlociraptor/varlociraptor/issues/233)) ([4318e5f](https://www.github.com/varlociraptor/varlociraptor/commit/4318e5fc2ce3231264772a708a365aa854aea010))
* fixed update of possible events when switching between variant types. ([f9037f7](https://www.github.com/varlociraptor/varlociraptor/commit/f9037f764ad9c8af9bc5fb2c9812ee50bcff30e6))
* simplified but more accurate homopolymer error model ([#230](https://www.github.com/varlociraptor/varlociraptor/issues/230)) ([7f8cbd8](https://www.github.com/varlociraptor/varlociraptor/commit/7f8cbd885b5d397799d23fcc7362be7afbda4342))
* understandable error message if BCF record contains invalid positions ([#228](https://www.github.com/varlociraptor/varlociraptor/issues/228)) ([6533a70](https://www.github.com/varlociraptor/varlociraptor/commit/6533a70b881282adc8a7e857f91b506e6dc548e1))

## [4.6.0](https://www.github.com/varlociraptor/varlociraptor/compare/v4.5.0...v4.6.0) (2021-11-02)


### Features

* better homopolymer error detection ([#226](https://www.github.com/varlociraptor/varlociraptor/issues/226)) ([b95d03d](https://www.github.com/varlociraptor/varlociraptor/commit/b95d03d99519cb00c0a59aa61ab834a56e5f3777))
* eliminate false positives and false negatives: MAPQ adjustment (adjust MAPQ by mean pileup MAPQ; this gets rid of stochastically inflated MAPQs), do not estimate strand bias if all reads come from the same strand, infer wildtype forward strand rate from non-alt reads if possible ([#222](https://www.github.com/varlociraptor/varlociraptor/issues/222)) ([d54ccb0](https://www.github.com/varlociraptor/varlociraptor/commit/d54ccb07b02788d026ee1defe540a84587d63998))


### Bug Fixes

* properly handle reads with both allele probs -inf ([6b04fd6](https://www.github.com/varlociraptor/varlociraptor/commit/6b04fd6a118524eb27ba469309c620771e4a7760))

## [4.5.0](https://www.github.com/varlociraptor/varlociraptor/compare/v4.4.3...v4.5.0) (2021-10-07)


### Features

* added inference of classical genotypes from Varlociraptor's AF field (see varlociraptor genotype --help) ([b12c0f5](https://www.github.com/varlociraptor/varlociraptor/commit/b12c0f5eec405dd48089b98c202fb16ad25d8d7c))


### Bug Fixes

* Add check for order of records in bcf (before, coordinate sorting was implicitly assumed but not checked, leading to overflow errors upon violation) ([#218](https://www.github.com/varlociraptor/varlociraptor/issues/218)) ([c3d20a9](https://www.github.com/varlociraptor/varlociraptor/commit/c3d20a95b0616901978cdec81c1d80b8e87786b6))

### [4.4.3](https://www.github.com/varlociraptor/varlociraptor/compare/v4.4.2...v4.4.3) (2021-09-30)


### Bug Fixes

* fix out of bounds error when evaluating replacement records ([#216](https://www.github.com/varlociraptor/varlociraptor/issues/216)) ([7d5fde1](https://www.github.com/varlociraptor/varlociraptor/commit/7d5fde1e47b8ff7e1a514d16a830890112976449))

### [4.4.2](https://www.github.com/varlociraptor/varlociraptor/compare/v4.4.1...v4.4.2) (2021-09-14)


### Performance Improvements

* do not consider empty pileups for bias preprocessing ([#214](https://www.github.com/varlociraptor/varlociraptor/issues/214)) ([e9beec3](https://www.github.com/varlociraptor/varlociraptor/commit/e9beec3db1ed7cceb6a302f2fc1ecc0e7e612ad0))

### [4.4.1](https://www.github.com/varlociraptor/varlociraptor/compare/v4.4.0...v4.4.1) (2021-09-10)


### Performance Improvements

* various speed improvements for scenario evaluations ([#212](https://www.github.com/varlociraptor/varlociraptor/issues/212)) ([f865e1b](https://www.github.com/varlociraptor/varlociraptor/commit/f865e1bc5ae7e4369d4c9a4ece104fc6f33d6447))

## [4.4.0](https://www.github.com/varlociraptor/varlociraptor/compare/v4.3.0...v4.4.0) (2021-09-06)


### Features

* overhauled somatic prior ([#209](https://www.github.com/varlociraptor/varlociraptor/issues/209)) ([7530469](https://www.github.com/varlociraptor/varlociraptor/commit/7530469d2eb838baac897f27be2b7755c393ef54))

## [4.3.0](https://www.github.com/varlociraptor/varlociraptor/compare/v4.2.0...v4.3.0) (2021-09-02)


### Features

* adaptive integration (speeding up the evaluation of complex scenarios) ([#199](https://www.github.com/varlociraptor/varlociraptor/issues/199)) ([8043f94](https://www.github.com/varlociraptor/varlociraptor/commit/8043f94702402c20c69ede78fb3bed981fb9915b))


### Bug Fixes

* fixed record count in preprocess progress logging; fixed alt allele allocation error with certain longer replacements ([#207](https://www.github.com/varlociraptor/varlociraptor/issues/207)) ([ba66c36](https://www.github.com/varlociraptor/varlociraptor/commit/ba66c363b4d549c92227b85b30080b01c9a9b7b4))

## [4.2.0](https://www.github.com/varlociraptor/varlociraptor/compare/v4.1.3...v4.2.0) (2021-08-30)


### Features

* added VAF log fold change and comparison operators to the scenario event grammar ([#202](https://www.github.com/varlociraptor/varlociraptor/issues/202)) ([0f6da9e](https://www.github.com/varlociraptor/varlociraptor/commit/0f6da9e0b7cd337b39035ed3d7f5803a3ec51416))

### [4.1.3](https://www.github.com/varlociraptor/varlociraptor/compare/v4.1.2...v4.1.3) (2021-08-13)


### Bug Fixes

* erroneous integral boundary adjustment when having small VAF intervals and only very few observations. ([97f0124](https://www.github.com/varlociraptor/varlociraptor/commit/97f01244d97d24cd1eb22f89f1317da33b278a23))

### [4.1.2](https://www.github.com/varlociraptor/varlociraptor/compare/v4.1.1...v4.1.2) (2021-08-10)


### Bug Fixes

* improved divindel bias estimation precision (getting rid of false negatives due to erroneous divindel bias estimates) ([#196](https://www.github.com/varlociraptor/varlociraptor/issues/196)) ([6a0ac5d](https://www.github.com/varlociraptor/varlociraptor/commit/6a0ac5d2716c309cc88e1906412253e26bd4e5a2))

### [4.1.1](https://www.github.com/varlociraptor/varlociraptor/compare/v4.1.0...v4.1.1) (2021-07-22)


### Bug Fixes

* display 3 digits in normalized formulae ([f5e5261](https://www.github.com/varlociraptor/varlociraptor/commit/f5e5261b82499beb14ec6f5df7f3d8d63992643e))
* Removed check for coding variants ([#193](https://www.github.com/varlociraptor/varlociraptor/issues/193)) ([2f0f653](https://www.github.com/varlociraptor/varlociraptor/commit/2f0f653cc07479384f62e4f28ab75b347a40df6a))

## [4.1.0](https://www.github.com/varlociraptor/varlociraptor/compare/v4.0.1...v4.1.0) (2021-07-09)


### Features

* add divindel bias for detecting artifacts caused by diverging indels associated with the alt allele ([#189](https://www.github.com/varlociraptor/varlociraptor/issues/189)) ([5ab2620](https://www.github.com/varlociraptor/varlociraptor/commit/5ab2620fe2e0702a4ff57f54b7c5ce0251739e04))

### [4.0.1](https://www.github.com/varlociraptor/varlociraptor/compare/v4.0.0...v4.0.1) (2021-07-07)


### Bug Fixes

* testcase 71 sometimes fails because formula normalization via BDD is not deterministic ([#187](https://www.github.com/varlociraptor/varlociraptor/issues/187)) ([bdeddcc](https://www.github.com/varlociraptor/varlociraptor/commit/bdeddcc7594e49ebea4812e13303686fbec6ef69))

## [4.0.0](https://www.github.com/varlociraptor/varlociraptor/compare/v3.5.0...v4.0.0) (2021-07-07)


### ⚠ BREAKING CHANGES

* Rename estimate tmb subcommand into estimate mutational-burden (#182)

### Features

* release automation ([b6cc00f](https://www.github.com/varlociraptor/varlociraptor/commit/b6cc00ff5cc77c2760e12891540f36c326ca0c1b))
* improved formula simplification, detection of overlapping events ([#177](https://www.github.com/varlociraptor/varlociraptor/issues/177)) ([7e875a5](https://www.github.com/varlociraptor/varlociraptor/commit/7e875a5ff4a4afee9b690dd5ea05d6bee762ac47))

### Code Refactoring

* Rename estimate tmb subcommand into estimate mutational-burden ([#182](https://www.github.com/varlociraptor/varlociraptor/issues/182)) ([9fa5a1f](https://www.github.com/varlociraptor/varlociraptor/commit/9fa5a1f6cb58d37c9fdf5c235ae035cf096b10a0))

## [3.5.0] - 2021-07-02
- Improvements to SNV and MNV calling: in case of reads with indels, perform a realignment to properly assess evidences (@johanneskoester).
- Improved error messages (@johanneskoester).
- Add subcommand for scatterplotting VAFs between multiple samples (@jafors).

## [3.4.0] - 2021-05-25
- Allow to control local instead of global FDR (`varlociraptor filter-calls control-fdr --local`).
- Allow to configure mutation rate reduction factors for MNVs, Indels, and SVs.


## [3.3.2] - 2021-05-22
- Stop decode-phred subcommand from removing quoting from description field in header lines of decoded fields.

## [3.3.1] - 2021-05-22
- Fix debug output accidentally printing to STDOUT, thereby invalidating BCF.

## [3.3.0] - 2021-05-20
- Add softclip bias for detecting SNV or MNV artifacts induced by alignment issues (@johanneskoester).
- Generalize mendelian inheritance towards arbitrary ploidies (@johanneskoester).
- Add option for anonymizing testcases (`--testcase-anonymous`) (@johanneskoester)

## [3.2.0] - 2021-05-18
- Added SOBS field, showing simplified observations that omit all information but the alt allele evidence (@johanneskoester).
- Enable to use any event as an expression via $eventname in other events. This can be used to e.g. formulate "otherwise" events, via defining a negated disjuction of all other events (e.g. `otherwise: !($denovo | $absent)`).
- Add BDD based automatic simplification of all event formulas (@johanneskoester).
- Automatically omit unlikely biases: for a bias to be evaluated, it needs to be supported by at least two third of the observations which are strongly (by their Kass-Raftery-Score) supporting the alt allele. This is safe, because bias likelihoods will otherwise be astronomically small. Since the three biases (strand, read position, read orientation) add 2 + 1 + 2 = 5 additional formula tree evaluations, this should make `varlociraptor call` 6 times faster on loci without a bias and three times faster on loci with a bias (@johanneskoester).
- Several fixes for corner cases in the implementation of the prior (@johanneskoester).

## [3.1.0] - 2021-05-04
- Generalized TMB plots to mutational burden plots (not tumor specific). Added a multi-sample variant using barplots (@jafors).
- Added testcase for and handling of missing CIGAR operation lenghts observations (@nh13, @dlaehnemann, @johanneskoester).
- Use hash map instead of BTree for managing event space. This should result in some speed-up (@johanneskoester).
- Fixed bugs with sex chromosome handling and prior distribution (@christopher-schroeder, @dawidkrzeciesa, @johanneskoester).

## [3.0.1] - 2021-04-22
- Fixed a bug in the new prior implementation that led to almost infinite filling of a result cache, leading to a memory leak.

## [3.0.0] - 2021-04-12
- Added a prior distribution that is fully configurable via the variant calling grammar, covering population genetics (heterozygosity), mendelian inheritance, and tumor evolution (Williams et al.).

## [2.6.5] - 2021-03-19
- Fix posterior probabilities when having zero observations and no bias estimation

## [2.6.4] - 2021-03-19
- Fix strum dependency version

## [2.6.3] - 2021-03-19
- Fix probabilities when all biases are omitted.

## [2.6.1] - 2021-02-10
- Fix out of bounds error when replacement ends at the end of a contig.

## [2.6.0] - 2021-01-20
- Added read position bias into model.
- Performance improvements for exploration of bias events.
- Fixed accidental reporting of only the last variant in multi-allelic records.
- Fixed node-selector evaluation (e.g. G>A in scenario) for non SNV alleles.
- Fixed bug in replacement evaluation that could lead to artificially small probabilities due to incorrectly assembled alt alleles.

## [2.5.4] - 2021-01-04
- Added support for aux per base strand information as provided by rust-bio-tools.
- Fixed read orientation detection for overlapping reads.
- Better error messages for formulas.
- Fixed parsing of expression usage in formulas.

## [2.5.3] - 2020-11-23
- Update to latest hts-sys, containing a fix for macOS.

## [2.5.2] - 2020-11-20
- Fixed handling of missing insert size information. This is now detected automatically.
- CLI fixes.
- Adapt to htslib 0.35, thereby fixing some potential memory issues previously caused by unsafe reuse of buffers.

## [2.5.1] - 2020-11-17
- Various performance improvements (using jemalloc, avoiding bound checks that cannot fail, enabling more compiler optimizations, avoid fetching irrelevant reads).

## [2.5.0] - 2020-11-11
- Allow definition of re-usable expressions in scenarions (via a new key "expressions:", see https://varlociraptor.github.io).
- Remove ability to parallelize Varlociraptor via --threads. This never properly saturated the given cores and caused some overhead. Instead, we recommend to parallelize in scatter/gather style via `rbt vcf-split`, see https://varlociraptor.github.io).
- `resolution:` can now be skipped in scenarios. For continuous universes, this will then assume a resolution of 100, for discrete universes resolution isn't used anyway.

## [2.4.0] - 2020-11-05
- Allow scenarios to contain samples for which no BAM files are available. This allows to e.g. model tumor/normal from just the tumor sample with known contamination. Resulting probabilities will properly reflect the uncertainty about whether a variant is somatic or germline.
- Speed up SNV and MNV computations by precomputed call and miscall likelihoods.
- Add a flag --pairhmm-mode [exact|fast], that allows to instruct varlociraptor to only compute the optimal path in the pairHMM. This should be much faster in practice, but can come with some wrong likelihoods in rare extreme cases. Advice: only use on large cohorts or where exact allele frequencies do not matter.
- Fix insert size handling on single end samples (@dlaehnemann).

## [2.3.0] - 2020-09-09
- Include read orientation bias into the model.
- Exlcude softclipped and non-standard orientation reads from SNV and MNV calling as they are indicative of SVs and often cause artifact substitutions while sometimes not being reflected via a higher uncertainty in the MAPQ. Not considering them is the conservative choice.

## [2.2.1] - 2020-08-24
- Allow to set reference buffer size (--reference-buffer-size), for improved parallelization when calling SVs.
- Fix breakend handling getting confused between events when calling BCF with multiple breakends and more than a single thread.

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
