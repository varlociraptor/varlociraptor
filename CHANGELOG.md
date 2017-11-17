# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

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
