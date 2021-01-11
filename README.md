# Varlociraptor

[![Bioconda](https://img.shields.io/conda/dn/bioconda/varlociraptor?label=bioconda%20downloads)](https://bioconda.github.io/recipes/varlociraptor/README.html)
![GitHub Workflow Status](https://img.shields.io/github/workflow/status/varlociraptor/varlociraptor/CI/master?label=tests)
[![Codecov](https://img.shields.io/codecov/c/github/varlociraptor/varlociraptor/master.svg?label=test%20coverage)](https://codecov.io/gh/varlociraptor/varlociraptor)
[![API docs](https://img.shields.io/badge/API-documentation-blue.svg)](https://docs.rs/varlociraptor)
[![Commitizen friendly](https://img.shields.io/badge/commitizen-friendly-brightgreen.svg)](http://commitizen.github.io/cz-cli/)

Varlociraptor implements a novel, unified fully uncertainty-aware approach to genomic variant calling in arbitrary scenarios. 

### Key features

* Calls SNVs, MNVs, indels, inversions, duplications, replacements and breakends in all length ranges (from small to structural) with a unified statistical model.
* The statistical model encompasses all possible sources of uncertainty and biases.
* Resulting variant calls can be filtered by false discovery rate. No parameter tuning necessary.
* Unbiased, maximum a posteriori allele frequency estimates are provided with each call.

### Calling modes

* Tumor-normal-calling, classifying variants as somatic in tumor, somatic in normal, germline, or absent.
* Generic, grammar based configuration of the statistical model, allowing to classify arbitrary scenarios, from poplation genetics, to pedigrees, complex tumor scenarios and arbitrary combinations thereof.

**For details, see the homepage: https://varlociraptor.github.io**
