# varlociraptor

[![Travis](https://img.shields.io/travis/varlociraptor/varlociraptor/master.svg?maxAge=2592000)](https://travis-ci.org/varlociraptor/varlociraptor)
[![Codecov](https://img.shields.io/codecov/c/github/varlociraptor/varlociraptor/master.svg)](https://codecov.io/gh/varlociraptor/varlociraptor)
[![API docs](https://img.shields.io/badge/API-documentation-blue.svg)](https://docs.rs/varlociraptor)
[![Crates.io](https://img.shields.io/crates/d/varlociraptor.svg)](https://crates.io/crates/varlociraptor)
[![Commitizen friendly](https://img.shields.io/badge/commitizen-friendly-brightgreen.svg)](http://commitizen.github.io/cz-cli/)

A Rust library for calling genomic variants using a latent variable model. This library is under active development. A preprint describing an early version of the model can be found here: https://doi.org/10.1101/121954.

## Authors

### Model

* [Johannes Köster](https://github.com/johanneskoester)
* [Louis Dijkstra](https://github.com/louisdijkstra)
* [David Lähnemann](https://github.com/dlaehnemann)
* [Tobias Marschall](https://github.com/tobiasmarschall)
* [Alexander Schönhuth](https://github.com/aschoen)

### Implementation

* [Johannes Köster](https://github.com/johanneskoester) (framework, software architecture, model, tumor-normal variant calling, generic variant calling, filtration, CNV calling)
* [David Lähnemann](https://github.com/dlaehnemann) (single cell variant calling, FDR on Event sets)

### Supervision

* Supervision of David Lähnemann: [Alice McHardy](https://github.com/alicemchardy), [Alexander Schönhuth](https://github.com/aschoen), [Johannes Köster](https://github.com/johanneskoester)

### Affiliations during work on the project

* [Life Sciences and Health group](https://www.cwi.nl/research/groups/life-sciences-and-health), Centrum Wiskunde & Informatica, Amsterdam, The Netherlands: Louis Dijkstra, Johannes Köster, Alexander Schönhuth
* [Computational Biology of Infection Research Group](https://www.helmholtz-hzi.de/en/research/research_topics/bacterial_and_viral_pathogens/computational_biology_of_infection_research/our_research/), Helmholtz Centre for Infection Research, Braunschweig, Germany: David Lähnemann, Alice McHardy
* [Algorithms for reproducible bioinformatics lab](https://koesterlab.github.io/), Institute of Human Genetics, University Hospital Essen, University of Duisburg-Essen, Germany: Johannes Köster, David Lähnemann
* [Department of Pediatric Oncology, Hematology, and Clinical Immunology](https://www.uniklinik-duesseldorf.de/en/unternehmen/kliniken/department-of-paediatric-oncology-haematology-and-immunology/), University Children’s Hospital, Medical Faculty, Heinrich Heine University, Düsseldorf, Germany: David Lähnemann

