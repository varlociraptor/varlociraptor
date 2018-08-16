# libprosic

[![Travis](https://img.shields.io/travis/PROSIC/libprosic.svg?maxAge=2592000?style=flat-square)](https://travis-ci.org/PROSIC/libprosic)
[![API docs](https://img.shields.io/badge/API-documentation-blue.svg)](https://docs.rs/libprosic)
[![Crates.io](https://img.shields.io/crates/d/libprosic.svg)](https://crates.io/crates/libprosic)
[![Commitizen friendly](https://img.shields.io/badge/commitizen-friendly-brightgreen.svg)](http://commitizen.github.io/cz-cli/)

A Rust library for calling genomic variants using a latent variable model. This library is under active development and can break at any time.

## Authors

### Model

* [Johannes Köster](https://github.com/johanneskoester) (sampling probability, indel allele probabilities, FDR control, prior probabilities)
* [Louis Dijkstra](https://github.com/louisdijkstra) (latent variable model)
* [David Lähnemann](https://github.com/dlaehnemann) (single cell whole genome amplification model, single cell & bulk joint calling model / event setup)
* [Alexander Schönhuth](https://github.com/aschoen) (latent variable model, single cell whole genome amplification model, single cell & bulk joint calling model / event setup)

### Implementation

* [Johannes Köster](https://github.com/johanneskoester) (framework, software architecture, model, tumor-normal variant calling)
* [David Lähnemann](https://github.com/dlaehnemann) (single cell variant calling, FDR on Event sets) with strong support from Johannes Köster

### Supervision

* Supervision of David Lähnemann: [Alice McHardy](https://github.com/alicemchardy) and [Alexander Schönhuth](https://github.com/aschoen)

### Affiliations

Affiliations during work on the project:

* [Life Sciences and Health group](https://www.cwi.nl/research/groups/life-sciences-and-health), Centrum Wiskunde & Informatica, Amsterdam, The Netherlands: Louis Dijkstra, Johannes Köster, Alexander Schönhuth
* [Computational Biology of Infection Research Group](https://www.helmholtz-hzi.de/en/research/research_topics/bacterial_and_viral_pathogens/computational_biology_of_infection_research/our_research/), Helmholtz Centre for Infection Research, Braunschweig, Germany: David Lähnemann, Alice McHardy
* [Algorithms for reproducible bioinformatics lab](https://koesterlab.github.io/), Institute of Human Genetics, University Hospital Essen, University of Duisburg-Essen, Germany: Johannes Köster
* [Department of Pediatric Oncology, Hematology, and Clinical Immunology](https://www.uniklinik-duesseldorf.de/en/unternehmen/kliniken/department-of-paediatric-oncology-haematology-and-immunology/), University Children’s Hospital, Medical Faculty, Heinrich Heine University, Düsseldorf, Germany: David Lähnemann
