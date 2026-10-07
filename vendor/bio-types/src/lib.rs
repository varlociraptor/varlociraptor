#[cfg(feature = "serde")]
extern crate serde;
#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate derive_new;

extern crate regex;

#[cfg(feature = "phylogeny")]
extern crate petgraph;

pub mod alignment;
pub mod annot;
pub mod genome;
#[cfg(feature = "phylogeny")]
pub mod phylogeny;
pub mod sequence;
pub mod strand;
pub mod variant;
