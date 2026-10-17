// Copyright 2020 Franklin Delehelle
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A phylogenetic tree is represented as a directed graph.
//! Each node is a taxon, identified as a string.
//! The edges are weighted by the phylogenetic distance if it was defined, or f32::NAN otherwise.

use petgraph::{graph::Graph, Directed};

pub type Taxon = String;
pub type Proximity = f32;

pub type TreeGraph = Graph<Taxon, Proximity, Directed>;
pub struct Tree {
    pub g: TreeGraph,
}
