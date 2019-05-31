use std::collections::HashSet;
use vec_map::VecMap;

use itertools::Itertools;

use crate::grammar::{formula::NormalizedFormula, formula::VAFUniverse, VAFSpectrum};
use crate::model::AlleleFreq;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct VAFTree {
    inner: Vec<Box<Node>>,
}

impl VAFTree {
    pub fn absent(n_samples: usize) -> Self {
        assert!(n_samples > 0, "bug: n_samples must be > 0");

        fn absent(sample: usize, n_samples: usize) -> Box<Node> {
            let children = if sample == n_samples - 1 {
                Vec::new()
            } else {
                vec![absent(sample + 1, n_samples)]
            };

            Box::new(Node {
                sample,
                vafs: VAFSpectrum::singleton(AlleleFreq(0.0)),
                children,
            })
        }

        VAFTree {
            inner: vec![absent(0, n_samples)],
        }
    }
}

#[derive(new, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Getters)]
#[get = "pub"]
pub struct Node {
    sample: usize,
    vafs: VAFSpectrum,
    #[new(default)]
    children: Vec<Box<Node>>,
}

impl Node {
    pub fn leafs<'a>(&'a self) -> Vec<&'a Node> {
        fn collect_leafs<'a>(node: &'a Node, leafs: &mut Vec<&'a Node>) {
            if node.children.is_empty() {
                leafs.push(node);
            } else {
                for child in &node.children {
                    collect_leafs(child, leafs);
                }
            }
        }

        let mut leafs = Vec::new();
        collect_leafs(self, &mut leafs);
        leafs
    }

    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }

    pub fn is_branching(&self) -> bool {
        self.children.len() > 1
    }
}

impl VAFTree {
    fn new(formula: &NormalizedFormula<usize>, vaf_universe: &VecMap<VAFUniverse>) -> Self {
        fn from(formula: &NormalizedFormula<usize>) -> Vec<Box<Node>> {
            match formula {
                NormalizedFormula::Atom { sample, vafs } => {
                    vec![Box::new(Node::new(*sample, vafs.clone()))]
                }
                NormalizedFormula::Disjunction { operands } => operands
                    .iter()
                    .map(|operand| from(operand))
                    .flatten()
                    .collect_vec(),
                NormalizedFormula::Conjunction { operands } => {
                    // sort disjunctions to the end
                    let operands = operands.iter().sorted_by_key(|o| match o.as_ref() {
                        NormalizedFormula::Disjunction { .. } => 1,
                        _ => 0,
                    });
                    let mut leafs: Option<Vec<&Node>> = None;
                    let mut roots = None;
                    for operand in operands {
                        let subtrees = from(operand);
                        if let Some(leafs) = leafs {
                            for leaf in leafs {
                                leaf.children = subtrees.clone();
                            }
                        }
                        leafs = Some(
                            subtrees
                                .iter()
                                .map(|subtree| subtree.leafs())
                                .flatten()
                                .collect_vec(),
                        );
                        if roots.is_none() {
                            roots = Some(subtrees);
                        }
                    }
                    roots.unwrap()
                }
            }
        }

        fn add_missing_samples(
            node: &mut Node,
            seen: HashSet<usize>,
            vaf_universe: &VecMap<VAFUniverse>,
        ) {
            if node.is_leaf() {
                // leaf, add missing samples
                for (sample, universe) in vaf_universe {
                    if !seen.contains(&sample) {
                        seen.insert(sample);
                        node.children = universe
                            .iter()
                            .map(|vafs| Box::new(Node::new(sample, vafs.clone())))
                            .collect();
                        add_missing_samples(node, seen, vaf_universe);
                    }
                }
            } else {
                if node.is_branching() {
                    for child in &mut node.children[1..] {
                        add_missing_samples(child.as_mut(), seen.clone(), vaf_universe);
                    }
                }
                add_missing_samples(node.children[0].as_mut(), seen, vaf_universe);
            }
        }

        let mut inner = from(formula);
        for node in &mut inner {
            let seen = HashSet::new();
            add_missing_samples(node, seen, vaf_universe);
        }

        VAFTree { inner }
    }

    pub fn iter<'a>(&'a self) -> impl Iterator<Item = &'a Node> {
        self.inner.iter().map(|node| node.as_ref())
    }
}
