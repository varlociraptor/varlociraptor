use std::collections::HashSet;
use std::collections::HashMap;

use itertools::Itertools;

use crate::grammar::{formula::NormalizedFormula, VAFSpectrum, Scenario};
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
    fn new(formula: &NormalizedFormula, scenario: &Scenario) -> Self {
        let sample_idx: HashMap<_, _> = scenario.samples().keys().enumerate().map(|(i, s)| (s.as_str(), i)).collect();

        fn from(formula: &NormalizedFormula, sample_idx: &HashMap<&str, usize>) -> Vec<Box<Node>> {
            match formula {
                NormalizedFormula::Atom { sample, vafs } => {
                    vec![Box::new(Node::new(*sample_idx.get(sample.as_str()).unwrap(), vafs.clone()))]
                }
                NormalizedFormula::Disjunction { operands } => operands
                    .iter()
                    .map(|operand| from(operand, sample_idx))
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
                        let subtrees = from(operand, sample_idx);
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

        fn add_missing_samples<'a>(
            node: &mut Node,
            seen: HashSet<&'a str>,
            scenario: &'a Scenario,
            sample_idx: &HashMap<&str, usize>
        ) {
            if node.is_leaf() {
                // leaf, add missing samples
                for (name, sample) in scenario.samples() {
                    if !seen.contains(name.as_str()) {
                        seen.insert(name.as_str());
                        node.children = sample
                            .universe()
                            .iter()
                            .map(|vafs| Box::new(Node::new(*sample_idx.get(name.as_str()).unwrap(), vafs.clone())))
                            .collect();
                        add_missing_samples(node, seen, scenario, sample_idx);
                    }
                }
            } else {
                if node.is_branching() {
                    for child in &mut node.children[1..] {
                        add_missing_samples(child.as_mut(), seen.clone(), scenario, sample_idx);
                    }
                }
                add_missing_samples(node.children[0].as_mut(), seen, scenario, sample_idx);
            }
        }

        let mut inner = from(formula, &sample_idx);
        for node in &mut inner {
            let seen = HashSet::new();
            add_missing_samples(node, seen, scenario, &sample_idx);
        }

        VAFTree { inner }
    }

    pub fn iter<'a>(&'a self) -> impl Iterator<Item = &'a Node> {
        self.inner.iter().map(|node| node.as_ref())
    }
}
