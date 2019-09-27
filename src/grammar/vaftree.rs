use std::collections::HashSet;

use itertools::Itertools;

use crate::errors;
use crate::grammar::{formula::NormalizedFormula, Scenario, VAFSpectrum};
use crate::model::AlleleFreq;
use crate::errors::Result;

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
    pub fn leafs<'a>(&'a mut self) -> Vec<&'a mut Node> {
        fn collect_leafs<'a>(node: &'a mut Node, leafs: &mut Vec<&'a mut Node>) {
            if node.children.is_empty() {
                leafs.push(node);
            } else {
                for child in &mut node.children {
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
    pub fn new(
        formula: &NormalizedFormula,
        scenario: &Scenario,
    ) -> Result<Self> {
        fn from(
            formula: &NormalizedFormula,
            scenario: &Scenario,
        ) -> Result<Vec<Box<Node>>> {
            match formula {
                NormalizedFormula::Atom { sample, vafs } => {
                    let sample = scenario.idx(sample.as_str()).ok_or_else(|| {
                        errors::Error::InvalidSampleName {
                            name: sample.to_owned(),
                        }
                    })?;
                    Ok(vec![Box::new(Node::new(sample, vafs.clone()))])
                }
                NormalizedFormula::Disjunction { operands } => {
                    let mut subtrees = Vec::new();
                    for operand in operands {
                        for subtree in from(operand, scenario)? {
                            subtrees.push(subtree);
                        }
                    }
                    Ok(subtrees)
                }
                NormalizedFormula::Conjunction { operands } => {
                    // sort disjunctions to the end
                    let operands = operands
                        .iter()
                        .sorted_by_key(|o| match o.as_ref() {
                            NormalizedFormula::Disjunction { .. } => 1,
                            _ => 0,
                        })
                        .collect_vec();
                    let mut roots = from(&operands[0], scenario)?;
                    for operand in &operands[1..] {
                        let subtrees = from(operand, scenario)?;
                        for subtree in &mut roots {
                            for leaf in subtree.leafs() {
                                leaf.children = subtrees.clone();
                            }
                        }
                    }
                    Ok(roots)
                }
            }
        }

        fn add_missing_samples<'a>(
            node: &mut Node,
            seen: &mut HashSet<usize>,
            scenario: &'a Scenario,
        ) {
            seen.insert(node.sample);
            if node.is_leaf() {
                // leaf, add missing samples
                for (name, sample) in scenario.samples() {
                    let idx = scenario.idx(name).unwrap();
                    if !seen.contains(&idx) {
                        seen.insert(idx);
                        node.children = sample
                            .universe()
                            .iter()
                            .map(|vafs| Box::new(Node::new(idx, vafs.clone())))
                            .collect();
                        add_missing_samples(node, seen, scenario);
                        break;
                    }
                }
            } else {
                if node.is_branching() {
                    for child in &mut node.children[1..] {
                        add_missing_samples(child.as_mut(), &mut seen.clone(), scenario);
                    }
                }
                add_missing_samples(node.children[0].as_mut(), seen, scenario);
            }
        }

        let mut inner = from(formula, scenario)?;
        for node in &mut inner {
            let mut seen = HashSet::new();
            add_missing_samples(node, &mut seen, scenario);
        }

        Ok(VAFTree { inner })
    }

    pub fn iter<'a>(&'a self) -> impl Iterator<Item = &'a Node> {
        self.inner.iter().map(|node| node.as_ref())
    }
}
