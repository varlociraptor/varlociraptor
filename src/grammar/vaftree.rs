use std::collections::HashSet;

use anyhow::Result;
use itertools::Itertools;

use crate::errors;
use crate::grammar::{formula::NormalizedFormula, formula::IUPAC, Scenario, VAFSpectrum};
use crate::variants::model::AlleleFreq;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct VAFTree {
    inner: Vec<Node>,
}

impl VAFTree {
    pub(crate) fn absent(n_samples: usize) -> Self {
        assert!(n_samples > 0, "bug: n_samples must be > 0");

        fn absent(sample: usize, n_samples: usize) -> Node {
            let children = if sample == n_samples - 1 {
                Vec::new()
            } else {
                vec![absent(sample + 1, n_samples)]
            };

            Node {
                kind: NodeKind::Sample {
                    sample,
                    vafs: VAFSpectrum::singleton(AlleleFreq(0.0)),
                },
                children,
            }
        }

        VAFTree {
            inner: vec![absent(0, n_samples)],
        }
    }
}

impl<'a> IntoIterator for &'a VAFTree {
    type Item = &'a Node;
    type IntoIter = std::slice::Iter<'a, Node>;

    fn into_iter(self) -> Self::IntoIter {
        (&self.inner).into_iter()
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum NodeKind {
    Variant {
        refbase: IUPAC,
        altbase: IUPAC,
        positive: bool,
    },
    Sample {
        sample: usize,
        vafs: VAFSpectrum,
    },
    False,
}

#[derive(new, Clone, Debug, PartialEq, Eq, Getters, Hash)]
#[get = "pub"]
pub(crate) struct Node {
    kind: NodeKind,
    #[new(default)]
    children: Vec<Node>,
}

impl Node {
    pub(crate) fn leafs<'a>(&'a mut self) -> Vec<&'a mut Node> {
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

    pub(crate) fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }

    pub(crate) fn is_branching(&self) -> bool {
        self.children.len() > 1
    }
}

impl VAFTree {
    pub(crate) fn new(
        formula: &NormalizedFormula,
        scenario: &Scenario,
        contig: &str,
    ) -> Result<Self> {
        fn from(formula: &NormalizedFormula, scenario: &Scenario) -> Result<Vec<Node>> {
            match formula {
                NormalizedFormula::Atom { sample, vafs } => {
                    let sample = scenario.idx(sample.as_str()).ok_or_else(|| {
                        errors::Error::InvalidSampleName {
                            name: sample.to_owned(),
                        }
                    })?;
                    Ok(vec![Node::new(NodeKind::Sample {
                        sample,
                        vafs: vafs.clone(),
                    })])
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
                        .sorted_by_key(|o| match o {
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
                &NormalizedFormula::Variant {
                    positive,
                    refbase,
                    altbase,
                } => Ok(vec![Node::new(NodeKind::Variant {
                    positive,
                    refbase,
                    altbase,
                })]),
                NormalizedFormula::False => Ok(vec![Node::new(NodeKind::False)]),
            }
        }

        fn add_missing_samples<'a>(
            node: &mut Node,
            seen: &mut HashSet<usize>,
            scenario: &'a Scenario,
            contig: &str,
        ) -> Result<()> {
            if let NodeKind::False = node.kind {
                // METHOD: no need to add further missing samples as the formula is false anyways
                return Ok(());
            }

            if let NodeKind::Sample { sample, .. } = node.kind {
                seen.insert(sample);
            }

            if node.is_leaf() {
                // leaf, add missing samples
                for (name, sample) in scenario.samples() {
                    let idx = scenario.idx(name).unwrap();
                    if !seen.contains(&idx) {
                        seen.insert(idx);

                        node.children = sample
                            .contig_universe(contig, scenario.species())?
                            .iter()
                            .map(|vafs| {
                                Node::new(NodeKind::Sample {
                                    sample: idx,
                                    vafs: vafs.clone(),
                                })
                            })
                            .collect();
                        add_missing_samples(node, seen, scenario, contig)?;
                        break;
                    }
                }
            } else {
                if node.is_branching() {
                    for child in &mut node.children[1..] {
                        add_missing_samples(child, &mut seen.clone(), scenario, contig)?;
                    }
                }
                add_missing_samples(&mut node.children[0], seen, scenario, contig)?;
            }

            Ok(())
        }

        let mut inner = from(formula, scenario)?;
        for node in &mut inner {
            let mut seen = HashSet::new();
            add_missing_samples(node, &mut seen, scenario, contig)?;
        }

        Ok(VAFTree { inner })
    }
}
