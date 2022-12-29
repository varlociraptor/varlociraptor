use std::collections::HashSet;

use anyhow::Result;
use itertools::Itertools;

use crate::errors;
use crate::grammar::{formula::Iupac, formula::NormalizedFormula, Scenario, VAFSpectrum};
use crate::utils::log2_fold_change::{Log2FoldChange, Log2FoldChangePredicate};
use crate::variants::model::modes::generic::{LikelihoodOperands, VafLfc};
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

    pub(crate) fn contains(&self, operands: &LikelihoodOperands) -> bool {
        self.inner.iter().any(|node| {
            let mut lfcs = operands.lfcs().iter().collect();
            node.contains(operands, &mut lfcs)
        })
    }
}

impl<'a> IntoIterator for &'a VAFTree {
    type Item = &'a Node;
    type IntoIter = std::slice::Iter<'a, Node>;

    fn into_iter(self) -> Self::IntoIter {
        (&self.inner).iter()
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum NodeKind {
    Variant {
        refbase: Iupac,
        altbase: Iupac,
        positive: bool,
    },
    Sample {
        sample: usize,
        vafs: VAFSpectrum,
    },
    Log2FoldChange {
        sample_a: usize,
        sample_b: usize,
        predicate: Log2FoldChangePredicate,
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
    pub(crate) fn leafs(&mut self) -> Vec<&mut Node> {
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

    pub(crate) fn contains(&self, operands: &LikelihoodOperands, lfcs: &mut Vec<&VafLfc>) -> bool {
        let contained = match &self.kind {
            NodeKind::Sample { sample, vafs } => {
                vafs.contains(operands.events().get(*sample).unwrap().allele_freq)
            }
            NodeKind::Log2FoldChange {
                sample_a,
                sample_b,
                predicate,
            } => {
                let mut lfc_found = false;
                lfcs.retain(|lfc| {
                    let found = lfc.sample_a() == sample_a
                        && lfc.sample_b() == sample_b
                        && lfc.predicate() == predicate;
                    lfc_found |= found;
                    !found
                });
                lfc_found
            }
            NodeKind::False => false,
            NodeKind::Variant { .. } => true,
        };
        if self.children.is_empty() {
            // leaf, hence all given lfcs have to be already visited, otherwise they aren't contained in the path
            contained && lfcs.is_empty()
        } else {
            contained
                && self.children.iter().any(|node| {
                    if self.children.len() == 1 {
                        node.contains(operands, lfcs)
                    } else {
                        let mut lfcs = lfcs.clone();
                        node.contains(operands, &mut lfcs)
                    }
                })
        }
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
                    let mut roots = from(operands[0], scenario)?;
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
                NormalizedFormula::Log2FoldChange {
                    sample_a,
                    sample_b,
                    predicate,
                } => {
                    let sample_a = scenario.idx(sample_a.as_str()).ok_or_else(|| {
                        errors::Error::InvalidSampleName {
                            name: sample_a.to_owned(),
                        }
                    })?;
                    let sample_b = scenario.idx(sample_b.as_str()).ok_or_else(|| {
                        errors::Error::InvalidSampleName {
                            name: sample_b.to_owned(),
                        }
                    })?;
                    Ok(vec![Node::new(NodeKind::Log2FoldChange {
                        sample_a,
                        sample_b,
                        predicate: *predicate,
                    })])
                }
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
