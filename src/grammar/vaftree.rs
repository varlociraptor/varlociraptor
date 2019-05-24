use crate::grammar::{VAFSpectrum, NormalizedFormula};

pub struct VAFTree {
    inner: Vec<Box<Node>>
}

#[derive(new)]
pub struct Node {
    sample: usize, vafs: VAFSpectrum, children: Vec<Box<Node>>
}

impl Node {
    pub fn leafs(&self) -> Vec<&Node> {
        fn leafs(node: &Node, leafs: &mut Vec<Node>) {
            if node.children.is_empty() {
                leafs.push(node);
            } else {
                for child in &node.children {
                    leafs(child, leafs);
                }
            }
        }

        let mut leafs = Vec::new();
        leafs(self, leafs);
        leafs
    }
}

impl VAFTree {
    fn new(formula: &NormalizedFormula<usize>) -> Self {
        fn from(formula: &NormalizedFormula<usize>) -> Vec<Box<Node>> {
            match formula {
                NormalizedFormula::Atom { sample, vafs } => {
                    vec![Box::new(Node { sample: sample, vafs: vafs.clone(), children: Vec::new() })]
                },
                NormalizedFormula::Disjunction { operands } => {
                    operands.iter().map(|operand| Box::new(from(operand))).collect()
                },
                NormalizedFormula::Conjunction { operands } => {
                    // sort disjunctions to the end
                    let operands = operands.sorted_by_key(
                        |o| match o { NormalizedFormula::Disjunction { .. } => 1, _ => 0 }
                    );
                    let mut leafs = None;
                    let mut roots = None;
                    for operand in operands {
                        let subtrees = from(operand);
                        if leafs.is_some() {
                            for leaf in leafs {
                                leaf.children = subtrees.clone();
                            }
                        }
                        leafs = Some(subtrees.iter().map(|subtree| subtree.leafs()).flatten().collect_vec());
                        if roots.is_none() {
                            roots = subtrees;
                        }
                    }
                    roots.unwrap()
                }
            }
        }
    }
}
