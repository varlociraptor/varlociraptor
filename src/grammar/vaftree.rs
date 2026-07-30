use std::collections::btree_map::Entry;
use std::collections::{BTreeMap, BTreeSet, HashSet};

use anyhow::Result;
use itertools::Itertools;

use crate::errors;
use crate::grammar::{
    formula::Iupac, formula::NormalizedFormula, Scenario, VAFSpectrum, VAFUniverse,
};
use crate::utils::log2_fold_change::Log2FoldChangePredicate;
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

    pub(crate) fn contains(
        &self,
        operands: &LikelihoodOperands,
        exclude_sample: Option<usize>,
    ) -> bool {
        self.inner.iter().any(|node| {
            let mut lfcs = operands.lfcs().iter().collect();
            node.contains(operands, &mut lfcs, exclude_sample)
        })
    }
}

impl<'a> IntoIterator for &'a VAFTree {
    type Item = &'a Node;
    type IntoIter = std::slice::Iter<'a, Node>;

    fn into_iter(self) -> Self::IntoIter {
        self.inner.iter()
    }
}

/// A single conjunctive clause of a `VAFTree`, i.e. one root-to-leaf path. It captures exactly the
/// set of assignments (per-sample VAF plus variant / log2-fold-change context) for which that path
/// evaluates to true. The disjunction of all clauses of a tree is the event itself.
#[derive(Clone, Debug)]
pub(crate) struct Clause {
    /// Per-sample VAF constraint. Every sample of the scenario is present on every live path
    /// (see `VAFTree::new`'s `add_missing_samples`); a sample constrained more than once along the
    /// path carries the intersection of its constraints.
    vafs: BTreeMap<usize, VAFSpectrum>,
    /// Variant terminals (`refbase>altbase`, negated if `positive` is false) required on this path.
    variants: BTreeSet<(Iupac, Iupac, bool)>,
    /// Log2-fold-change terminals required on this path.
    lfcs: BTreeSet<(usize, usize, Log2FoldChangePredicate)>,
}

impl Clause {
    fn empty() -> Self {
        Clause {
            vafs: BTreeMap::new(),
            variants: BTreeSet::new(),
            lfcs: BTreeSet::new(),
        }
    }
}

/// Depth-first enumeration of a node's clauses, extending `acc` with the current node. Paths
/// through a `False` node or an empty VAF intersection are pruned (they can never be true).
fn collect_clauses(node: &Node, acc: &Clause, out: &mut Vec<Clause>) {
    let mut acc = acc.clone();
    match node.kind() {
        NodeKind::Sample { sample, vafs } => match acc.vafs.entry(*sample) {
            Entry::Occupied(mut entry) => {
                let intersection = entry.get().intersect(vafs);
                if intersection.is_empty() {
                    return;
                }
                *entry.get_mut() = intersection;
            }
            Entry::Vacant(entry) => {
                if vafs.is_empty() {
                    return;
                }
                entry.insert(vafs.clone());
            }
        },
        NodeKind::Variant {
            refbase,
            altbase,
            positive,
        } => {
            acc.variants.insert((*refbase, *altbase, *positive));
        }
        NodeKind::Log2FoldChange {
            sample_a,
            sample_b,
            predicate,
        } => {
            acc.lfcs.insert((*sample_a, *sample_b, *predicate));
        }
        NodeKind::True => {}
        NodeKind::False => return,
    }

    if node.is_leaf() {
        out.push(acc);
    } else {
        for child in node.children() {
            collect_clauses(child, &acc, out);
        }
    }
}

impl VAFTree {
    /// Enumerate the conjunctive clauses (root-to-leaf paths) of this tree, pruning unsatisfiable
    /// paths.
    fn clauses(&self) -> Vec<Clause> {
        let mut out = Vec::new();
        let root = Clause::empty();
        for node in &self.inner {
            collect_clauses(node, &root, &mut out);
        }
        out
    }

    /// If `self` and `other` are *not* disjoint, return a human-readable witness: a concrete VAF
    /// assignment (drawn from the universes) together with the shared variant / log2-fold-change
    /// context for which both events are true simultaneously. Returns `None` if the two events are
    /// disjoint.
    ///
    /// Two clauses can be jointly true only if they impose exactly the same variant and
    /// log2-fold-change terminals: this matches the model's own containment semantics
    /// (`VAFTree::contains`), where a leaf is satisfied only when the provided log2-fold-change
    /// facts are consumed exactly, and where a genomic locus carries a single, definite variant.
    /// For events that differ only in such terminals this is conservative (it never reports a
    /// spurious overlap), while for pure allele-frequency events the check is exact and complete.
    pub(crate) fn overlap_witness(
        &self,
        other: &VAFTree,
        universes: &[(String, VAFUniverse)],
    ) -> Option<String> {
        let own = self.clauses();
        let others = other.clauses();
        for a in &own {
            for b in &others {
                if a.variants != b.variants || a.lfcs != b.lfcs {
                    continue;
                }
                if let Some(witness) = clause_witness(a, b, universes) {
                    return Some(witness);
                }
            }
        }
        None
    }
}

/// Build a witness for a pair of clauses that share the same variant/lfc context, or `None` if for
/// some sample the two VAF constraints have no common value within that sample's universe.
fn clause_witness(a: &Clause, b: &Clause, universes: &[(String, VAFUniverse)]) -> Option<String> {
    let mut assignment = Vec::with_capacity(universes.len());
    for (sample, (name, universe)) in universes.iter().enumerate() {
        // Combine the (possibly absent) constraints of both clauses for this sample.
        let combined = match (a.vafs.get(&sample), b.vafs.get(&sample)) {
            (Some(x), Some(y)) => Some(x.intersect(y)),
            (Some(x), None) => Some(x.clone()),
            (None, Some(y)) => Some(y.clone()),
            (None, None) => None,
        };
        if let Some(ref spectrum) = combined {
            if spectrum.is_empty() {
                return None;
            }
        }

        // Clip against the universe (a disjunction of spectra) and pick a representative value. If
        // no universe value satisfies the combined constraint, this VAF combination can never
        // occur and the two events do not actually overlap here.
        let mut universe_specs: Vec<&VAFSpectrum> = universe.iter().collect();
        universe_specs.sort();
        let representative = universe_specs.iter().find_map(|u| match &combined {
            Some(spectrum) => spectrum.intersect(u).representative(),
            None => u.representative(),
        })?;
        assignment.push(format!("{}={}", name, representative));
    }

    let mut witness = assignment.join(", ");
    let context = clause_context(a, universes);
    if !context.is_empty() {
        witness = format!("{}, {}", witness, context.join(", "));
    }
    Some(witness)
}

/// Human-readable description of the variant and log2-fold-change terminals shared by a clause.
fn clause_context(clause: &Clause, universes: &[(String, VAFUniverse)]) -> Vec<String> {
    let name = |idx: usize| {
        universes
            .get(idx)
            .map(|(name, _)| name.clone())
            .unwrap_or_else(|| format!("sample{}", idx))
    };
    let mut context = Vec::new();
    for (refbase, altbase, positive) in &clause.variants {
        context.push(format!(
            "{}({}>{})",
            if *positive { "" } else { "!" },
            refbase,
            altbase
        ));
    }
    for (sample_a, sample_b, predicate) in &clause.lfcs {
        context.push(format!(
            "l2fc({}, {}) {}",
            name(*sample_a),
            name(*sample_b),
            predicate
        ));
    }
    context
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
    True,
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

    pub(crate) fn contains(
        &self,
        operands: &LikelihoodOperands,
        lfcs: &mut Vec<&VafLfc>,
        exclude_sample: Option<usize>,
    ) -> bool {
        let contained = match &self.kind {
            NodeKind::Sample { sample, vafs } => {
                if exclude_sample.map_or(false, |s| s == *sample) {
                    // if sample is excluded, always return true since it may not
                    // influence the containment check
                    return true;
                }
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
            NodeKind::True => true,
            NodeKind::Variant { .. } => true,
        };
        if self.children.is_empty() {
            // leaf, hence all given lfcs have to be already visited, otherwise they aren't contained in the path
            contained && lfcs.is_empty()
        } else {
            contained
                && self.children.iter().any(|node| {
                    if self.children.len() == 1 {
                        node.contains(operands, lfcs, exclude_sample)
                    } else {
                        let mut lfcs = lfcs.clone();
                        node.contains(operands, &mut lfcs, exclude_sample)
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
                NormalizedFormula::True => Ok(vec![Node::new(NodeKind::True)]),
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

        fn add_missing_samples(
            node: &mut Node,
            seen: &mut HashSet<usize>,
            scenario: &Scenario,
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
