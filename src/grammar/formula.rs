use std::cmp::{Ord, Ordering, PartialOrd};
use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};
use std::fmt;
use std::ops;

use anyhow::Result;
use boolean_expression::Expr;
use itertools::Itertools;
use pest::iterators::{Pair, Pairs};
use pest::Parser;
use serde::de;
use serde::Deserialize;

use crate::errors;
use crate::grammar::{ExpressionIdentifier, LogFoldChangePredicate, Scenario};
use crate::variants::model::AlleleFreq;

#[derive(Shrinkwrap, Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct Iupac(u8);

impl Iupac {
    pub(crate) fn contains(&self, base: u8) -> bool {
        if base == **self {
            return true;
        }
        match **self {
            b'R' if base == b'A' || base == b'G' => true,
            b'Y' if base == b'C' || base == b'T' => true,
            b'S' if base == b'G' || base == b'C' => true,
            b'W' if base == b'A' || base == b'T' => true,
            b'K' if base == b'G' || base == b'T' => true,
            b'M' if base == b'A' || base == b'C' => true,
            b'B' if base == b'C' || base == b'G' || base == b'T' => true,
            b'D' if base == b'A' || base == b'G' || base == b'T' => true,
            b'H' if base == b'A' || base == b'C' || base == b'T' => true,
            b'V' if base == b'A' || base == b'C' || base == b'G' => true,
            b'N' => true,
            _ => false,
        }
    }
}

#[derive(Parser)]
#[grammar = "grammar/formula.pest"]
pub(crate) struct FormulaParser;

#[derive(PartialEq, Eq, Clone, Debug, Hash, PartialOrd, Ord)]
pub(crate) enum Formula {
    Conjunction { operands: Vec<Formula> },
    Disjunction { operands: Vec<Formula> },
    Negation { operand: Box<Formula> },
    Terminal(FormulaTerminal),
}

impl From<NormalizedFormula> for Formula {
    fn from(formula: NormalizedFormula) -> Self {
        fn from_normalized(formula: NormalizedFormula) -> Formula {
            match formula {
                NormalizedFormula::Conjunction { operands } => Formula::Conjunction {
                    operands: operands.into_iter().map(from_normalized).collect(),
                },
                NormalizedFormula::Disjunction { operands } => Formula::Disjunction {
                    operands: operands.into_iter().map(from_normalized).collect(),
                },
                NormalizedFormula::Atom { sample, vafs } => {
                    Formula::Terminal(FormulaTerminal::Atom { sample, vafs })
                }
                NormalizedFormula::Variant {
                    altbase,
                    positive,
                    refbase,
                } => Formula::Terminal(FormulaTerminal::Variant {
                    altbase,
                    positive,
                    refbase,
                }),
                NormalizedFormula::False => Formula::Terminal(FormulaTerminal::False),
            }
        }
        from_normalized(formula)
    }
}

#[derive(PartialEq, Eq, Clone, Debug, Hash, PartialOrd, Ord)]
pub(crate) enum FormulaTerminal {
    Atom {
        sample: String,
        vafs: VAFSpectrum,
    },
    Variant {
        positive: bool,
        refbase: Iupac,
        altbase: Iupac,
    },
    Expression {
        identifier: ExpressionIdentifier,
        negated: bool,
    },
    LogFoldChange {
        sample_a: String,
        sample_b: String,
        value: LogFoldChangePredicate,
    },
    False,
}

impl FormulaTerminal {
    fn merge_conjunctions(&mut self, other: &FormulaTerminal) {
        match (self, other) {
            (
                FormulaTerminal::Atom {
                    sample: sample_a,
                    vafs: ref mut vafs_a,
                },
                FormulaTerminal::Atom {
                    sample: sample_b,
                    vafs: vafs_b,
                },
            ) if sample_a == sample_b => match (&vafs_a, &vafs_b) {
                (VAFSpectrum::Range(a), VAFSpectrum::Range(b)) => {
                    *vafs_a = VAFSpectrum::Range(a & b)
                }
                (VAFSpectrum::Range(a), VAFSpectrum::Set(b)) => {
                    *vafs_a = VAFSpectrum::Set(
                        b.iter().filter(|vaf| a.contains(**vaf)).cloned().collect(),
                    );
                }
                (VAFSpectrum::Set(a), VAFSpectrum::Range(b)) => {
                    *vafs_a = VAFSpectrum::Set(
                        a.iter().filter(|vaf| b.contains(**vaf)).cloned().collect(),
                    );
                }
                (VAFSpectrum::Set(a), VAFSpectrum::Set(b)) => {
                    *vafs_a = VAFSpectrum::Set(a.intersection(b).cloned().collect());
                }
            },
            _ => {
                panic!("bug: trying to merge FormulaTerminals that are not both atoms and for the same sample")
            }
        }
    }

    /// Try building the union of two instances of VAFSpectrum.
    /// Returns None in case ranges do not overlap or a set is not fully contained within a range,
    /// otherwise returns the union.
    fn try_merge_disjunction(&self, other: &FormulaTerminal) -> Option<VAFSpectrum> {
        match (self, other) {
            (
                FormulaTerminal::Atom {
                    sample: sample_a,
                    vafs: vafs_a,
                },
                FormulaTerminal::Atom {
                    sample: sample_b,
                    vafs: vafs_b,
                },
            ) if sample_a == sample_b => match (&vafs_a, &vafs_b) {
                (VAFSpectrum::Range(a), VAFSpectrum::Range(b)) => match a.overlap(b) {
                    VAFRangeOverlap::None => None,
                    _ => Some(VAFSpectrum::Range((a | b).0)),
                },
                (VAFSpectrum::Range(a), VAFSpectrum::Set(b)) => {
                    if b.iter().all(|v| a.contains(*v)) {
                        Some(VAFSpectrum::Range(a.clone()))
                    } else {
                        None
                    }
                }
                (VAFSpectrum::Set(a), VAFSpectrum::Range(b)) => {
                    if a.iter().all(|v| b.contains(*v)) {
                        Some(VAFSpectrum::Range(b.clone()))
                    } else {
                        None
                    }
                }
                (VAFSpectrum::Set(a), VAFSpectrum::Set(b)) => {
                    Some(VAFSpectrum::Set(a.union(b).cloned().collect()))
                }
            },
            _ => {
                panic!("bug: trying to merge FormulaTerminals that are not both atoms and for the same sample")
            }
        }
    }
}

impl std::fmt::Display for Formula {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let fmt_operand = |formula: &Formula| match formula {
            Formula::Terminal(_) => format!("{}", formula),
            Formula::Negation { operand } if operand.is_terminal() => format!("{}", formula),
            _ => format!("({})", formula),
        };

        let formatted = match self {
            Formula::Terminal(FormulaTerminal::False) => "false".to_owned(),
            Formula::Terminal(FormulaTerminal::Atom {
                sample,
                vafs: VAFSpectrum::Set(vafs),
            }) => {
                if self.is_terminal_false() {
                    "false".to_owned()
                } else {
                    vafs.iter()
                        .map(|vaf| format!("{}:{}", sample, vaf))
                        .join("|")
                }
            }
            Formula::Terminal(FormulaTerminal::Atom {
                sample,
                vafs: VAFSpectrum::Range(vafrange),
            }) => {
                let left_bracket = if vafrange.left_exclusive { ']' } else { '[' };
                let right_bracket = if vafrange.right_exclusive { '[' } else { ']' };
                format!(
                    "{}:{}{},{}{}",
                    sample, left_bracket, vafrange.start, vafrange.end, right_bracket
                )
            }
            Formula::Terminal(FormulaTerminal::Variant {
                positive,
                refbase,
                altbase,
            }) => format!(
                "{negate}({refbase}>{altbase})",
                negate = if *positive { "" } else { "!" },
                refbase = **refbase,
                altbase = **altbase,
            ),
            Formula::Terminal(FormulaTerminal::Expression {
                identifier,
                negated,
            }) => format!(
                "{negate}${expr}",
                negate = if *negated { "!" } else { "" },
                expr = **identifier
            ),
            Formula::Negation { operand } => format!("!{operand}", operand = fmt_operand(operand)),
            Formula::Conjunction { operands } => operands.iter().map(&fmt_operand).join(" & "),
            Formula::Disjunction { operands } => operands.iter().map(&fmt_operand).join(" | "),
        };
        write!(f, "{}", formatted)
    }
}

impl From<Expr<FormulaTerminal>> for Formula {
    fn from(expr: Expr<FormulaTerminal>) -> Self {
        match expr {
            Expr::Terminal(terminal) => Formula::Terminal(terminal),
            Expr::Not(a) => Formula::Negation {
                operand: Box::new((*a).into()),
            },
            Expr::And(a, b) => match ((*a).into(), (*b).into()) {
                (
                    Formula::Conjunction {
                        operands: mut left_operands,
                    },
                    Formula::Conjunction {
                        operands: right_operands,
                    },
                ) => {
                    left_operands.extend(right_operands);
                    Formula::Conjunction {
                        operands: left_operands,
                    }
                }
                (
                    Formula::Conjunction {
                        operands: mut left_operands,
                    },
                    right,
                ) => {
                    left_operands.push(right);
                    Formula::Conjunction {
                        operands: left_operands,
                    }
                }
                (
                    left,
                    Formula::Conjunction {
                        operands: mut right_operands,
                    },
                ) => {
                    right_operands.push(left);
                    Formula::Conjunction {
                        operands: right_operands,
                    }
                }
                (left, right) => Formula::Conjunction {
                    operands: vec![left, right],
                },
            },
            Expr::Or(a, b) => match ((*a).into(), (*b).into()) {
                (
                    Formula::Disjunction {
                        operands: mut left_operands,
                    },
                    Formula::Disjunction {
                        operands: right_operands,
                    },
                ) => {
                    left_operands.extend(right_operands);
                    Formula::Disjunction {
                        operands: left_operands,
                    }
                }
                (
                    Formula::Disjunction {
                        operands: mut left_operands,
                    },
                    right,
                ) => {
                    left_operands.push(right);
                    Formula::Disjunction {
                        operands: left_operands,
                    }
                }
                (
                    left,
                    Formula::Disjunction {
                        operands: mut right_operands,
                    },
                ) => {
                    right_operands.push(left);
                    Formula::Disjunction {
                        operands: right_operands,
                    }
                }
                (left, right) => Formula::Disjunction {
                    operands: vec![left, right],
                },
            },
            Expr::Const(false) => Formula::Terminal(FormulaTerminal::False),
            _ => panic!("bug: unexpected boolean expression containing constant"),
        }
    }
}

impl Into<Expr<FormulaTerminal>> for Formula {
    fn into(self) -> Expr<FormulaTerminal> {
        if self.is_terminal_false() {
            return Expr::Const(false);
        }
        match self {
            Formula::Terminal(terminal) => Expr::Terminal(terminal),
            Formula::Conjunction { mut operands } => {
                let mut expr = operands.pop().unwrap().into();
                for operand in operands {
                    expr &= operand.into();
                }
                expr
            }
            Formula::Disjunction { mut operands } => {
                let mut expr = operands.pop().unwrap().into();
                for operand in operands {
                    expr |= operand.into();
                }
                expr
            }
            Formula::Negation { operand } => Expr::Not(Box::new((*operand).into())),
        }
    }
}

impl Formula {
    pub(crate) fn is_terminal(&self) -> bool {
        matches!(self, Formula::Terminal(_))
    }

    pub(crate) fn to_terminal(&self) -> Option<&FormulaTerminal> {
        if let Formula::Terminal(terminal) = self {
            Some(terminal)
        } else {
            None
        }
    }

    /// Return true if this formula is a terminal that is always false (an empty VAF set).
    pub(crate) fn is_terminal_false(&self) -> bool {
        match self {
            Formula::Terminal(FormulaTerminal::Atom {
                vafs: VAFSpectrum::Set(vafs),
                ..
            }) => vafs.is_empty(),
            Formula::Terminal(FormulaTerminal::False) => true,
            _ => false,
        }
    }

    pub(crate) fn into_terminal(self) -> Option<FormulaTerminal> {
        if let Formula::Terminal(terminal) = self {
            Some(terminal)
        } else {
            None
        }
    }

    /// Generate formula representing the absent event
    pub(crate) fn absent(scenario: &Scenario) -> Self {
        Formula::Conjunction {
            operands: scenario
                .samples()
                .keys()
                .map(|sample| {
                    Formula::Terminal(FormulaTerminal::Atom {
                        sample: sample.to_owned(),
                        vafs: VAFSpectrum::singleton(AlleleFreq(0.0)),
                    })
                })
                .collect(),
        }
    }

    pub(crate) fn normalize(&self, scenario: &Scenario, contig: &str) -> Result<NormalizedFormula> {
        // METHOD: Expand all expressions and move negations down to atoms. Then, simplify via BDDs,
        // merge atoms (VAF intervals) of same sample in the same conjuction, and simplify again.
        let mut simplified = self
            .expand_expressions(scenario)?
            .apply_negations(scenario, contig)?
            .simplify()
            .merge_atoms()
            .simplify();
        simplified.strip_false();
        Ok(simplified.into_normalized_formula())
    }

    fn expand_expressions(&self, scenario: &Scenario) -> Result<Self> {
        Ok(match self {
            Formula::Conjunction { operands } => Formula::Conjunction {
                operands: operands
                    .iter()
                    .map(|operand| operand.expand_expressions(scenario))
                    .collect::<Result<Vec<_>>>()?,
            },
            Formula::Disjunction { operands } => Formula::Disjunction {
                operands: operands
                    .iter()
                    .map(|operand| operand.expand_expressions(scenario))
                    .collect::<Result<Vec<_>>>()?,
            },
            Formula::Negation { operand } => Formula::Negation {
                operand: Box::new(operand.expand_expressions(scenario)?),
            },
            Formula::Terminal(terminal) => {
                if let FormulaTerminal::Expression {
                    identifier,
                    negated,
                } = terminal
                {
                    if let Some(formula) = scenario.expressions().get(identifier) {
                        if *negated {
                            Formula::Negation {
                                operand: Box::new(formula.clone()),
                            }
                        } else {
                            formula.clone()
                        }
                    } else {
                        return Err(errors::Error::UndefinedExpression {
                            identifier: identifier.to_string(),
                        }
                        .into());
                    }
                } else {
                    Formula::Terminal(terminal.clone())
                }
            }
        })
    }

    fn into_normalized_formula(&self) -> NormalizedFormula {
        match self {
            Formula::Terminal(FormulaTerminal::Atom { sample, vafs }) => NormalizedFormula::Atom {
                sample: sample.to_owned(),
                vafs: vafs.to_owned(),
            },
            Formula::Conjunction { operands } => NormalizedFormula::Conjunction {
                operands: operands
                    .iter()
                    .map(|o| o.into_normalized_formula())
                    .collect(),
            },
            Formula::Disjunction { operands } => NormalizedFormula::Disjunction {
                operands: operands
                    .iter()
                    .map(|o| o.into_normalized_formula())
                    .collect(),
            },
            &Formula::Terminal(FormulaTerminal::Variant {
                positive,
                refbase,
                altbase,
            }) => NormalizedFormula::Variant {
                positive,
                refbase,
                altbase,
            },
            &Formula::Terminal(FormulaTerminal::Expression {
                identifier: _,
                negated: _,
            }) => {
                panic!("bug: expressions should be expanded before normalization");
            }
            Formula::Negation { operand: _ } => {
                panic!("bug: negations should have been applied before normalization")
            }
            Formula::Terminal(FormulaTerminal::False) => NormalizedFormula::False,
        }
    }

    fn merge_atoms(&self) -> Self {
        let group_operands = |operands: &Vec<Formula>| -> HashMap<Option<String>, Vec<Formula>> {
            operands.iter().cloned().into_group_map_by(|operand| {
                if let Formula::Terminal(FormulaTerminal::Atom { sample, vafs: _ }) = operand {
                    Some(sample.to_owned())
                } else {
                    // group all non-atoms together
                    None
                }
            })
        };
        match self {
            Formula::Conjunction { operands } => {
                // collect statements per sample
                let mut grouped_operands = group_operands(operands);

                // merge atoms of the same sample
                for (sample, statements) in &mut grouped_operands {
                    if let Some(_sample) = sample {
                        let mut merged_statement =
                            statements.pop().unwrap().into_terminal().unwrap();
                        for statement in statements.iter() {
                            merged_statement.merge_conjunctions(statement.to_terminal().unwrap());
                        }
                        *statements = vec![Formula::Terminal(merged_statement)];
                    } else {
                        continue;
                    }
                }

                Formula::Conjunction {
                    operands: grouped_operands
                        .into_iter()
                        .map(|(_, statements)| statements)
                        .flatten()
                        .collect(),
                }
            }
            Formula::Disjunction { operands } => {
                // collect statements per sample
                let mut grouped_operands = group_operands(operands);

                // merge atoms of the same sample
                for (sample, statements) in &mut grouped_operands {
                    if let Some(sample) = sample {
                        // Sort by start position of VAFRange or minimum of VAFSet.
                        // The idea is to try to keep merging neighbouring ranges/sets, greedily.
                        statements.sort_unstable_by_key(|stmt| {
                            if let Formula::Terminal(FormulaTerminal::Atom { sample: _, vafs }) =
                                stmt
                            {
                                match vafs {
                                    VAFSpectrum::Set(s) => s.iter().min().copied(),
                                    VAFSpectrum::Range(r) => Some(r.start),
                                }
                            } else {
                                None
                            }
                        });
                        let mut merged_statements = vec![];
                        // Pick off the first terminal from the list of statements..
                        let mut current_statement = statements.remove(0).into_terminal().unwrap();
                        for statement in statements.iter() {
                            // then look at the (current) next one if the merge was
                            // successful. Otherwise, `other_statement` will be the last statement
                            // that could *not* be merged with `current_statement`
                            let other_statement = statement.to_terminal().unwrap();

                            // Try merging (i.e. unionise!) the two spectra
                            if let Some(merged) =
                                current_statement.try_merge_disjunction(other_statement)
                            {
                                // if it succeeds, use the merge result as `current_statement`
                                current_statement = FormulaTerminal::Atom {
                                    sample: sample.into(),
                                    vafs: merged,
                                };
                            } else {
                                // if it fails, stash our `current_statement` and replace it with
                                // the other statement which couldn't be merged with
                                // `current_statement`.
                                merged_statements.push(Formula::Terminal(current_statement));
                                current_statement = other_statement.clone();
                            }
                        }
                        // don't forget to push the last result.
                        merged_statements.push(Formula::Terminal(current_statement));
                        *statements = merged_statements;
                    } else {
                        continue;
                    }
                }
                Formula::Disjunction {
                    operands: grouped_operands
                        .into_iter()
                        .map(|(_, statements)| statements)
                        .flatten()
                        .collect(),
                }
            }
            Formula::Negation { operand } => Formula::Negation {
                operand: Box::new(operand.merge_atoms()),
            },
            terminal => terminal.clone(),
        }
    }

    fn strip_false(&mut self) {
        if let Formula::Disjunction { ref mut operands } = self {
            *operands = operands
                .iter()
                .filter(|operand| {
                    if operand.is_terminal_false() {
                        false
                    } else if let Formula::Conjunction { operands } = operand {
                        !operands.iter().any(|operand| operand.is_terminal_false())
                    } else {
                        true
                    }
                })
                .cloned()
                .collect();
        }
    }

    /// Simplify formula via a BDD
    fn simplify(self) -> Self {
        let expr: Expr<FormulaTerminal> = self.into();
        let simplified = expr.simplify_via_bdd();
        simplified.into()
    }

    /// Negate formula.
    fn negate(&self, scenario: &Scenario, contig: &str) -> Result<Self> {
        Ok(match self {
            Formula::Terminal(FormulaTerminal::False) => {
                panic!("bug: negation not implemented for false terminal (this is unexpected since the grammar does not allow to specify false).")
            }
            Formula::Conjunction { operands } => Formula::Disjunction {
                operands: operands
                    .iter()
                    .map(|o| o.negate(scenario, contig))
                    .collect::<Result<Vec<Formula>>>()?,
            },
            Formula::Disjunction { operands } => Formula::Conjunction {
                operands: operands
                    .iter()
                    .map(|o| o.negate(scenario, contig))
                    .collect::<Result<Vec<Formula>>>()?,
            },
            Formula::Negation { operand } => operand.as_ref().clone(),
            &Formula::Terminal(FormulaTerminal::Variant {
                positive,
                refbase,
                altbase,
            }) => Formula::Terminal(FormulaTerminal::Variant {
                positive: !positive,
                refbase,
                altbase,
            }),
            Formula::Terminal(FormulaTerminal::Expression {
                identifier,
                negated,
            }) => Formula::Terminal(FormulaTerminal::Expression {
                identifier: identifier.clone(),
                negated: !negated,
            }),
            Formula::Terminal(FormulaTerminal::Atom { sample, vafs }) => {
                let universe = scenario
                    .samples()
                    .get(sample)
                    .ok_or_else(|| errors::Error::InvalidSampleName {
                        name: sample.to_owned(),
                    })?
                    .contig_universe(contig, scenario.species())?;

                let mut disjunction = Vec::new();
                match vafs {
                    VAFSpectrum::Set(vafs) => {
                        let mut uvaf_stack: VecDeque<_> = universe.iter().cloned().collect();
                        while let Some(uvafs) = uvaf_stack.pop_front() {
                            match uvafs {
                                VAFSpectrum::Set(uvafs) => {
                                    let difference: BTreeSet<_> =
                                        uvafs.difference(vafs).cloned().collect();
                                    if !difference.is_empty() {
                                        disjunction.push(VAFSpectrum::Set(difference));
                                    }
                                }
                                VAFSpectrum::Range(urange) => {
                                    for &vaf in vafs {
                                        if urange.contains(vaf) {
                                            let (left_urange, right_urange) = urange.split_at(vaf);
                                            if let Some(right_urange) = right_urange {
                                                uvaf_stack.push_back(right_urange);
                                            }
                                            if let Some(left_urange) = left_urange {
                                                disjunction.push(left_urange);
                                            }
                                        } else {
                                            disjunction.push(VAFSpectrum::Range(urange.clone()));
                                        }
                                    }
                                }
                            }
                        }
                    }
                    VAFSpectrum::Range(range) => {
                        for uvafs in universe.iter() {
                            match uvafs {
                                VAFSpectrum::Set(uvafs) => {
                                    let set: BTreeSet<_> = uvafs
                                        .iter()
                                        .filter(|uvaf| !range.contains(**uvaf))
                                        .cloned()
                                        .collect();
                                    if !set.is_empty() {
                                        disjunction.push(VAFSpectrum::Set(set));
                                    }
                                }
                                VAFSpectrum::Range(urange) => match range.overlap(urange) {
                                    VAFRangeOverlap::Equal => {
                                        // range is already covered entirely, nothing to add
                                    }
                                    VAFRangeOverlap::Contained => {
                                        if let Some(left) = urange.split_at(range.start).0 {
                                            disjunction.push(left);
                                        }
                                        if let Some(right) = urange.split_at(range.end).1 {
                                            disjunction.push(right);
                                        }
                                    }
                                    VAFRangeOverlap::End => {
                                        if let Some(spec) = urange.split_at(range.end).1 {
                                            disjunction.push(spec);
                                        }
                                    }
                                    VAFRangeOverlap::Start => {
                                        if let Some(spec) = urange.split_at(range.start).0 {
                                            disjunction.push(spec);
                                        }
                                    }
                                    VAFRangeOverlap::None => {
                                        disjunction.push(VAFSpectrum::Range(urange.clone()))
                                    }
                                    VAFRangeOverlap::Contains => (),
                                },
                            }
                        }
                    }
                }

                if disjunction.is_empty() {
                    // impossible, return empty set
                    Formula::Terminal(FormulaTerminal::Atom {
                        sample: sample.to_owned(),
                        vafs: VAFSpectrum::empty(),
                    })
                } else {
                    Formula::Disjunction {
                        operands: disjunction
                            .into_iter()
                            .map(|vafs| {
                                Formula::Terminal(FormulaTerminal::Atom {
                                    sample: sample.clone(),
                                    vafs,
                                })
                            })
                            .collect(),
                    }
                }
            }
        })
    }

    fn apply_negations(&self, scenario: &Scenario, contig: &str) -> Result<Self> {
        Ok(match self {
            Formula::Negation { operand } => operand
                .negate(scenario, contig)?
                .apply_negations(scenario, contig)?,
            Formula::Terminal(FormulaTerminal::Atom { sample, vafs }) => {
                Formula::Terminal(FormulaTerminal::Atom {
                    sample: sample.to_owned(),
                    vafs: vafs.to_owned(),
                })
            }
            Formula::Conjunction { operands } => {
                let operands = operands
                    .iter()
                    .map(|o| o.apply_negations(scenario, contig))
                    .collect::<Result<Vec<Formula>>>()?;

                Formula::Conjunction { operands }
            }
            Formula::Disjunction { operands } => Formula::Disjunction {
                operands: operands
                    .iter()
                    .map(|o| o.apply_negations(scenario, contig))
                    .collect::<Result<Vec<Formula>>>()?,
            },
            &Formula::Terminal(FormulaTerminal::Variant {
                positive,
                refbase,
                altbase,
            }) => Formula::Terminal(FormulaTerminal::Variant {
                positive,
                refbase,
                altbase,
            }),
            &Formula::Terminal(FormulaTerminal::Expression {
                identifier: _,
                negated: _,
            }) => {
                panic!("bug: expressions should be expanded before applying negations");
            }
            Formula::Terminal(FormulaTerminal::False) => {
                panic!("bug: false terminals may not appear in formula to be negated because this is not allowed in the grammar");
            }
        })
    }
}

#[derive(PartialEq, Eq, Clone, Debug, Hash, PartialOrd, Ord)]
pub(crate) enum NormalizedFormula {
    Conjunction {
        operands: Vec<NormalizedFormula>,
    },
    Disjunction {
        operands: Vec<NormalizedFormula>,
    },
    Atom {
        sample: String,
        vafs: VAFSpectrum,
    },
    Variant {
        positive: bool,
        refbase: Iupac,
        altbase: Iupac,
    },
    LogFoldChange {
        sample_a: String,
        sample_b: String,
        predicate: LogFoldChangePredicate,
    },
    False,
}

impl std::fmt::Display for NormalizedFormula {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let fmt_operand = |formula: &NormalizedFormula| match formula {
            NormalizedFormula::Atom { .. } | NormalizedFormula::Variant { .. } => {
                format!("{}", formula)
            }
            _ => format!("({})", formula),
        };

        let formatted = match self {
            NormalizedFormula::Atom {
                sample,
                vafs: VAFSpectrum::Set(vafs),
            } => match vafs.len() {
                1 => format!("{}:{}", sample, vafs.iter().next().unwrap()),
                x if x > 1 => format!(
                    "{}:{{{}}}",
                    sample,
                    vafs.iter().map(|vaf| format!("{:.3}", vaf)).join(", "),
                ),
                _ => "false".to_owned(),
            },
            NormalizedFormula::Atom {
                sample,
                vafs: VAFSpectrum::Range(vafrange),
            } => {
                let left_bracket = if vafrange.left_exclusive { ']' } else { '[' };
                let right_bracket = if vafrange.right_exclusive { '[' } else { ']' };
                format!(
                    "{}:{}{:.3},{:.3}{}",
                    sample, left_bracket, vafrange.start, vafrange.end, right_bracket
                )
            }
            NormalizedFormula::Variant {
                positive,
                refbase,
                altbase,
            } => format!(
                "{negate}({refbase}>{altbase})",
                negate = if *positive { "" } else { "!" },
                refbase = **refbase,
                altbase = **altbase,
            ),
            NormalizedFormula::Conjunction { operands } => {
                operands.iter().map(&fmt_operand).join(" & ")
            }
            NormalizedFormula::Disjunction { operands } => {
                operands.iter().map(&fmt_operand).join(" | ")
            }
            NormalizedFormula::False => "false".to_owned(),
        };
        write!(f, "{}", formatted)
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) enum VAFSpectrum {
    Set(BTreeSet<AlleleFreq>),
    Range(VAFRange),
}

impl VAFSpectrum {
    pub(crate) fn singleton(vaf: AlleleFreq) -> Self {
        let mut set = BTreeSet::new();
        set.insert(vaf);
        VAFSpectrum::Set(set)
    }

    pub(crate) fn empty() -> Self {
        VAFSpectrum::Set(BTreeSet::new())
    }

    pub(crate) fn contains(&self, vaf: AlleleFreq) -> bool {
        match self {
            VAFSpectrum::Set(ref set) => set.contains(&vaf),
            VAFSpectrum::Range(ref range) => range.contains(vaf),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, TypedBuilder, Hash)]
pub(crate) struct VAFRange {
    inner: ops::Range<AlleleFreq>,
    left_exclusive: bool,
    right_exclusive: bool,
}

#[derive(Debug, PartialEq, Eq, Copy, Clone, Hash)]
pub(crate) enum VAFRangeOverlap {
    Contained,
    Contains,
    End,
    Start,
    Equal,
    None,
}

impl VAFRange {
    pub(crate) fn empty() -> Self {
        Self {
            inner: AlleleFreq(0.0)..AlleleFreq(0.0),
            left_exclusive: true,
            right_exclusive: true,
        }
    }

    pub(crate) fn contains(&self, vaf: AlleleFreq) -> bool {
        match (self.left_exclusive, self.right_exclusive) {
            (true, true) => self.start < vaf && self.end > vaf,
            (true, false) => self.start < vaf && self.end >= vaf,
            (false, true) => self.start <= vaf && self.end > vaf,
            (false, false) => self.start <= vaf && self.end >= vaf,
        }
    }

    pub(crate) fn split_at(&self, vaf: AlleleFreq) -> (Option<VAFSpectrum>, Option<VAFSpectrum>) {
        assert!(
            self.contains(vaf),
            "bug: split_at is only defined if given VAF is contained in range"
        );
        let left = VAFRange {
            inner: self.start..vaf,
            left_exclusive: self.left_exclusive,
            right_exclusive: true,
        };
        let right = VAFRange {
            inner: vaf..self.end,
            left_exclusive: true,
            right_exclusive: self.right_exclusive,
        };

        let to_spectrum = |range: VAFRange| {
            if range.start == range.end {
                if !(range.left_exclusive && self.right_exclusive) {
                    Some(VAFSpectrum::singleton(range.start))
                } else {
                    None
                }
            } else {
                Some(VAFSpectrum::Range(range))
            }
        };

        (to_spectrum(left), to_spectrum(right))
    }

    pub(crate) fn overlap(&self, vafs: &VAFRange) -> VAFRangeOverlap {
        if self == vafs {
            return VAFRangeOverlap::Equal;
        }
        let range = self;
        let other_range = vafs;
        let start_is_right_of_start = match (self.left_exclusive, other_range.left_exclusive) {
            (true, true) => self.start > other_range.start,
            (true, false) => self.start >= other_range.start,
            (false, true) => self.start > other_range.start,
            (false, false) => self.start > other_range.start,
        };
        let end_is_left_of_end = match (self.right_exclusive, other_range.right_exclusive) {
            (true, true) => range.end < other_range.end,
            (true, false) => range.end <= other_range.end,
            (false, true) => range.end < other_range.end,
            (false, false) => range.end < other_range.end,
        };
        if (range.end < other_range.start || range.start > other_range.end)
            || (range.end <= other_range.start
                && (range.right_exclusive || other_range.left_exclusive))
            || (range.start >= other_range.end
                && (range.left_exclusive || other_range.right_exclusive))
        {
            VAFRangeOverlap::None
        } else {
            match (start_is_right_of_start, end_is_left_of_end) {
                (true, true) => VAFRangeOverlap::Contained,
                (true, false) => VAFRangeOverlap::Start,
                (false, true) => VAFRangeOverlap::End,
                (false, false) => VAFRangeOverlap::Contains,
            }
        }
    }

    pub(crate) fn observable_min(&self, n_obs: usize) -> AlleleFreq {
        let min_vaf = if n_obs < 10 {
            self.start
        } else {
            let obs_count = Self::expected_observation_count(self.start, n_obs);
            let adjust_allelefreq = |obs_count: f64| AlleleFreq(obs_count.ceil() / n_obs as f64);

            if self.left_exclusive && obs_count % 1.0 == 0.0 {
                // We are left exclusive and need to find a supremum from the right.

                let adjusted_end = self.observable_max(n_obs);

                for offset in &[1.0, 0.0] {
                    let adjusted_obs_count = obs_count + offset;
                    let adjusted_start = adjust_allelefreq(adjusted_obs_count);
                    if *adjusted_start <= 1.0 && adjusted_start <= adjusted_end {
                        return adjusted_start;
                    }
                }
            }

            adjust_allelefreq(obs_count)
        };
        if min_vaf >= self.observable_max(n_obs) {
            // If the adjustments destroys the order of the boundaries, we don't do it.
            // This can happen if the two boundaries are close together and we have only few observations.
            self.start
        } else {
            min_vaf
        }
    }

    pub(crate) fn observable_max(&self, n_obs: usize) -> AlleleFreq {
        assert!(
            *self.end != 0.0,
            "bug: observable_max may not be called if end=0.0."
        );
        if n_obs < 10 {
            self.end
        } else {
            let mut obs_count = Self::expected_observation_count(self.end, n_obs);
            if self.right_exclusive && obs_count % 1.0 == 0.0 {
                obs_count -= 1.0;
            }
            obs_count = obs_count.floor();
            if obs_count == 0.0 {
                // too few observations to handle exclusiveness
                self.end
            } else {
                AlleleFreq(obs_count.floor() / n_obs as f64)
            }
        }
    }

    fn expected_observation_count(freq: AlleleFreq, n_obs: usize) -> f64 {
        n_obs as f64 * *freq
    }
}

use auto_ops::impl_op_ex;

impl_op_ex!(&|a: &VAFRange, b: &VAFRange| -> VAFRange {
    match a.overlap(b) {
        VAFRangeOverlap::Contained => a.clone(),
        VAFRangeOverlap::Contains => b.clone(),
        VAFRangeOverlap::Start => VAFRange {
            inner: a.inner.start..b.inner.end,
            left_exclusive: a.left_exclusive,
            right_exclusive: b.right_exclusive,
        },
        VAFRangeOverlap::End => VAFRange {
            inner: b.inner.start..a.inner.end,
            left_exclusive: b.left_exclusive,
            right_exclusive: a.right_exclusive,
        },
        VAFRangeOverlap::Equal => a.clone(),
        VAFRangeOverlap::None => VAFRange::empty(),
    }
});

impl_op_ex!(| |a: &VAFRange, b: &VAFRange| -> (VAFRange, Option<VAFRange>) {
    match a.overlap(b) {
        VAFRangeOverlap::Contained => (b.clone(), None),
        VAFRangeOverlap::Contains => (a.clone(), None),
        VAFRangeOverlap::Start => (VAFRange {
            inner: b.inner.start..a.inner.end,
            left_exclusive: b.left_exclusive,
            right_exclusive: a.right_exclusive,
        }, None),
        VAFRangeOverlap::End => (VAFRange {
            inner: a.inner.start..b.inner.end,
            left_exclusive: a.left_exclusive,
            right_exclusive: b.right_exclusive,
        }, None),
        VAFRangeOverlap::Equal => (a.clone(), None),
        VAFRangeOverlap::None => (a.clone(), Some(b.clone()))
    }
});

impl ops::Deref for VAFRange {
    type Target = ops::Range<AlleleFreq>;

    fn deref(&self) -> &ops::Range<AlleleFreq> {
        &self.inner
    }
}

impl Ord for VAFRange {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.start.cmp(&other.start) {
            Ordering::Equal => self.end.cmp(&other.end),
            ord => ord,
        }
    }
}

impl PartialOrd for VAFRange {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone, Default, Derefable)]
pub(crate) struct VAFUniverse(#[deref(mutable)] HashSet<VAFSpectrum>);

impl VAFUniverse {
    pub(crate) fn contains(&self, vaf: AlleleFreq) -> bool {
        for atom in &**self {
            if atom.contains(vaf) {
                return true;
            }
        }
        false
    }
}

impl<'de> Deserialize<'de> for VAFUniverse {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de>,
    {
        deserializer.deserialize_string(VAFUniverseVisitor)
    }
}

struct VAFUniverseVisitor;

impl<'de> de::Visitor<'de> for VAFUniverseVisitor {
    type Value = VAFUniverse;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str(
            "a disjunction of possible VAFs (see https://varlociraptor.github.io/docs/calling)",
        )
    }

    fn visit_str<E>(self, v: &str) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        let res = FormulaParser::parse(Rule::universe, v);
        match res {
            Ok(pairs) => {
                let mut operands = HashSet::new();
                for pair in pairs {
                    match pair.as_rule() {
                        Rule::vaf => {
                            operands.insert(parse_vaf(pair));
                        }
                        Rule::vafrange => {
                            let inner = pair.into_inner();
                            operands.insert(parse_vafrange(inner));
                        }
                        Rule::EOI => (),
                        _ => unreachable!(),
                    }
                }
                Ok(VAFUniverse(operands))
            }
            Err(e) => {
                eprintln!("{}", e);
                Err(de::Error::invalid_value(
                    serde::de::Unexpected::Other("invalid VAF formula"),
                    &self,
                ))
            }
        }
    }
}

impl<'de> Deserialize<'de> for Formula {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de>,
    {
        deserializer.deserialize_string(FormulaVisitor)
    }
}

struct FormulaVisitor;

impl<'de> de::Visitor<'de> for FormulaVisitor {
    type Value = Formula;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter
            .write_str("a valid VAF formula (see https://varlociraptor.github.io/docs/calling)")
    }

    fn visit_str<E>(self, v: &str) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        let res = FormulaParser::parse(Rule::formula, v);
        match res {
            Ok(mut pairs) => {
                let pair = pairs.next().expect("bug: expecting formula");
                parse_formula(pair)
            }
            Err(e) => Err(de::Error::invalid_value(
                serde::de::Unexpected::Other(&format!("invalid VAF formula:\n{}", e)),
                &self,
            )),
        }
    }
}

fn parse_vaf(pair: Pair<Rule>) -> VAFSpectrum {
    let vaf = pair.as_str().parse().expect("bug: unable to parse VAF");
    VAFSpectrum::singleton(AlleleFreq(vaf))
}

fn parse_vafrange(mut inner: Pairs<Rule>) -> VAFSpectrum {
    let left = inner.next().unwrap().as_str();
    let lower = inner.next().unwrap().as_str().parse().unwrap();
    let upper = inner.next().unwrap().as_str().parse().unwrap();
    let right = inner.next().unwrap().as_str();

    let range = VAFRange {
        inner: lower..upper,
        left_exclusive: left == "]",
        right_exclusive: right == "[",
    };

    VAFSpectrum::Range(range)
}

fn parse_formula<E>(pair: Pair<Rule>) -> Result<Formula, E>
where
    E: de::Error,
{
    Ok(match pair.as_rule() {
        Rule::expression => {
            let mut inner = pair.into_inner();
            let identifier = inner.next().unwrap().as_str();
            Formula::Terminal(FormulaTerminal::Expression {
                identifier: ExpressionIdentifier(identifier.to_owned()),
                negated: false,
            })
        }
        Rule::variant => {
            let mut inner = pair.into_inner();
            let refbase = inner.next().unwrap().as_str().as_bytes()[0];
            let altbase = inner.next().unwrap().as_str().as_bytes()[0];
            Formula::Terminal(FormulaTerminal::Variant {
                refbase: Iupac(refbase),
                altbase: Iupac(altbase),
                positive: true,
            })
        }
        Rule::sample_vaf => {
            let mut inner = pair.into_inner();
            let sample = inner.next().unwrap().as_str().to_owned();
            Formula::Terminal(FormulaTerminal::Atom {
                sample,
                vafs: parse_vaf(inner.next().unwrap()),
            })
        }
        Rule::sample_vafrange => {
            let mut inner = pair.into_inner();
            let sample = inner.next().unwrap().as_str().to_owned();
            Formula::Terminal(FormulaTerminal::Atom {
                sample,
                vafs: parse_vafrange(inner.next().unwrap().into_inner()),
            })
        }
        Rule::conjunction => {
            let inner = pair.into_inner();
            let mut operands = Vec::new();
            for operand in inner {
                operands.push(parse_formula(operand)?);
            }
            Formula::Conjunction { operands }
        }
        Rule::disjunction => {
            let inner = pair.into_inner();
            let mut operands = Vec::new();
            for operand in inner {
                operands.push(parse_formula(operand)?);
            }
            Formula::Disjunction { operands }
        }
        Rule::negation => {
            let mut inner = pair.into_inner();
            Formula::Negation {
                operand: Box::new(parse_formula(inner.next().unwrap())?),
            }
        }
        Rule::formula => unreachable!(),
        Rule::subformula => unreachable!(),
        Rule::vafdef => unreachable!(),
        Rule::bound => unreachable!(),
        Rule::universe => unreachable!(),
        Rule::vafrange => unreachable!(),
        Rule::identifier => unreachable!(),
        Rule::vaf => unreachable!(),
        Rule::sample_vafdef => unreachable!(),
        Rule::EOI => unreachable!(),
        Rule::WHITESPACE => unreachable!(),
        Rule::COMMENT => unreachable!(),
        Rule::iupac => unreachable!(),
    })
}

#[cfg(test)]
mod test {
    use crate::grammar::Scenario;
    use crate::grammar::{Formula, VAFRange};
    use crate::variants::model::AlleleFreq;

    #[test]
    fn test_vaf_range_overlap() {
        let r1 = VAFRange {
            inner: AlleleFreq(0.0)..AlleleFreq(0.7),
            left_exclusive: false,
            right_exclusive: true,
        };
        let r2 = VAFRange {
            inner: AlleleFreq(0.3)..AlleleFreq(1.0),
            left_exclusive: false,
            right_exclusive: false,
        };
        let expected = VAFRange {
            inner: AlleleFreq(0.3)..AlleleFreq(0.7),
            left_exclusive: false,
            right_exclusive: true,
        };
        assert_eq!(expected, r1 & r2);
    }

    #[test]
    fn test_range_conjunction() {
        let scenario: Scenario = serde_yaml::from_str(
            r#"samples:
  normal:
    resolution: 100
    universe: "[0.0,1.0]"
events:
  full: "normal:[0.0,1.0]"
  part1: "normal:[0.0,0.7]"
  part2: "normal:[0.3,1.0]"
  expected: "normal:[0.3,0.7]""#,
        )
        .unwrap();
        let expected = scenario.events["expected"].clone();
        let full = scenario.events["full"].clone();
        let part1 = scenario.events["part1"].clone();
        let part2 = scenario.events["part2"].clone();
        let conjunction = Formula::Conjunction {
            operands: vec![part1, part2],
        }
        .normalize(&scenario, "all")
        .unwrap();
        assert_eq!(conjunction, expected.normalize(&scenario, "all").unwrap());
        assert_ne!(conjunction, full.normalize(&scenario, "all").unwrap());
    }

    #[test]
    fn test_nested_range_disjunction() {
        let scenario: Scenario = serde_yaml::from_str(
            r#"samples:
  normal:
    resolution: 100
    universe: "[0.0,1.0]"
events:
  full: "(normal:[0.0, 0.25] | normal:[0.5,0.75]) | (normal:[0.25,0.5] | normal:[0.75,1.0]) | normal:[0.1,0.4] | normal:0.1"
  expected: "normal:[0.0,1.0]""#,
        )
            .unwrap();
        let expected = scenario.events["expected"].clone();
        let full = scenario.events["full"].clone();
        let full = full.normalize(&scenario, "all").unwrap();
        assert_eq!(full, expected.normalize(&scenario, "all").unwrap());
    }

    #[test]
    fn test_two_separate_range_disjunctions() {
        let scenario: Scenario = serde_yaml::from_str(
            r#"samples:
  normal:
    resolution: 100
    universe: "[0.0,1.0]"
events:
  full: "(normal:[0.0, 0.25] | normal:[0.5,0.6]) | ((normal:[0.25,0.5] | normal:[0.7,0.9]) | normal:[0.9,1.0]) | normal:]0.8,0.9[ | normal:0.75"
  expected: "normal:[0.0,0.6] | normal:[0.7,1.0]""#,
        )
            .unwrap();
        let expected = scenario.events["expected"].clone();
        let full = scenario.events["full"].clone();
        let full = full.normalize(&scenario, "all").unwrap();
        assert_eq!(full, expected.normalize(&scenario, "all").unwrap());
    }
}
