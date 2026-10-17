// boolean_expression: a Rust crate for Boolean expressions and BDDs.
//
// Copyright (c) 2016 Chris Fallin <cfallin@c1f.net>. Released under the MIT
// License.
//

use itertools::Itertools;
use std::cmp;
use std::collections::hash_map::Entry as HashEntry;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fmt;
use std::fmt::Debug;
use std::hash::Hash;
use std::usize;

use cubes::{Cube, CubeList, CubeVar};
use idd::*;
use Expr;

/// A `BDDFunc` is a function index within a particular `BDD`. It must only
/// be used with the `BDD` instance which produced it.
pub type BDDFunc = usize;

/// A special terminal `BDDFunc` which is constant `false` (zero).
pub const BDD_ZERO: BDDFunc = usize::MAX;
/// A special terminal `BDDFunc` which is constant `true` (one).
pub const BDD_ONE: BDDFunc = usize::MAX - 1;

pub(crate) type BDDLabel = usize;

#[derive(Clone, PartialEq, Eq, Hash)]
pub(crate) struct BDDNode {
    pub label: BDDLabel,
    pub lo: BDDFunc,
    pub hi: BDDFunc,
    pub varcount: usize,
}

fn bdd_func_str(b: BDDFunc) -> String {
    if b == BDD_ZERO {
        "ZERO".to_owned()
    } else if b == BDD_ONE {
        "ONE".to_owned()
    } else {
        format!("{}", b)
    }
}

impl fmt::Debug for BDDNode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "BDDNode(label = {}, lo = {}, hi = {})",
            self.label,
            bdd_func_str(self.lo),
            bdd_func_str(self.hi)
        )
    }
}

#[derive(Clone, Debug)]
pub(crate) struct LabelBDD {
    pub nodes: Vec<BDDNode>,
    dedup_hash: HashMap<BDDNode, BDDFunc>,
}

impl LabelBDD {
    pub fn new() -> LabelBDD {
        LabelBDD {
            nodes: Vec::new(),
            dedup_hash: HashMap::new(),
        }
    }

    fn get_node(&mut self, label: BDDLabel, lo: BDDFunc, hi: BDDFunc) -> BDDFunc {
        if lo == hi {
            return lo;
        }
        let n = BDDNode {
            label: label,
            lo: lo,
            hi: hi,
            varcount: cmp::min(self.sat_varcount(lo), self.sat_varcount(hi) + 1),
        };
        match self.dedup_hash.entry(n.clone()) {
            HashEntry::Occupied(o) => *o.get(),
            HashEntry::Vacant(v) => {
                let idx = self.nodes.len() as BDDFunc;
                self.nodes.push(n);
                v.insert(idx);
                idx
            }
        }
    }

    fn sat_varcount(&self, f: BDDFunc) -> usize {
        if f == BDD_ZERO || f == BDD_ONE {
            0
        } else {
            self.nodes[f].varcount
        }
    }

    pub fn terminal(&mut self, label: BDDLabel) -> BDDFunc {
        self.get_node(label, BDD_ZERO, BDD_ONE)
    }

    pub fn constant(&mut self, value: bool) -> BDDFunc {
        if value {
            BDD_ONE
        } else {
            BDD_ZERO
        }
    }

    /// Restrict: fundamental building block of logical combinators. Takes a
    /// Shannon cofactor: i.e., returns a new function based on `f` but with the
    /// given label forced to the given value.
    pub fn restrict(&mut self, f: BDDFunc, label: BDDLabel, val: bool) -> BDDFunc {
        if f == BDD_ZERO {
            return BDD_ZERO;
        }
        if f == BDD_ONE {
            return BDD_ONE;
        }

        let node = self.nodes[f].clone();
        if label < node.label {
            f
        } else if label == node.label {
            if val {
                node.hi
            } else {
                node.lo
            }
        } else {
            let lo = self.restrict(node.lo, label, val);
            let hi = self.restrict(node.hi, label, val);
            self.get_node(node.label, lo, hi)
        }
    }

    fn min_label(&self, f: BDDFunc) -> Option<BDDLabel> {
        if f == BDD_ZERO || f == BDD_ONE {
            None
        } else {
            Some(self.nodes[f].label)
        }
    }

    /// If-then-else: fundamental building block of logical combinators. Works
    /// by divide-and-conquer: split on the lowest appearing label, take Shannon
    /// cofactors for the two cases, recurse, and recombine with a new node.
    pub fn ite(&mut self, i: BDDFunc, t: BDDFunc, e: BDDFunc) -> BDDFunc {
        if i == BDD_ONE {
            t
        } else if i == BDD_ZERO {
            e
        } else if t == e {
            t
        } else if t == BDD_ONE && e == BDD_ZERO {
            i
        } else {
            let i_var = self.min_label(i).unwrap_or(usize::MAX);
            let t_var = self.min_label(t).unwrap_or(usize::MAX);
            let e_var = self.min_label(e).unwrap_or(usize::MAX);
            let split = cmp::min(i_var, cmp::min(t_var, e_var));
            assert!(split != usize::MAX);
            let i_lo = self.restrict(i, split, false);
            let t_lo = self.restrict(t, split, false);
            let e_lo = self.restrict(e, split, false);
            let i_hi = self.restrict(i, split, true);
            let t_hi = self.restrict(t, split, true);
            let e_hi = self.restrict(e, split, true);
            let lo = self.ite(i_lo, t_lo, e_lo);
            let hi = self.ite(i_hi, t_hi, e_hi);
            self.get_node(split, lo, hi)
        }
    }

    pub fn not(&mut self, n: BDDFunc) -> BDDFunc {
        self.ite(n, BDD_ZERO, BDD_ONE)
    }

    pub fn and(&mut self, a: BDDFunc, b: BDDFunc) -> BDDFunc {
        self.ite(a, b, BDD_ZERO)
    }

    pub fn or(&mut self, a: BDDFunc, b: BDDFunc) -> BDDFunc {
        self.ite(a, BDD_ONE, b)
    }

    pub fn xor(&mut self, a: BDDFunc, b: BDDFunc) -> BDDFunc {
        let not_b = self.not(b);
        self.ite(a, not_b, b)
    }

    pub fn implies(&mut self, a: BDDFunc, b: BDDFunc) -> BDDFunc {
        let not_a = self.not(a);
        self.or(not_a, b)
    }

    pub fn evaluate(&self, func: BDDFunc, inputs: &[bool]) -> Option<bool> {
        let mut f = func;
        for (i, val) in inputs.iter().enumerate() {
            if f == BDD_ZERO || f == BDD_ONE {
                break;
            }
            let node = &self.nodes[f];
            if node.label > i {
                continue;
            } else if node.label == i {
                f = if *val { node.hi } else { node.lo };
            }
        }
        match f {
            BDD_ZERO => Some(false),
            BDD_ONE => Some(true),
            _ => None,
        }
    }

    fn compute_cubelist(&self, memoize_vec: &mut Vec<Option<CubeList>>, n: BDDFunc, nvars: usize) {
        if memoize_vec[n].is_some() {
            return;
        }
        let label = self.nodes[n].label;
        let lo = self.nodes[n].lo;
        let hi = self.nodes[n].hi;
        let lo_list = match lo {
            BDD_ZERO => CubeList::new(),
            BDD_ONE => CubeList::from_list(&[Cube::true_cube(nvars)])
                .with_var(label as usize, CubeVar::False),
            _ => {
                self.compute_cubelist(memoize_vec, lo, nvars);
                memoize_vec[lo]
                    .as_ref()
                    .unwrap()
                    .with_var(label as usize, CubeVar::False)
            }
        };
        let hi_list = match hi {
            BDD_ZERO => CubeList::new(),
            BDD_ONE => CubeList::from_list(&[Cube::true_cube(nvars)])
                .with_var(label as usize, CubeVar::True),
            _ => {
                self.compute_cubelist(memoize_vec, hi, nvars);
                memoize_vec[hi]
                    .as_ref()
                    .unwrap()
                    .with_var(label as usize, CubeVar::True)
            }
        };
        let new_list = lo_list.merge(&hi_list);
        memoize_vec[n] = Some(new_list);
    }

    fn cube_to_expr(&self, c: &Cube) -> Expr<BDDLabel> {
        c.vars()
            .enumerate()
            .flat_map(|(i, v)| match v {
                &CubeVar::False => Some(Expr::not(Expr::Terminal(i))),
                &CubeVar::True => Some(Expr::Terminal(i)),
                &CubeVar::DontCare => None,
            })
            .fold1(|a, b| Expr::and(a, b))
            .unwrap_or(Expr::Const(true))
    }

    fn cubelist_to_expr(&self, c: &CubeList) -> Expr<BDDLabel> {
        c.cubes()
            .map(|c| self.cube_to_expr(c))
            .fold1(|a, b| Expr::or(a, b))
            .unwrap_or(Expr::Const(false))
    }

    pub fn to_expr(&self, func: BDDFunc, nvars: usize) -> Expr<BDDLabel> {
        if func == BDD_ZERO {
            Expr::Const(false)
        } else if func == BDD_ONE {
            Expr::Const(true)
        } else {
            // At each node, we construct a cubelist, starting from the roots.
            let mut cubelists: Vec<Option<CubeList>> = Vec::with_capacity(self.nodes.len());
            cubelists.resize(self.nodes.len(), None);
            self.compute_cubelist(&mut cubelists, func, nvars);
            self.cubelist_to_expr(cubelists[func].as_ref().unwrap())
        }
    }

    /// Returns a function that is true whenever the maximum number of
    /// functions in `funcs` are true.
    pub fn max_sat(&mut self, funcs: &[BDDFunc]) -> BDDFunc {
        // First, construct an IDD function for each BDD function,
        // with value 1 if true and 0 if false. Then add these
        // together to obtain a single IDD function whose value is the
        // number of satisfied (true) BDD functions.
        let mut idd = LabelIDD::from_bdd(self);
        let idd_funcs: Vec<_> = funcs.iter().map(|f| idd.from_bdd_func(*f)).collect();
        let satisfied_count = idd_funcs
            .iter()
            .cloned()
            .fold1(|a, b| idd.add(a.clone(), b.clone()))
            .unwrap();

        // Now, find the maximum reachable count.
        let max_count = idd.max_value(satisfied_count.clone());

        // Finally, return a boolean function that is true when the
        // maximal number of functions are satisfied.
        let c = idd.constant(max_count);
        idd.eq(satisfied_count, c, self)
    }
}

/// A `BDD` is a Binary Decision Diagram, an efficient way to represent a
/// Boolean function in a canonical way. (It is actually a "Reduced Ordered
/// Binary Decision Diagram", which gives it its canonicity assuming terminals
/// are ordered consistently.)
///
/// A BDD is built up from terminals (free variables) and constants, combined
/// with the logical combinators AND, OR, and NOT. It may be evaluated with
/// certain terminal assignments.
///
/// The major advantage of a BDD is that its logical operations are performed,
/// it will "self-simplify": i.e., taking the OR of `And(a, b)` and `And(a,
/// Not(b))` will produce `a` without any further simplification step. Furthermore,
/// the `BDDFunc` representing this value is canonical: if two different
/// expressions are produced within the same BDD and they both result in
/// (simplify down to) `a`, then the `BDDFunc` values will be equal. The
/// tradeoff is that logical operations may be expensive: they are linear in
/// BDD size, but BDDs may have exponential size (relative to terminal count)
/// in the worst case.
#[derive(Clone, Debug)]
pub struct BDD<T>
where
    T: Clone + Debug + Eq + Hash,
{
    bdd: LabelBDD,
    labels: HashMap<T, BDDLabel>,
    rev_labels: Vec<T>,
}

impl<T> BDD<T>
where
    T: Clone + Debug + Eq + Hash,
{
    /// Produce a new, empty, BDD.
    pub fn new() -> BDD<T> {
        BDD {
            bdd: LabelBDD::new(),
            labels: HashMap::new(),
            rev_labels: Vec::new(),
        }
    }

    fn label(&mut self, t: T) -> BDDLabel {
        match self.labels.entry(t.clone()) {
            HashEntry::Occupied(o) => *o.get(),
            HashEntry::Vacant(v) => {
                let next_id = self.rev_labels.len() as BDDLabel;
                v.insert(next_id);
                self.rev_labels.push(t);
                next_id
            }
        }
    }

    /// Produce a function within the BDD representing the terminal `t`. If
    /// this terminal has been used in the BDD before, the same `BDDFunc` will be
    /// returned.
    pub fn terminal(&mut self, t: T) -> BDDFunc {
        let l = self.label(t);
        self.bdd.terminal(l)
    }

    /// Produce a function within the BDD representing the constant value `val`.
    pub fn constant(&mut self, val: bool) -> BDDFunc {
        self.bdd.constant(val)
    }

    /// Produce a function within the BDD representing the logical if-then-else
    /// of the functions `i`, `t`, and `e`
    pub fn ite(&mut self, i: BDDFunc, t: BDDFunc, e: BDDFunc) -> BDDFunc {
        self.bdd.ite(i, t, e)
    }

    /// Produce a function within the BDD representing the logical complement
    /// of the function `n`.
    pub fn not(&mut self, n: BDDFunc) -> BDDFunc {
        self.bdd.not(n)
    }

    /// Produce a function within the BDD representing the logical AND of the
    /// functions `a` and `b`.
    pub fn and(&mut self, a: BDDFunc, b: BDDFunc) -> BDDFunc {
        self.bdd.and(a, b)
    }

    /// Produce a function within the BDD representing the logical OR of the
    /// functions `a` and `b`.
    pub fn or(&mut self, a: BDDFunc, b: BDDFunc) -> BDDFunc {
        self.bdd.or(a, b)
    }

    /// Produce a function within the BDD representing the logical XOR of the
    /// functions `a` and `b`.
    pub fn xor(&mut self, a: BDDFunc, b: BDDFunc) -> BDDFunc {
        self.bdd.xor(a, b)
    }

    /// Produce a function within the BDD representing the logical implication `a` -> `b`.
    pub fn implies(&mut self, a: BDDFunc, b: BDDFunc) -> BDDFunc {
        self.bdd.implies(a, b)
    }

    /// Check whether the function `f` within the BDD is satisfiable.
    pub fn sat(&self, f: BDDFunc) -> bool {
        match f {
            BDD_ZERO => false,
            _ => true,
        }
    }

    /// Return a new function based on `f` but with the given label forced to the given value.
    pub fn restrict(&mut self, f: BDDFunc, t: T, val: bool) -> BDDFunc {
        self.bdd.restrict(f, self.labels[&t], val)
    }

    /// Produce a function within the BDD representing the given expression
    /// `e`, which may contain ANDs, ORs, NOTs, terminals, and constants.
    pub fn from_expr(&mut self, e: &Expr<T>) -> BDDFunc {
        match e {
            &Expr::Terminal(ref t) => self.terminal(t.clone()),
            &Expr::Const(val) => self.constant(val),
            &Expr::Not(ref x) => {
                let xval = self.from_expr(&**x);
                self.not(xval)
            }
            &Expr::And(ref a, ref b) => {
                let aval = self.from_expr(&**a);
                let bval = self.from_expr(&**b);
                self.and(aval, bval)
            }
            &Expr::Or(ref a, ref b) => {
                let aval = self.from_expr(&**a);
                let bval = self.from_expr(&**b);
                self.or(aval, bval)
            }
        }
    }

    /// Evaluate the function `f` in the BDD with the given terminal
    /// assignments. Any terminals not specified in `values` default to `false`.
    pub fn evaluate(&self, f: BDDFunc, values: &HashMap<T, bool>) -> bool {
        let size = self.rev_labels.len();
        let mut valarray = Vec::with_capacity(size);
        valarray.resize(size, false);
        for (t, l) in &self.labels {
            valarray[*l as usize] = *values.get(t).unwrap_or(&false);
        }
        self.bdd.evaluate(f, &valarray).unwrap()
    }

    /// Compute an assignment for terminals which satisfies 'f'.  If
    /// satisfiable, this function returns a HashMap with the
    /// assignments (true, false) for terminals unless a terminal's
    /// assignment does not matter for satisfiability. If 'f' is not
    /// satisfiable, returns None.
    ///
    /// Example: for the boolean function "a or b", this function
    /// could return one of the following two HashMaps: {"a" -> true}
    /// or {"b" -> true}.
    pub fn sat_one(&self, f: BDDFunc) -> Option<HashMap<T, bool>> {
        let mut h = HashMap::new();
        if self.sat_one_internal(f, &mut h) {
            Some(h)
        } else {
            None
        }
    }

    fn sat_one_internal(&self, f: BDDFunc, assignments: &mut HashMap<T, bool>) -> bool {
        match f {
            BDD_ZERO => false,
            BDD_ONE => true,
            _ => {
                let hi = self.bdd.nodes[f].hi;
                let lo = self.bdd.nodes[f].lo;
                if hi != BDD_ZERO
                    && (lo == BDD_ZERO || self.bdd.sat_varcount(hi) < self.bdd.sat_varcount(lo))
                {
                    assignments.insert(self.rev_labels[self.bdd.nodes[f].label].clone(), true);
                    self.sat_one_internal(hi, assignments);
                } else {
                    assignments.insert(self.rev_labels[self.bdd.nodes[f].label].clone(), false);
                    self.sat_one_internal(lo, assignments);
                }
                true
            }
        }
    }

    /// Convert the BDD to a minimized sum-of-products expression.
    pub fn to_expr(&self, f: BDDFunc) -> Expr<T> {
        self.bdd
            .to_expr(f, self.rev_labels.len())
            .map(|t: &BDDLabel| self.rev_labels[*t as usize].clone())
    }

    /// Export BDD to `dot` format (from the graphviz package) to enable visualization.
    pub fn to_dot(&self, f: BDDFunc) -> String {
        // the algorithm starts at the f BDDfunction and then recursively collects all BDDNodes
        // until BDD_ZERO and BDD_ONE. The output for each node is straightforward: just a single
        // `dot` node.
        let mut out = String::from("digraph bdd {\n");

        let mut nodes = BTreeSet::new();
        self.reachable_nodes(f, &mut nodes);
        for func in nodes {
            if func <= f {
                out.push_str(&format!(
                    "n{} [label = {:?}];\n",
                    func, self.rev_labels[self.bdd.nodes[func].label]
                ));
                out.push_str(&format!(
                    "n{} -> n{} [style=dotted];\n",
                    func, self.bdd.nodes[func].lo
                ));
                out.push_str(&format!("n{} -> n{};\n", func, self.bdd.nodes[func].hi));
            }
        }

        out.push_str(&format!("n{} [label=\"0\"];\n", BDD_ZERO));
        out.push_str(&format!("n{} [label=\"1\"];\n", BDD_ONE));

        out.push_str("}");

        out.to_string()
    }

    /// Produce an ordered set of nodes in the BDD function `f`: the transitive closure of
    /// reachable nodes.
    fn reachable_nodes(&self, f: BDDFunc, s: &mut BTreeSet<BDDFunc>) {
        if f != BDD_ZERO && f != BDD_ONE {
            // we use a BTreeSet instead of a HashSet since its order is stable.
            if s.insert(f) {
                self.reachable_nodes(self.bdd.nodes[f].hi, s);
                self.reachable_nodes(self.bdd.nodes[f].lo, s);
            }
        }
    }

    /// Produce a function that is true when the maximal number of
    /// given input functions are true.
    pub fn max_sat(&mut self, funcs: &[BDDFunc]) -> BDDFunc {
        self.bdd.max_sat(funcs)
    }

    /// Return a vector of all labels in the BDD.
    pub fn labels(&self) -> Vec<T> {
        self.labels.keys().cloned().collect()
    }
}

/// The `BDDOutput` trait provides an interface to inform a listener about new
/// BDD nodes that are created. It allows the user to persist a BDD to a stream
/// (e.g., a log or trace file) as a long-running process executes. A
/// `BDDOutput` instance may be provided to all BDD operations.
pub trait BDDOutput<T, E> {
    fn write_label(&self, label: T, label_id: u64) -> Result<(), E>;
    fn write_node(
        &self,
        node_id: BDDFunc,
        label_id: u64,
        lo: BDDFunc,
        hi: BDDFunc,
    ) -> Result<(), E>;
}

/// A `PersistedBDD` is a wrapper around a `BDD` that provides a means to write
/// BDD labels and nodes out to a `BDDOutput`. It tracks how much of the BDD
/// has already been writen out, and writes out new nodes and labels as
/// required when its `persist()` or `persist_all()` method is called.
pub struct PersistedBDD<T>
where
    T: Clone + Debug + Eq + Hash,
{
    bdd: BDD<T>,
    next_output_func: BDDFunc,
    next_output_label: BDDLabel,
}

impl<T> PersistedBDD<T>
where
    T: Clone + Debug + Eq + Hash,
{
    /// Create a new `PersistedBDD`.
    pub fn new() -> PersistedBDD<T> {
        PersistedBDD {
            bdd: BDD::new(),
            next_output_func: 0,
            next_output_label: 0,
        }
    }

    /// Return the inner BDD.
    pub fn bdd(&self) -> &BDD<T> {
        &self.bdd
    }

    /// Return the inner BDD.
    pub fn bdd_mut(&mut self) -> &mut BDD<T> {
        &mut self.bdd
    }

    /// Persist (at least) all labels and nodes in the BDD necessary to fully
    /// describe BDD function `f`. More records than strictly necessary may be
    /// written out.
    pub fn persist<E>(&mut self, f: BDDFunc, out: &dyn BDDOutput<T, E>) -> Result<(), E> {
        if f == BDD_ZERO || f == BDD_ONE {
            // No need to persist the terminal constant nodes!
            return Ok(());
        }
        while self.next_output_label < self.bdd.rev_labels.len() {
            let id = self.next_output_label;
            let t = self.bdd.rev_labels[id].clone();
            out.write_label(t, id as u64)?;
            self.next_output_label += 1;
        }
        while self.next_output_func <= f {
            let id = self.next_output_func;
            let node = &self.bdd.bdd.nodes[id];
            out.write_node(id, node.label as u64, node.lo, node.hi)?;
            self.next_output_func += 1;
        }
        Ok(())
    }

    /// Persist all labels and nodes in the BDD.
    pub fn persist_all<E>(&mut self, out: &dyn BDDOutput<T, E>) -> Result<(), E> {
        if self.bdd.bdd.nodes.len() > 0 {
            let last_f = self.bdd.bdd.nodes.len() - 1;
            self.persist(last_f, out)
        } else {
            Ok(())
        }
    }
}

/// A `BDDLoader` provides a way to inject BDD nodes directly, as they were
/// previously dumped by a `PersistedBDD` to a `BDDOutput`. The user should
/// create a `BDDLoader` instance wrapped around a `BDD` and call
/// `inject_label` and `inject_node` as appropriate to inject labels and nodes.
pub struct BDDLoader<'a, T>
where
    T: Clone + Debug + Eq + Hash + 'a,
{
    bdd: &'a mut BDD<T>,
}

impl<'a, T> BDDLoader<'a, T>
where
    T: Clone + Debug + Eq + Hash + 'a,
{
    /// Create a new `BDDLoader` wrapping the given `bdd`. The `BDDLoader`
    /// holds a mutable reference to `bdd` until destroyed. `bdd` must be empty
    /// initially.
    pub fn new(bdd: &'a mut BDD<T>) -> BDDLoader<'a, T> {
        assert!(bdd.labels.len() == 0);
        assert!(bdd.rev_labels.len() == 0);
        assert!(bdd.bdd.nodes.len() == 0);
        BDDLoader { bdd: bdd }
    }

    /// Inject a new label into the BDD. The `id` must be the next consecutive
    /// `id`; i.e., labels must be injected in the order they were dumped to a
    /// `BDDOutput`.
    pub fn inject_label(&mut self, t: T, id: u64) {
        assert!(id == self.bdd.rev_labels.len() as u64);
        self.bdd.rev_labels.push(t.clone());
        self.bdd.labels.insert(t, id as BDDLabel);
    }

    /// Inject a new node into the BDD. The `id` must be the next consecutive
    /// `id`; i.e., nodes must be injected in the order they were dumped to a
    /// `BDDOutput`.
    pub fn inject_node(&mut self, id: BDDFunc, label_id: u64, lo: BDDFunc, hi: BDDFunc) {
        assert!(id == self.bdd.bdd.nodes.len() as BDDFunc);
        let n = BDDNode {
            label: label_id as BDDLabel,
            lo: lo,
            hi: hi,
            varcount: cmp::min(
                self.bdd.bdd.sat_varcount(lo),
                self.bdd.bdd.sat_varcount(hi) + 1,
            ),
        };
        self.bdd.bdd.nodes.push(n.clone());
        self.bdd.bdd.dedup_hash.insert(n, id);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::cell::RefCell;
    use std::collections::HashMap;
    use Expr;
    extern crate rand;
    extern crate rand_xorshift;
    use self::rand::{Rng, SeedableRng};
    use self::rand_xorshift::XorShiftRng;

    fn term_hashmap(vals: &[bool], h: &mut HashMap<u32, bool>) {
        h.clear();
        for (i, v) in vals.iter().enumerate() {
            h.insert(i as u32, *v);
        }
    }

    fn test_bdd(
        b: &BDD<u32>,
        f: BDDFunc,
        h: &mut HashMap<u32, bool>,
        inputs: &[bool],
        expected: bool,
    ) {
        term_hashmap(inputs, h);
        assert!(b.evaluate(f, h) == expected);
    }

    #[test]
    fn bdd_eval() {
        let mut h = HashMap::new();
        let mut b = BDD::new();
        let expr = Expr::or(
            Expr::and(Expr::Terminal(0), Expr::Terminal(1)),
            Expr::and(Expr::not(Expr::Terminal(2)), Expr::not(Expr::Terminal(3))),
        );
        let f = b.from_expr(&expr);
        test_bdd(&b, f, &mut h, &[false, false, true, true], false);
        test_bdd(&b, f, &mut h, &[true, false, true, true], false);
        test_bdd(&b, f, &mut h, &[true, true, true, true], true);
        test_bdd(&b, f, &mut h, &[false, false, false, true], false);
        test_bdd(&b, f, &mut h, &[false, false, false, false], true);
    }

    fn bits_to_hashmap(bits: usize, n: usize, h: &mut HashMap<u32, bool>) {
        for b in 0..bits {
            h.insert(b as u32, (n & (1 << b)) != 0);
        }
    }

    fn test_bdd_expr(e: Expr<u32>, nterminals: usize) {
        let mut b = BDD::new();
        let f = b.from_expr(&e);
        let mut terminal_values = HashMap::new();
        let mut expected_satisfiable = false;
        for v in 0..(1 << nterminals) {
            bits_to_hashmap(nterminals, v, &mut terminal_values);
            let expr_val = e.evaluate(&terminal_values);
            let bdd_val = b.evaluate(f, &terminal_values);
            assert!(expr_val == bdd_val);
            if expr_val {
                expected_satisfiable = true;
            }
        }
        // test sat_one
        let sat_result = b.sat_one(f);
        assert!(sat_result.is_some() == expected_satisfiable);
        if expected_satisfiable {
            assert!(b.evaluate(f, &sat_result.unwrap()));
        }
    }

    fn random_expr<R: Rng>(r: &mut R, nterminals: usize) -> Expr<u32> {
        match r.gen_range(0, 5) {
            0 => Expr::Terminal(r.gen_range(0, nterminals) as u32),
            1 => Expr::Const(r.gen::<bool>()),
            2 => Expr::Not(Box::new(random_expr(r, nterminals))),
            3 => Expr::And(
                Box::new(random_expr(r, nterminals)),
                Box::new(random_expr(r, nterminals)),
            ),
            4 => Expr::Or(
                Box::new(random_expr(r, nterminals)),
                Box::new(random_expr(r, nterminals)),
            ),
            _ => unreachable!(),
        }
    }

    #[test]
    fn bdd_exhaustive_exprs() {
        let mut rng = XorShiftRng::from_seed(rand::thread_rng().gen());
        for _ in 0..100 {
            let expr = random_expr(&mut rng, 6);
            test_bdd_expr(expr, 6);
        }
    }

    #[test]
    fn bdd_to_expr() {
        let mut b = BDD::new();
        let f_true = b.constant(true);
        assert!(b.to_expr(f_true) == Expr::Const(true));
        let f_false = b.constant(false);
        assert!(b.to_expr(f_false) == Expr::Const(false));
        let f_0 = b.terminal(0);
        let f_1 = b.terminal(1);
        let f_and = b.and(f_0, f_1);
        assert!(b.to_expr(f_and) == Expr::and(Expr::Terminal(0), Expr::Terminal(1)));
        let f_or = b.or(f_0, f_1);
        assert!(b.to_expr(f_or) == Expr::or(Expr::Terminal(1), Expr::Terminal(0)));
        let f_not = b.not(f_0);
        assert!(b.to_expr(f_not) == Expr::not(Expr::Terminal(0)));
        let f_2 = b.terminal(2);
        let f_1_or_2 = b.or(f_1, f_2);
        let f_0_and_1_or_2 = b.and(f_0, f_1_or_2);
        assert!(
            b.to_expr(f_0_and_1_or_2)
                == Expr::or(
                    Expr::and(Expr::Terminal(0), Expr::Terminal(2)),
                    Expr::and(Expr::Terminal(0), Expr::Terminal(1))
                )
        );
    }

    #[derive(Clone, Debug)]
    struct InMemoryBDDLog {
        labels: RefCell<Vec<(u64, String)>>,
        nodes: RefCell<Vec<(BDDFunc, u64, BDDFunc, BDDFunc)>>,
    }

    impl InMemoryBDDLog {
        pub fn new() -> InMemoryBDDLog {
            InMemoryBDDLog {
                labels: RefCell::new(Vec::new()),
                nodes: RefCell::new(Vec::new()),
            }
        }
    }

    impl BDDOutput<String, ()> for InMemoryBDDLog {
        fn write_label(&self, l: String, label_id: u64) -> Result<(), ()> {
            let mut labels = self.labels.borrow_mut();
            labels.push((label_id, l));
            Ok(())
        }

        fn write_node(
            &self,
            node_id: BDDFunc,
            label_id: u64,
            lo: BDDFunc,
            hi: BDDFunc,
        ) -> Result<(), ()> {
            let mut nodes = self.nodes.borrow_mut();
            nodes.push((node_id, label_id, lo, hi));
            Ok(())
        }
    }

    #[test]
    // the tests compare the dot output to a string which has been manually verified to be correct
    fn dot_output() {
        let mut bdd = BDD::new();
        let a = bdd.terminal("a");
        let b = bdd.terminal("b");
        let b_and_a = bdd.and(a, b);
        {
            let dot = bdd.to_dot(b_and_a);
            assert_eq!(
                dot,
                indoc!(
                    "
                            digraph bdd {
                            n1 [label = \"b\"];
                            n1 -> n18446744073709551615 [style=dotted];
                            n1 -> n18446744073709551614;
                            n2 [label = \"a\"];
                            n2 -> n18446744073709551615 [style=dotted];
                            n2 -> n1;
                            n18446744073709551615 [label=\"0\"];
                            n18446744073709551614 [label=\"1\"];
                            }"
                )
            );
        }
        let c = bdd.terminal("c");
        let c_or_a = bdd.or(c, a);
        {
            let dot = bdd.to_dot(c_or_a);
            assert_eq!(
                dot,
                indoc!(
                    "
                            digraph bdd {
                            n3 [label = \"c\"];
                            n3 -> n18446744073709551615 [style=dotted];
                            n3 -> n18446744073709551614;
                            n4 [label = \"a\"];
                            n4 -> n3 [style=dotted];
                            n4 -> n18446744073709551614;
                            n18446744073709551615 [label=\"0\"];
                            n18446744073709551614 [label=\"1\"];
                            }"
                )
            );
        }
        let not_c_and_b = bdd.not(b_and_a);
        let c_or_a_and_not_b_and_a = bdd.and(c_or_a, not_c_and_b);
        {
            let dot = bdd.to_dot(c_or_a_and_not_b_and_a);
            assert_eq!(
                dot,
                indoc!(
                    "
                            digraph bdd {
                            n3 [label = \"c\"];
                            n3 -> n18446744073709551615 [style=dotted];
                            n3 -> n18446744073709551614;
                            n5 [label = \"b\"];
                            n5 -> n18446744073709551614 [style=dotted];
                            n5 -> n18446744073709551615;
                            n7 [label = \"a\"];
                            n7 -> n3 [style=dotted];
                            n7 -> n5;
                            n18446744073709551615 [label=\"0\"];
                            n18446744073709551614 [label=\"1\"];
                            }"
                )
            );
        }
        {
            let new_a = bdd.terminal("a");
            let d = bdd.terminal("d");
            let new_a_or_d = bdd.or(new_a, d);
            let dot = bdd.to_dot(new_a_or_d);
            assert_eq!(
                dot,
                indoc!(
                    "
                            digraph bdd {
                            n8 [label = \"d\"];
                            n8 -> n18446744073709551615 [style=dotted];
                            n8 -> n18446744073709551614;
                            n9 [label = \"a\"];
                            n9 -> n8 [style=dotted];
                            n9 -> n18446744073709551614;
                            n18446744073709551615 [label=\"0\"];
                            n18446744073709551614 [label=\"1\"];
                            }"
                )
            );
        }
    }

    #[test]
    fn sat_one() {
        let mut bdd = BDD::new();

        // empty bdds
        assert!(bdd.sat_one(BDD_ONE).is_some());
        assert!(bdd.sat_one(BDD_ZERO).is_none());

        let a = bdd.terminal("a");
        let b = bdd.terminal("b");
        let b_and_a = bdd.and(a, b);
        let result = bdd.sat_one(b_and_a);
        assert!(result.is_some());
        assert!(bdd.evaluate(b_and_a, &result.unwrap()));

        let c = bdd.terminal("c");
        let not_c = bdd.not(c);
        let b_and_a_or_not_c = bdd.or(b_and_a, not_c);
        let result = bdd.sat_one(b_and_a_or_not_c);
        assert!(result.is_some());
        assert!(bdd.evaluate(b_and_a_or_not_c, &result.unwrap()));

        // unsatisfiable formula
        let c_and_not_c = bdd.and(c, not_c);
        assert!(bdd.sat_one(c_and_not_c).is_none());
    }

    #[test]
    fn max_sat() {
        let mut bdd = BDD::new();
        // Test: a, a+b, a+c, c', c, bd, ad, d'
        let a = bdd.terminal(0);
        let b = bdd.terminal(1);
        let c = bdd.terminal(2);
        let d = bdd.terminal(3);
        let cnot = bdd.not(c);
        let dnot = bdd.not(d);
        let ab = bdd.and(a, b);
        let ac = bdd.and(a, c);
        let bd = bdd.and(b, d);
        let ad = bdd.and(a, d);
        let max_sat = bdd.max_sat(&[a, ab, ac, cnot, c, bd, ad, dnot]);
        let abc = bdd.and(ab, c);
        let abcd = bdd.and(abc, d);
        assert!(max_sat == abcd);
    }

    #[test]
    fn max_sat_min_var() {
        let mut bdd = BDD::new();
        // Test: a, a+b, a+c, c', c, bd, d'
        let a = bdd.terminal(0);
        let b = bdd.terminal(1);
        let c = bdd.terminal(2);
        let d = bdd.terminal(3);
        let cnot = bdd.not(c);
        let dnot = bdd.not(d);
        let ab = bdd.and(a, b);
        let ac = bdd.and(a, c);
        let bd = bdd.and(b, d);
        let max_sat = bdd.max_sat(&[a, ab, ac, cnot, c, bd, dnot]);
        let abc = bdd.and(ab, c);
        assert!(max_sat == abc);
        let assn = bdd.sat_one(max_sat).unwrap();
        assert!(assn.get(&0) == Some(&true));
        assert!(assn.get(&1) == Some(&true));
        assert!(assn.get(&2) == Some(&true));
        assert!(assn.get(&3) == None);
    }

    #[test]
    fn persist_bdd() {
        let out = InMemoryBDDLog::new();
        let mut p = PersistedBDD::new();
        let term_a = p.bdd_mut().terminal("A".to_owned());
        let term_b = p.bdd_mut().terminal("B".to_owned());
        let term_c = p.bdd_mut().terminal("C".to_owned());
        let ab = p.bdd_mut().and(term_a, term_b);
        let ab_or_c = p.bdd_mut().or(ab, term_c);
        p.persist(ab_or_c, &out).unwrap();
        assert!(
            *out.labels.borrow()
                == vec![
                    (0, "A".to_owned()),
                    (1, "B".to_owned()),
                    (2, "C".to_owned()),
                ]
        );
        assert!(
            *out.nodes.borrow()
                == vec![
                    (0, 0, BDD_ZERO, BDD_ONE),
                    (1, 1, BDD_ZERO, BDD_ONE),
                    (2, 2, BDD_ZERO, BDD_ONE),
                    (3, 0, BDD_ZERO, 1),
                    (4, 1, 2, BDD_ONE),
                    (5, 0, 2, 4),
                ]
        );
    }

    #[test]
    fn load_bdd() {
        let mut bdd = BDD::new();
        {
            let mut loader = BDDLoader::new(&mut bdd);
            loader.inject_label("A".to_owned(), 0);
            loader.inject_label("B".to_owned(), 1);
            loader.inject_node(0, 1, BDD_ZERO, BDD_ONE);
            loader.inject_node(1, 0, BDD_ZERO, 0);
        }
        let mut h = HashMap::new();
        h.insert("A".to_owned(), true);
        h.insert("B".to_owned(), true);
        assert!(bdd.evaluate(1, &h) == true);
    }

    #[test]
    fn persist_empty_bdd() {
        let out = InMemoryBDDLog::new();
        let mut p = PersistedBDD::new();
        p.persist(BDD_ZERO, &out).unwrap();
    }
}
