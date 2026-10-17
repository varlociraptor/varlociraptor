// boolean_expression: a Rust crate for Boolean expressions and BDDs.
//
// Copyright (c) 2016 Chris Fallin <cfallin@c1f.net>. Released under the MIT
// License.
//

use smallvec::SmallVec;
use std::cmp::{Eq, Ord, Ordering, PartialEq, PartialOrd};
use std::collections::VecDeque;
use std::iter;
use std::slice;

/// A variable assignment in a cube.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum CubeVar {
    /// This variable must be false.
    False,
    /// This variable must be true.
    True,
    /// This variable may be true or false.
    DontCare,
}

const CUBE_ALLOCED_SIZE: usize = 16;

/// A `Cube` is one (multidimensional) cube in the variable space: i.e., a
/// particular set of variable assignments, where each variable is assigned
/// either true, false, or "don't care".
#[derive(Clone, Debug)]
pub struct Cube(SmallVec<[CubeVar; CUBE_ALLOCED_SIZE]>);

/// The result of attempting to merge two cubes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CubeMergeResult {
    /// The cubes could not be merged.
    None,
    /// The left cube was canceled because it is completely covered by the right cube.
    CancelLeft,
    /// The right cube was canceled because it is completely covered by the left cube.
    CancelRight,
    /// The two cubes merge into one.
    Merge(Cube),
    /// The left cube may be expanded (increase its number of `DontCare`s) by
    /// overlapping with the right cube.
    ExpandLeft(Cube),
    /// The right cube may be expanded (increase its number of `DontCare`s) by
    /// overlapping with the left cube.
    ExpandRight(Cube),
}

impl Cube {
    /// Construct an always-true cube (all variables are `DontCare`) for `vars`
    /// variables.
    pub fn true_cube(vars: usize) -> Cube {
        Cube(iter::repeat(CubeVar::DontCare).take(vars).collect())
    }

    /// Return an iterator over variable assignments.
    pub fn vars(&self) -> slice::Iter<CubeVar> {
        self.0.iter()
    }

    /// Attempt to merge this cube with another, which may cancel one or the
    /// other (if completely covered), expand one or the other (if possible, to
    /// increase the number of `DontCare`s and thus simplify the final
    /// expression), or merge the two into a single cube, or do nothing.
    pub fn merge_with(&self, other: &Cube) -> CubeMergeResult {
        // Cubes that differ in exactly one position can merge.
        // A cube that matches another except has fixed values where the other
        // has don't-cares can be eliminated.
        if self.0.len() != other.0.len() {
            CubeMergeResult::None
        } else if self == other {
            CubeMergeResult::CancelRight
        } else {
            let mut mismatches = 0;
            let mut mismatch_pos = 0;
            let mut left_covered = 0;
            let mut right_covered = 0;
            for (i, (lvar, rvar)) in self.0.iter().zip(other.0.iter()).enumerate() {
                match (lvar, rvar) {
                    (&CubeVar::False, &CubeVar::True) | (&CubeVar::True, &CubeVar::False) => {
                        mismatches += 1;
                        mismatch_pos = i;
                    }
                    (&CubeVar::False, &CubeVar::DontCare)
                    | (&CubeVar::True, &CubeVar::DontCare) => {
                        left_covered += 1;
                    }
                    (&CubeVar::DontCare, &CubeVar::False)
                    | (&CubeVar::DontCare, &CubeVar::True) => {
                        right_covered += 1;
                    }
                    _ => {}
                }
            }
            if mismatches == 0 && left_covered > 0 && right_covered == 0 {
                CubeMergeResult::CancelLeft
            } else if mismatches == 0 && right_covered > 0 && left_covered == 0 {
                CubeMergeResult::CancelRight
            } else if mismatches == 1 && right_covered == 0 && left_covered == 0 {
                CubeMergeResult::Merge(self.with_var(mismatch_pos, CubeVar::DontCare))
            } else if mismatches == 1 && right_covered > 0 && left_covered == 0 {
                CubeMergeResult::ExpandRight(other.with_var(mismatch_pos, CubeVar::DontCare))
            } else if mismatches == 1 && right_covered == 0 && left_covered > 0 {
                CubeMergeResult::ExpandLeft(self.with_var(mismatch_pos, CubeVar::DontCare))
            } else {
                CubeMergeResult::None
            }
        }
    }

    /// Return a new cube equal to this cube, but with the particular variable
    // assigned the given value.
    pub fn with_var(&self, idx: usize, val: CubeVar) -> Cube {
        Cube(
            self.0
                .iter()
                .enumerate()
                .map(|(i, var)| if i == idx { val.clone() } else { var.clone() })
                .collect(),
        )
    }
}

impl PartialEq for Cube {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for Cube {}
impl PartialOrd for Cube {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Cube {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.iter().cmp(other.0.iter())
    }
}

const CUBELIST_ALLOCED_SIZE: usize = 4;

/// A `CubeList` is a sum (OR'd list) of cubes. It represents a Boolean
/// expression in sum-of-products form, by construction.
///
/// The `CubeList` abstraction supports only indexed anonymous variables
/// (variable 0, 1, ...), and does not (yet) have a wrapper supporting an
/// arbitrary terminal type `T`. This may be implemented in the future.
///
/// The `CubeList` abstraction is used internally to convert from a `BDD`
/// to a quasi-minimized Boolean expression.
#[derive(Clone, Debug)]
pub struct CubeList(SmallVec<[Cube; CUBE_ALLOCED_SIZE]>);

impl PartialEq for CubeList {
    fn eq(&self, other: &Self) -> bool {
        self.0.iter().cmp(other.0.iter()) == Ordering::Equal
    }
}

impl CubeList {
    /// Construct a new, empty, cube list (equivalent to a constant `false`).
    pub fn new() -> CubeList {
        CubeList(SmallVec::new())
    }

    /// Construct a cube list from a list of `Cube` structs.
    pub fn from_list(cubes: &[Cube]) -> CubeList {
        CubeList(cubes.iter().map(|c| c.clone()).collect())
    }

    /// Return an iterator over all cubes.
    pub fn cubes(&self) -> slice::Iter<Cube> {
        self.0.iter()
    }

    /// Merge this cube list with another, canceling or merging cubes where
    /// possible. The resulting cube list is not guaranteed to be minimal (that
    /// is the set-cover problem, which is NP-Complete), but is reduced somewhat
    /// by a very simple reduction/merging algorithm.
    pub fn merge(&self, other: &CubeList) -> CubeList {
        let mut out: SmallVec<[Cube; CUBE_ALLOCED_SIZE]> = SmallVec::new();
        let mut canceled: SmallVec<[bool; CUBE_ALLOCED_SIZE]> = SmallVec::new();
        for cube in self.0.iter().chain(other.0.iter()) {
            out.push(cube.clone());
            canceled.push(false);
        }

        let mut worklist = VecDeque::new();
        for i in 0..out.len() {
            worklist.push_back(i);
        }
        while let Some(i) = worklist.pop_front() {
            if canceled[i] {
                continue;
            }
            for j in 0..out.len() {
                if i == j {
                    continue;
                }
                if canceled[i] {
                    break;
                }
                if canceled[j] {
                    continue;
                }
                match out[i].merge_with(&out[j]) {
                    CubeMergeResult::None => {}
                    CubeMergeResult::CancelLeft => {
                        canceled[i] = true;
                    }
                    CubeMergeResult::CancelRight => {
                        canceled[j] = true;
                    }
                    CubeMergeResult::Merge(n) => {
                        out[i] = n;
                        worklist.push_back(i);
                        canceled[j] = true;
                    }
                    CubeMergeResult::ExpandLeft(n) => {
                        out[i] = n;
                        worklist.push_back(i);
                    }
                    CubeMergeResult::ExpandRight(n) => {
                        out[j] = n;
                        worklist.push_back(j);
                    }
                }
            }
        }

        let out = out
            .into_iter()
            .zip(canceled.iter())
            .filter(|&(_, flag)| !flag)
            .map(|(cube, _)| cube)
            .collect();
        CubeList(out)
    }

    pub fn with_var(&self, idx: usize, val: CubeVar) -> CubeList {
        CubeList(
            self.0
                .iter()
                .map(|c| c.with_var(idx, val.clone()))
                .collect(),
        )
    }
}

mod test {
    use super::*;

    fn make_cube(elems: &[u32]) -> Cube {
        Cube(
            elems
                .iter()
                .map(|i| match *i {
                    0 => CubeVar::False,
                    1 => CubeVar::True,
                    _ => CubeVar::DontCare,
                })
                .collect(),
        )
    }

    #[test]
    fn cube_eq() {
        assert!(make_cube(&[1, 0, 0]) == make_cube(&[1, 0, 0]));
        assert!(make_cube(&[1, 0, 0]) != make_cube(&[1, 0, 1]));
    }

    #[test]
    fn cube_ord() {
        assert!(make_cube(&[1, 0, 0]) < make_cube(&[1, 1, 0]));
        assert!(make_cube(&[1, 0, 2]) > make_cube(&[1, 0, 1]));
    }

    #[test]
    fn cube_with_var() {
        assert!(make_cube(&[0, 1, 0]).with_var(2, CubeVar::True) == make_cube(&[0, 1, 1]));
    }

    #[test]
    fn cube_merge() {
        // Cubes of mismatching dimension: no cancelation.
        assert!(make_cube(&[0, 0]).merge_with(&make_cube(&[0])) == CubeMergeResult::None);
        // Identical cubes: cancelation (our implementation: cancel right).
        assert!(make_cube(&[0, 0]).merge_with(&make_cube(&[0, 0])) == CubeMergeResult::CancelRight);
        // Irredundant cubes: no cancelation.
        assert!(make_cube(&[1, 0]).merge_with(&make_cube(&[0, 1])) == CubeMergeResult::None);
        // Irredundant cubes with some overlap: no cancelation.
        assert!(make_cube(&[1, 2]).merge_with(&make_cube(&[2, 1])) == CubeMergeResult::None);
        // Left cube covers right cube: cancel right.
        assert!(
            make_cube(&[1, 2, 2]).merge_with(&make_cube(&[1, 1, 0]))
                == CubeMergeResult::CancelRight
        );
        // Right cube may be expanded to overlap with left cube.
        assert!(
            make_cube(&[1, 1, 2]).merge_with(&make_cube(&[1, 0, 0]))
                == CubeMergeResult::ExpandRight(make_cube(&[1, 2, 0]))
        );
        // The above, with right and left flipped.
        assert!(
            make_cube(&[1, 1, 0]).merge_with(&make_cube(&[1, 2, 2])) == CubeMergeResult::CancelLeft
        );
        assert!(
            make_cube(&[1, 0, 0]).merge_with(&make_cube(&[1, 1, 2]))
                == CubeMergeResult::ExpandLeft(make_cube(&[1, 2, 0]))
        );
        // Cubes with one mismatch: merge.
        assert!(
            make_cube(&[1, 1, 0]).merge_with(&make_cube(&[1, 1, 1]))
                == CubeMergeResult::Merge(make_cube(&[1, 1, 2]))
        );
        // Cubes with more than one mismatch: no merge.
        assert!(make_cube(&[1, 1, 0]).merge_with(&make_cube(&[1, 0, 1])) == CubeMergeResult::None);
    }

    #[test]
    fn cubelist_merge() {
        // No merges.
        assert!(
            CubeList::new().merge(&CubeList::from_list(&[
                make_cube(&[1, 0, 0]),
                make_cube(&[0, 1, 1])
            ])) == CubeList::from_list(&[make_cube(&[1, 0, 0]), make_cube(&[0, 1, 1])])
        );
        // Multiple-stage merge.
        assert!(
            CubeList::new().merge(&CubeList::from_list(&[
                make_cube(&[1, 0, 0]),
                make_cube(&[1, 1, 1]),
                make_cube(&[1, 0, 1]),
                make_cube(&[1, 1, 0])
            ])) == CubeList::from_list(&[make_cube(&[1, 2, 2])])
        );
        // Last term merges with first.
        assert!(
            CubeList::new().merge(&CubeList::from_list(&[
                make_cube(&[1, 0, 0]),
                make_cube(&[0, 1, 1]),
                make_cube(&[1, 0, 1])
            ])) == CubeList::from_list(&[make_cube(&[1, 0, 2]), make_cube(&[0, 1, 1])])
        );
        // New item cancels existing in list.
        assert!(
            CubeList::new().merge(&CubeList::from_list(&[
                make_cube(&[1, 0, 0]),
                make_cube(&[1, 2, 2])
            ])) == CubeList::from_list(&[make_cube(&[1, 2, 2])])
        );
        // Existing list item cancels new item.
        assert!(
            CubeList::new().merge(&CubeList::from_list(&[
                make_cube(&[1, 2, 2]),
                make_cube(&[1, 0, 0])
            ])) == CubeList::from_list(&[make_cube(&[1, 2, 2])])
        );
    }
}
