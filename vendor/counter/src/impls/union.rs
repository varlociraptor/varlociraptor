use crate::Counter;

use num_traits::Zero;

use std::hash::Hash;
use std::ops::{BitOr, BitOrAssign};

impl<T, N> BitOr for Counter<T, N>
where
    T: Hash + Eq,
    N: Ord + Zero,
{
    type Output = Counter<T, N>;

    /// Returns the union of `self` and `rhs` as a new `Counter`.
    ///
    /// `out = c | d;` -> `out[x] == max(c[x], d[x])`
    ///
    /// ```rust
    /// # use counter::Counter;
    /// # use std::collections::HashMap;
    /// let c = "aaab".chars().collect::<Counter<_>>();
    /// let d = "abb".chars().collect::<Counter<_>>();
    ///
    /// let e = c | d;
    ///
    /// let expect = [('a', 3), ('b', 2)].iter().cloned().collect::<HashMap<_, _>>();
    /// assert_eq!(e.into_map(), expect);
    /// ```
    fn bitor(mut self, rhs: Counter<T, N>) -> Self::Output {
        for (key, rhs_value) in rhs.map {
            let entry = self.map.entry(key).or_insert_with(N::zero);
            // We want to update the value of the now occupied entry in `self` with the maximum of
            // its current value and `rhs_value`.  If that max is `rhs_value`, we can just update
            // the value of the entry.  If the max is the current value, we do nothing.  Note that
            // `Ord::max()` returns the second argument (here `rhs_value`) if its two arguments are
            // equal, justifying the use of the weak inequality below instead of a strict
            // inequality.
            //
            // Doing it this way with an inequality instead of actually using `std::cmp::max()`
            // lets us avoid trying (and failing) to move the non-copy value out of the entry in
            // order to pass it as an argument to `std::cmp::max()`, while still holding a mutable
            // reference to the value slot in the entry.
            //
            // And while using the inequality seemingly only requires the bound `N: PartialOrd`, we
            // nevertheless prefer to require `Ord` as though we were using `std::cmp::max()`
            // because the semantics of `BitOr` for `Counter` really do not make sense if there are
            // possibly non-comparable values of type `N`.
            if rhs_value >= *entry {
                *entry = rhs_value;
            }
        }
        self
    }
}

impl<T, N> BitOrAssign for Counter<T, N>
where
    T: Hash + Eq,
    N: Ord + Zero,
{
    /// Updates `self` with the union of `self` and `rhs`
    ///
    /// `c |= d;` -> `c[x] == max(c[x], d[x])`
    ///
    /// ```rust
    /// # use counter::Counter;
    /// # use std::collections::HashMap;
    /// let mut c = "aaab".chars().collect::<Counter<_>>();
    /// let d = "abb".chars().collect::<Counter<_>>();
    ///
    /// c |= d;
    ///
    /// let expect = [('a', 3), ('b', 2)].iter().cloned().collect::<HashMap<_, _>>();
    /// assert_eq!(c.into_map(), expect);
    /// ```
    fn bitor_assign(&mut self, mut rhs: Counter<T, N>) {
        for (key, rhs_count) in rhs.drain() {
            if rhs_count > self[&key] {
                self.map.insert(key, rhs_count);
            }
        }
    }
}
