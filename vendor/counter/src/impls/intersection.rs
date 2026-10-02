use crate::Counter;

use num_traits::Zero;

use std::hash::Hash;
use std::ops::{BitAnd, BitAndAssign};

impl<T, N> BitAnd for Counter<T, N>
where
    T: Hash + Eq,
    N: Ord + Zero,
{
    type Output = Counter<T, N>;

    /// Returns the intersection of `self` and `rhs` as a new `Counter`.
    ///
    /// `out = c & d;` -> `out[x] == min(c[x], d[x])`
    ///
    /// ```rust
    /// # use counter::Counter;
    /// # use std::collections::HashMap;
    /// let c = "aaab".chars().collect::<Counter<_>>();
    /// let d = "abb".chars().collect::<Counter<_>>();
    ///
    /// let e = c & d;
    ///
    /// let expect = [('a', 1), ('b', 1)].iter().cloned().collect::<HashMap<_, _>>();
    /// assert_eq!(e.into_map(), expect);
    /// ```
    fn bitand(self, mut rhs: Counter<T, N>) -> Self::Output {
        use std::cmp::min;

        let mut counter = Counter::new();
        for (key, lhs_count) in self.map {
            if let Some(rhs_count) = rhs.remove(&key) {
                let count = min(lhs_count, rhs_count);
                counter.map.insert(key, count);
            }
        }
        counter
    }
}

impl<T, N> BitAndAssign for Counter<T, N>
where
    T: Hash + Eq,
    N: Ord + Zero,
{
    /// Updates `self` with the intersection of `self` and `rhs`
    ///
    /// `c &= d;` -> `c[x] == min(c[x], d[x])`
    ///
    /// ```rust
    /// # use counter::Counter;
    /// # use std::collections::HashMap;
    /// let mut c = "aaab".chars().collect::<Counter<_>>();
    /// let d = "abb".chars().collect::<Counter<_>>();
    ///
    /// c &= d;
    ///
    /// let expect = [('a', 1), ('b', 1)].iter().cloned().collect::<HashMap<_, _>>();
    /// assert_eq!(c.into_map(), expect);
    /// ```
    fn bitand_assign(&mut self, mut rhs: Counter<T, N>) {
        for (key, rhs_count) in rhs.drain() {
            if rhs_count < self[&key] {
                self.map.insert(key, rhs_count);
            }
        }
    }
}
