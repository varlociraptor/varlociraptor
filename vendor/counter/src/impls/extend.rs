use crate::Counter;

use num_traits::{One, Zero};

use std::hash::Hash;
use std::ops::AddAssign;

impl<T, N> Extend<T> for Counter<T, N>
where
    T: Hash + Eq,
    N: AddAssign + Zero + One,
{
    /// Extend a `Counter` with an iterator of items.
    ///
    /// ```rust
    /// # use counter::Counter;
    /// # use std::collections::HashMap;
    /// let mut counter = "abbccc".chars().collect::<Counter<_>>();
    /// counter.extend("bccddd".chars());
    /// let expect = [('a', 1), ('b', 3), ('c', 5), ('d', 3)].iter().cloned().collect::<HashMap<_, _>>();
    /// assert_eq!(counter.into_map(), expect);
    /// ```
    fn extend<I: IntoIterator<Item = T>>(&mut self, iter: I) {
        self.update(iter);
    }
}

impl<T, N> Extend<(T, N)> for Counter<T, N>
where
    T: Hash + Eq,
    N: AddAssign + Zero,
{
    /// Extend a counter with `(item, count)` tuples.
    ///
    /// The counts of duplicate items are summed.
    /// ```rust
    /// # use counter::Counter;
    /// # use std::collections::HashMap;
    /// let mut counter = "abbccc".chars().collect::<Counter<_>>();
    /// counter.extend([('a', 1), ('b', 2), ('c', 3), ('a', 4)].iter().cloned());
    /// let expect = [('a', 6), ('b', 4), ('c', 6)].iter()
    ///     .cloned().collect::<HashMap<_, _>>();
    /// assert_eq!(counter.into_map(), expect);
    /// ```
    fn extend<I: IntoIterator<Item = (T, N)>>(&mut self, iter: I) {
        for (item, item_count) in iter {
            let entry = self.map.entry(item).or_insert_with(N::zero);
            *entry += item_count;
        }
    }
}

impl<'a, T: 'a, N: 'a> Extend<(&'a T, &'a N)> for Counter<T, N>
where
    T: Hash + Eq + Clone,
    N: AddAssign + Zero + Clone,
{
    /// Extend a counter with `(item, count)` tuples.
    ///
    /// You can extend a `Counter` with another `Counter`:
    /// ```rust
    /// # use counter::Counter;
    /// # use std::collections::HashMap;
    /// let mut counter = "abbccc".chars().collect::<Counter<_>>();
    /// let another = "bccddd".chars().collect::<Counter<_>>();
    /// counter.extend(&another);
    /// let expect = [('a', 1), ('b', 3), ('c', 5), ('d', 3)].iter()
    ///     .cloned().collect::<HashMap<_, _>>();
    /// assert_eq!(counter.into_map(), expect);
    /// ```
    fn extend<I: IntoIterator<Item = (&'a T, &'a N)>>(&mut self, iter: I) {
        for (item, item_count) in iter {
            let entry = self.map.entry(item.clone()).or_insert_with(N::zero);
            *entry += item_count.clone();
        }
    }
}
