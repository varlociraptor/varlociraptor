use crate::Counter;

use std::hash::Hash;

impl<'a, T, N> IntoIterator for &'a Counter<T, N>
where
    T: Hash + Eq,
{
    type Item = (&'a T, &'a N);
    type IntoIter = std::collections::hash_map::Iter<'a, T, N>;

    fn into_iter(self) -> Self::IntoIter {
        self.map.iter()
    }
}

impl<T, N> IntoIterator for Counter<T, N>
where
    T: Hash + Eq,
{
    type Item = (T, N);
    type IntoIter = std::collections::hash_map::IntoIter<T, N>;

    /// Consumes the `Counter` to produce an iterator that owns the values it returns.
    ///
    /// # Examples
    /// ```rust
    /// # use counter::Counter;
    ///
    /// let counter: Counter<_> = "aaab".chars().collect();
    ///
    /// let vec: Vec<_> = counter.into_iter().collect();
    ///
    /// for (item, count) in &vec {
    ///     if item == &'a' {
    ///         assert_eq!(count, &3);
    ///     }
    ///     if item == &'b' {
    ///         assert_eq!(count, &1);
    ///     }
    /// }
    /// ```

    fn into_iter(self) -> Self::IntoIter {
        self.map.into_iter()
    }
}

impl<'a, T, N> IntoIterator for &'a mut Counter<T, N>
where
    T: Hash + Eq,
{
    type Item = (&'a T, &'a mut N);
    type IntoIter = std::collections::hash_map::IterMut<'a, T, N>;

    /// Creates an iterator that provides mutable references to the counts, but keeps the keys immutable.
    ///
    /// # Examples
    /// ```rust
    /// # use counter::Counter;
    ///
    /// let mut counter: Counter<_> = "aaab".chars().collect();
    ///
    /// for (item, count) in &mut counter {
    ///     if *item == 'a' {
    ///         // 'a' is so great it counts as 2
    ///         *count *= 2;
    ///     }
    /// }
    ///
    /// assert_eq!(counter[&'a'], 6);
    /// assert_eq!(counter[&'b'], 1);
    /// ```
    fn into_iter(self) -> Self::IntoIter {
        self.map.iter_mut()
    }
}
