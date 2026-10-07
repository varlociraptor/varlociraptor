use crate::Counter;

use num_traits::{One, Zero};

use std::hash::Hash;
use std::ops::{Add, AddAssign};

impl<I, T, N> Add<I> for Counter<T, N>
where
    I: IntoIterator<Item = T>,
    T: Hash + Eq,
    N: AddAssign + Zero + One,
{
    type Output = Self;
    /// Consume `self` producing a `Counter` like `self` updated with the counts of
    /// the elements of `I`.
    ///
    /// ```rust
    /// # use counter::Counter;
    /// # use std::collections::HashMap;
    /// let counter = "abbccc".chars().collect::<Counter<_>>();
    ///
    /// let new_counter = counter + "aeeeee".chars();
    /// let expected: HashMap<char, usize> = [('a', 2), ('b', 2), ('c', 3), ('e', 5)]
    ///     .iter().cloned().collect();
    /// assert_eq!(new_counter.into_map(), expected);
    /// ```
    fn add(mut self, rhs: I) -> Self::Output {
        self.update(rhs);
        self
    }
}

impl<I, T, N> AddAssign<I> for Counter<T, N>
where
    I: IntoIterator<Item = T>,
    T: Hash + Eq,
    N: AddAssign + Zero + One,
{
    /// Directly add the counts of the elements of `I` to `self`.
    ///
    /// ```rust
    /// # use counter::Counter;
    /// # use std::collections::HashMap;
    /// let mut counter = "abbccc".chars().collect::<Counter<_>>();
    ///
    /// counter += "aeeeee".chars();
    /// let expected: HashMap<char, usize> = [('a', 2), ('b', 2), ('c', 3), ('e', 5)]
    ///     .iter().cloned().collect();
    /// assert_eq!(counter.into_map(), expected);
    /// ```
    fn add_assign(&mut self, rhs: I) {
        self.update(rhs);
    }
}
