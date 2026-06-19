use crate::Counter;

use num_traits::Zero;

use std::hash::Hash;
use std::ops::{Sub, SubAssign};

impl<T, N> Sub for Counter<T, N>
where
    T: Hash + Eq,
    N: PartialOrd + PartialEq + SubAssign + Zero,
{
    type Output = Counter<T, N>;

    /// Subtract (keeping only positive values).
    ///
    /// `out = c - d;` -> `out[x] == c[x] - d[x]` for all `x`,
    /// keeping only items with a value greater than [`N::zero()`].
    ///
    /// [`N::zero()`]:
    /// https://docs.rs/num-traits/latest/num_traits/identities/trait.Zero.html#tymethod.zero
    ///
    /// ```rust
    /// # use counter::Counter;
    /// # use std::collections::HashMap;
    /// let c = "aaab".chars().collect::<Counter<_>>();
    /// let d = "abb".chars().collect::<Counter<_>>();
    ///
    /// let e = c - d;
    ///
    /// let expect = [('a', 2)].iter().cloned().collect::<HashMap<_, _>>();
    /// assert_eq!(e.into_map(), expect);
    /// ```
    fn sub(mut self, rhs: Counter<T, N>) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<T, N> SubAssign for Counter<T, N>
where
    T: Hash + Eq,
    N: PartialOrd + PartialEq + SubAssign + Zero,
{
    /// Subtract (keeping only positive values).
    ///
    /// `c -= d;` -> `c[x] -= d[x]` for all `x`,
    /// keeping only items with a value greater than [`N::zero()`].
    ///
    /// [`N::zero()`]:
    /// https://docs.rs/num-traits/latest/num_traits/identities/trait.Zero.html#tymethod.zero
    ///
    /// ```rust
    /// # use counter::Counter;
    /// # use std::collections::HashMap;
    /// let mut c = "aaab".chars().collect::<Counter<_>>();
    /// let d = "abb".chars().collect::<Counter<_>>();
    ///
    /// c -= d;
    ///
    /// let expect = [('a', 2)].iter().cloned().collect::<HashMap<_, _>>();
    /// assert_eq!(c.into_map(), expect);
    /// ```
    fn sub_assign(&mut self, rhs: Self) {
        for (key, value) in rhs.map {
            let mut remove = false;
            if let Some(entry) = self.map.get_mut(&key) {
                if *entry >= value {
                    *entry -= value;
                } else {
                    remove = true;
                }
                if *entry == N::zero() {
                    remove = true;
                }
            }
            if remove {
                self.map.remove(&key);
            }
        }
    }
}
