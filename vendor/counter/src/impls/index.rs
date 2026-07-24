use crate::Counter;

use num_traits::Zero;

use std::borrow::Borrow;
use std::hash::Hash;
use std::ops::{Index, IndexMut};

impl<T, Q, N> Index<&'_ Q> for Counter<T, N>
where
    T: Hash + Eq + Borrow<Q>,
    Q: Hash + Eq,
    N: Zero,
{
    type Output = N;

    /// Index in immutable contexts.
    ///
    /// Returns a reference to a [`zero`] value for missing keys.
    ///
    /// [`zero`]:
    /// https://docs.rs/num-traits/latest/num_traits/identities/trait.Zero.html#tymethod.zero
    ///
    /// ```
    /// # use counter::Counter;
    /// let counter = "aabbcc".chars().collect::<Counter<_>>();
    /// assert_eq!(counter[&'a'], 2);
    /// assert_eq!(counter[&'b'], 2);
    /// assert_eq!(counter[&'c'], 2);
    /// assert_eq!(counter[&'d'], 0);
    /// ```
    ///
    /// Note that the [`zero`] is a struct field but not one of the values of the inner
    /// [`HashMap`].  This method does not modify any existing value.
    ///
    /// [`zero`]:
    /// https://docs.rs/num-traits/latest/num_traits/identities/trait.Zero.html#tymethod.zero
    /// [`HashMap`]: https://doc.rust-lang.org/stable/std/collections/struct.HashMap.html
    ///
    /// ```
    /// # use counter::Counter;
    /// let counter = "".chars().collect::<Counter<_>>();
    /// assert_eq!(counter[&'a'], 0);
    /// assert_eq!(counter.get(&'a'), None); // as `Deref<Target = HashMap<_, _>>`
    /// ```
    fn index(&self, key: &'_ Q) -> &N {
        self.map.get(key).unwrap_or(&self.zero)
    }
}

impl<T, Q, N> IndexMut<&'_ Q> for Counter<T, N>
where
    T: Hash + Eq + Borrow<Q>,
    Q: Hash + Eq + ToOwned<Owned = T>,
    N: Zero,
{
    /// Index in mutable contexts.
    ///
    /// If the given key is not present, creates a new entry and initializes it with a [`zero`]
    /// value.
    ///
    /// [`zero`]:
    /// https://docs.rs/num-traits/latest/num_traits/identities/trait.Zero.html#tymethod.zero
    ///
    /// ```
    /// # use counter::Counter;
    /// let mut counter = "aabbcc".chars().collect::<Counter<_>>();
    /// counter[&'c'] += 1;
    /// counter[&'d'] += 1;
    /// assert_eq!(counter[&'c'], 3);
    /// assert_eq!(counter[&'d'], 1);
    /// ```
    ///
    /// Unlike `Index::index`, the returned mutable reference to the [`zero`] is actually one of the
    /// values of the inner [`HashMap`].
    ///
    /// [`zero`]:
    /// https://docs.rs/num-traits/latest/num_traits/identities/trait.Zero.html#tymethod.zero
    /// [`HashMap`]: https://doc.rust-lang.org/stable/std/collections/struct.HashMap.html
    ///
    /// ```
    /// # use counter::Counter;
    /// let mut counter = "".chars().collect::<Counter<_>>();
    /// assert_eq!(counter.get(&'a'), None); // as `Deref<Target = HashMap<_, _>>`
    /// let _ = &mut counter[&'a'];
    /// assert_eq!(counter.get(&'a'), Some(&0));
    /// ```
    fn index_mut(&mut self, key: &'_ Q) -> &mut N {
        self.map.entry(key.to_owned()).or_insert_with(N::zero)
    }
}
