//! Counter counts recurrent elements of iterables. It is based on [the Python
//! implementation](https://docs.python.org/3/library/collections.html#collections.Counter).
//!
//! The struct [`Counter`](struct.Counter.html) is the entry-point type for this module.
//!
//! # Math Underpinnings
//!
//! Mathematically, a `Counter` implements a hash-based version of a [multiset],
//! or bag. This is simply an extension of the notion of a set to the idea that
//! we care not only about whether an entity exists within the set, but the number
//! of occurrences within the set. Normal set operations such as intersection,
//! union, etc. are of course still supported.
//!
//! [multiset]: https://en.wikipedia.org/wiki/Set_(abstract_data_type)#Multiset
//!
//! # Examples
//!
//! ## Just count an iterable
//!
//! ```rust
//! use counter::Counter;
//! let char_counts = "barefoot".chars().collect::<Counter<_>>();
//! let counts_counts = char_counts.values().collect::<Counter<_>>();
//! ```
//!
//! ## Update a count
//!
//! ```rust
//! # use counter::Counter;
//! let mut counts = "aaa".chars().collect::<Counter<_>>();
//! counts[&'a'] += 1;
//! counts[&'b'] += 1;
//! ```
//!
//! ```rust
//! # use counter::Counter;
//! let mut counts = "able babble table babble rabble table able fable scrabble"
//!     .split_whitespace().collect::<Counter<_>>();
//! // add or subtract an iterable of the same type
//! counts += "cain and abel fable table cable".split_whitespace();
//! // or add or subtract from another Counter of the same type
//! let other_counts = "scrabble cabbie fable babble"
//!     .split_whitespace().collect::<Counter<_>>();
//! let difference = counts - other_counts;
//! ```
//!
//! ## Extend a `Counter` with another `Counter`:
//! ```rust
//! # use counter::Counter;
//! # use std::collections::HashMap;
//! let mut counter = "abbccc".chars().collect::<Counter<_>>();
//! let another = "bccddd".chars().collect::<Counter<_>>();
//! counter.extend(&another);
//! let expect = [('a', 1), ('b', 3), ('c', 5), ('d', 3)].iter()
//!     .cloned().collect::<HashMap<_, _>>();
//! assert_eq!(counter.into_map(), expect);
//! ```
//! ## Get items with keys
//!
//! ```rust
//! # use counter::Counter;
//! let counts = "aaa".chars().collect::<Counter<_>>();
//! assert_eq!(counts[&'a'], 3);
//! assert_eq!(counts[&'b'], 0);
//! ```
//!
//! ## Get the most common items
//!
//! [`most_common_ordered()`] uses the natural ordering of keys which are [`Ord`].
//!
//! [`most_common_ordered()`]: Counter::most_common_ordered
//! [`Ord`]: https://doc.rust-lang.org/stable/std/cmp/trait.Ord.html
//!
//! ```rust
//! # use counter::Counter;
//! let by_common = "eaddbbccc".chars().collect::<Counter<_>>().most_common_ordered();
//! let expected = vec![('c', 3), ('b', 2), ('d', 2), ('a', 1), ('e', 1)];
//! assert!(by_common == expected);
//! ```
//!
//! [`k_most_common_ordered()`] takes an argument `k` of type `usize` and returns the top `k` most
//! common items.  This is functionally equivalent to calling `most_common_ordered()` and then
//! truncating the result to length `k`.  However, if `k` is smaller than the length of the counter
//! then `k_most_common_ordered()` can be more efficient, often much more so.
//!
//! ```rust
//! # use counter::Counter;
//! let by_common = "eaddbbccc".chars().collect::<Counter<_>>().k_most_common_ordered(2);
//! let expected = vec![('c', 3), ('b', 2)];
//! assert!(by_common == expected);
//! ```
//!
//! [`k_most_common_ordered()`]: Counter::k_most_common_ordered
//! [`most_common_ordered()`]: Counter::most_common_ordered
//!
//! ## Get the most common items using your own ordering
//!
//! For example, here we break ties reverse alphabetically.
//!
//! ```rust
//! # use counter::Counter;
//! let counter = "eaddbbccc".chars().collect::<Counter<_>>();
//! let by_common = counter.most_common_tiebreaker(|&a, &b| b.cmp(&a));
//! let expected = vec![('c', 3), ('d', 2), ('b', 2), ('e', 1), ('a', 1)];
//! assert!(by_common == expected);
//! ```
//!
//! ## Test counters against another
//!
//! Counters are multi-sets and so can be sub- or supersets of each other.
//!
//! A counter is a _subset_ of another if for all its elements, the other
//! counter has an equal or higher count. Test for this with [`is_subset()`]:
//!
//! ```rust
//! # use counter::Counter;
//! let counter = "aaabb".chars().collect::<Counter<_>>();
//! let superset = "aaabbbc".chars().collect::<Counter<_>>();
//! let not_a_superset = "aaae".chars().collect::<Counter<_>>();
//! assert!(counter.is_subset(&superset));
//! assert!(!counter.is_subset(&not_a_superset));
//! ```
//!
//! Testing for a _superset_ is the inverse, [`is_superset()`] is true if the counter can contain another counter in its entirety:
//!
//! ```rust
//! # use counter::Counter;
//! let counter = "aaabbbc".chars().collect::<Counter<_>>();
//! let subset = "aabbb".chars().collect::<Counter<_>>();
//! let not_a_subset = "aaae".chars().collect::<Counter<_>>();
//! assert!(counter.is_superset(&subset));
//! assert!(!counter.is_superset(&not_a_subset));
//! ```
//!
//! These relationships continue to work when [using a _signed_ integer type for the counter][signed]: all values in the subset must be equal or lower to the values in the superset. Negative
//! values are interpreted as 'missing' those values, and the subset would need to miss those
//! same elements, or be short more, to still be a subset:
//!
//! ```rust
//! # use counter::Counter;
//! let mut subset = "aaabb".chars().collect::<Counter<_, i8>>();
//! subset.insert('e', -2);  // short 2 'e's
//! subset.insert('f', -1);  // and 1 'f'
//! let mut superset = "aaaabbb".chars().collect::<Counter<_, i8>>();
//! superset.insert('e', -1);  // short 1 'e'
//! assert!(subset.is_subset(&superset));
//! assert!(superset.is_superset(&subset));
//! ```
//!
//! [`is_subset()`]: Counter::is_subset
//! [`is_superset()`]: Counter::is_superset
//! [signed]: #use-your-own-type-for-the-count
//!
//! ## Counter intersection and union
//!
//! You can intersect two counters, giving you the minimal counts of their
//! combined elements using the [`&` bitwise and operator][BitAnd], and produce
//! their union with the maximum counts using [`|` bitwise or][BitOr]:
//!
//! ```rust
//! # use counter::Counter;
//! let a = "aaabb".chars().collect::<Counter<_>>();
//! let b = "aabbbbe".chars().collect::<Counter<_>>();
//!
//! let intersection = a & b;
//! let expected_intersection = "aabb".chars().collect::<Counter<_>>();
//! assert_eq!(intersection, expected_intersection);
//!
//! let c = "aaabb".chars().collect::<Counter<_>>();
//! let d = "aabbbbe".chars().collect::<Counter<_>>();
//!
//! let union = c | d;
//! let expected_union = "aaabbbbe".chars().collect::<Counter<_>>();
//! assert_eq!(union, expected_union)
//! ```
//!
//! The in-place [`&=`] and [`|=`] operations are also supported.
//!
//! [BitAnd]: https://doc.rust-lang.org/std/ops/trait.BitAnd.html
//! [BitOr]: https://doc.rust-lang.org/std/ops/trait.BitOr.html
//! [`&=`]: https://doc.rust-lang.org/std/ops/trait.BitAndAssign.html
//! [`|=`]: https://doc.rust-lang.org/std/ops/trait.BitOrAssign.html
//!
//! ## Treat it like a `HashMap`
//!
//! `Counter<T, N>` implements [`Deref`]`<Target=HashMap<T, N>>` and
//! [`DerefMut`]`<Target=HashMap<T, N>>`, which means that you can perform any operations
//! on it which are valid for a [`HashMap`].
//!
//! [`HashMap`]: https://doc.rust-lang.org/std/collections/struct.HashMap.html
//! [`Deref`]: https://doc.rust-lang.org/stable/std/ops/trait.Deref.html
//! [`DerefMut`]: https://doc.rust-lang.org/stable/std/ops/trait.DerefMut.html
//!
//! ```rust
//! # use counter::Counter;
//! let mut counter = "aa-bb-cc".chars().collect::<Counter<_>>();
//! counter.remove(&'-');
//! assert!(counter == "aabbcc".chars().collect::<Counter<_>>());
//! ```
//!
//! Note that `Counter<T, N>` itself implements [`Index`]. `Counter::index` returns a reference to
//! a [`Zero::zero`] value for missing keys.
//!
//! [`Index`]: https://doc.rust-lang.org/stable/std/ops/trait.Index.html
//! [`Zero::zero`]: https://docs.rs/num-traits/latest/num_traits/identities/trait.Zero.html#tymethod.zero
//!
//! ```rust
//! # use counter::Counter;
//! let counter = "aaa".chars().collect::<Counter<_>>();
//! assert_eq!(counter[&'b'], 0);
//! // panics
//! // assert_eq!((*counter)[&'b'], 0);
//! ```
//!
//! # Advanced Usage
//!
//! ## Count any iterable which is `Hash + Eq`
//!
//! You can't use the `most_common*` functions unless `T` is also [`Clone`], but simple counting
//! works fine on a minimal data type.
//!
//! [`Clone`]: https://doc.rust-lang.org/stable/std/clone/trait.Clone.html
//!
//! ```rust
//! # use counter::Counter;
//! #[derive(Debug, Hash, PartialEq, Eq)]
//! struct Inty {
//!     i: usize,
//! }
//!
//! impl Inty {
//!     pub fn new(i: usize) -> Inty {
//!         Inty { i: i }
//!     }
//! }
//!
//! // <https://en.wikipedia.org/wiki/867-5309/Jenny>
//! let intys = vec![
//!     Inty::new(8),
//!     Inty::new(0),
//!     Inty::new(0),
//!     Inty::new(8),
//!     Inty::new(6),
//!     Inty::new(7),
//!     Inty::new(5),
//!     Inty::new(3),
//!     Inty::new(0),
//!     Inty::new(9),
//! ];
//!
//! let inty_counts = intys.iter().collect::<Counter<_>>();
//! println!("{:?}", inty_counts);
//! // {Inty { i: 8 }: 2, Inty { i: 0 }: 3, Inty { i: 9 }: 1, Inty { i: 3 }: 1,
//! //  Inty { i: 7 }: 1, Inty { i: 6 }: 1, Inty { i: 5 }: 1}
//! assert!(inty_counts.get(&Inty { i: 8 }) == Some(&2));
//! assert!(inty_counts.get(&Inty { i: 0 }) == Some(&3));
//! assert!(inty_counts.get(&Inty { i: 6 }) == Some(&1));
//! ```
//!
//! ## Use your own type for the count
//!
//! Sometimes [`usize`] just isn't enough. If you find yourself overflowing your
//! machine's native size, you can use your own type. Here, we use an [`i8`], but
//! you can use most numeric types, including bignums, as necessary.
//!
//! [`usize`]: https://doc.rust-lang.org/stable/std/primitive.usize.html
//! [`i8`]: https://doc.rust-lang.org/stable/std/primitive.i8.html
//!
//! ```rust
//! # use counter::Counter;
//! # use std::collections::HashMap;
//! let counter: Counter<_, i8> = "abbccc".chars().collect();
//! let expected: HashMap<char, i8> = [('a', 1), ('b', 2), ('c', 3)].iter().cloned().collect();
//! assert!(counter.into_map() == expected);
//! ```

#![allow(clippy::must_use_candidate)]
mod impls;

use num_traits::{One, Zero};

use std::collections::{BinaryHeap, HashMap};
use std::hash::Hash;
use std::iter;
use std::ops::{AddAssign, SubAssign};
#[cfg(test)]
mod unit_tests;

type CounterMap<T, N> = HashMap<T, N>;

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Counter<T: Hash + Eq, N = usize> {
    map: CounterMap<T, N>,
    // necessary for `Index::index` since we cannot declare generic `static` variables.
    zero: N,
}

impl<T, N> Counter<T, N>
where
    T: Hash + Eq,
{
    /// Consumes this counter and returns a [`HashMap`] mapping the items to the counts.
    ///
    /// [`HashMap`]: https://doc.rust-lang.org/stable/std/collections/struct.HashMap.html
    pub fn into_map(self) -> HashMap<T, N> {
        self.map
    }

    /// Returns the sum of the counts.
    ///
    /// Use [`len`] to get the number of elements in the counter and use `total` to get the sum of
    /// their counts.
    ///
    /// [`len`]: struct.Counter.html#method.len
    ///
    /// # Examples
    ///
    /// ```
    /// # use counter::Counter;
    /// let counter = "abracadabra".chars().collect::<Counter<_>>();
    /// assert_eq!(counter.total::<usize>(), 11);
    /// assert_eq!(counter.len(), 5);
    /// ```
    pub fn total<'a, S>(&'a self) -> S
    where
        S: iter::Sum<&'a N>,
    {
        self.map.values().sum()
    }
}

impl<T, N> Counter<T, N>
where
    T: Hash + Eq,
    N: AddAssign + Zero + One,
{
    /// Add the counts of the elements from the given iterable to this counter.
    pub fn update<I>(&mut self, iterable: I)
    where
        I: IntoIterator<Item = T>,
    {
        for item in iterable {
            let entry = self.map.entry(item).or_insert_with(N::zero);
            *entry += N::one();
        }
    }
}

impl<T, N> Counter<T, N>
where
    T: Hash + Eq,
    N: PartialOrd + SubAssign + Zero + One,
{
    /// Remove the counts of the elements from the given iterable to this counter.
    ///
    /// Non-positive counts are automatically removed.
    ///
    /// ```rust
    /// # use counter::Counter;
    /// # use std::collections::HashMap;
    /// let mut counter = "abbccc".chars().collect::<Counter<_>>();
    /// counter.subtract("abba".chars());
    /// let expect = [('c', 3)].iter().cloned().collect::<HashMap<_, _>>();
    /// assert_eq!(counter.into_map(), expect);
    /// ```
    pub fn subtract<I>(&mut self, iterable: I)
    where
        I: IntoIterator<Item = T>,
    {
        for item in iterable {
            let mut remove = false;
            if let Some(entry) = self.map.get_mut(&item) {
                if *entry > N::zero() {
                    *entry -= N::one();
                }
                remove = *entry == N::zero();
            }
            if remove {
                self.map.remove(&item);
            }
        }
    }
}

impl<T, N> Counter<T, N>
where
    T: Hash + Eq + Clone,
    N: Clone + Ord,
{
    /// Create a vector of `(elem, frequency)` pairs, sorted most to least common.
    ///
    /// ```rust
    /// # use counter::Counter;
    /// let mc = "pappaopolo".chars().collect::<Counter<_>>().most_common();
    /// let expected = vec![('p', 4), ('o', 3), ('a', 2), ('l', 1)];
    /// assert_eq!(mc, expected);
    /// ```
    ///
    /// Note that the ordering of duplicates is unstable.
    pub fn most_common(&self) -> Vec<(T, N)> {
        use std::cmp::Ordering;
        self.most_common_tiebreaker(|_a, _b| Ordering::Equal)
    }

    /// Create a vector of `(elem, frequency)` pairs, sorted most to least common.
    ///
    /// In the event that two keys have an equal frequency, use the supplied ordering function
    /// to further arrange the results.
    ///
    /// For example, we can sort reverse-alphabetically:
    ///
    /// ```rust
    /// # use counter::Counter;
    /// let counter = "eaddbbccc".chars().collect::<Counter<_>>();
    /// let by_common = counter.most_common_tiebreaker(|&a, &b| b.cmp(&a));
    /// let expected = vec![('c', 3), ('d', 2), ('b', 2), ('e', 1), ('a', 1)];
    /// assert_eq!(by_common, expected);
    /// ```
    pub fn most_common_tiebreaker<F>(&self, mut tiebreaker: F) -> Vec<(T, N)>
    where
        F: FnMut(&T, &T) -> ::std::cmp::Ordering,
    {
        let mut items = self
            .map
            .iter()
            .map(|(key, count)| (key.clone(), count.clone()))
            .collect::<Vec<_>>();
        items.sort_unstable_by(|(a_item, a_count), (b_item, b_count)| {
            b_count
                .cmp(a_count)
                .then_with(|| tiebreaker(a_item, b_item))
        });
        items
    }
}

impl<T, N> Counter<T, N>
where
    T: Hash + Eq + Clone + Ord,
    N: Clone + Ord,
{
    /// Create a vector of `(elem, frequency)` pairs, sorted most to least common.
    ///
    /// In the event that two keys have an equal frequency, use the natural ordering of the keys
    /// to further sort the results.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use counter::Counter;
    /// let mc = "abracadabra".chars().collect::<Counter<_>>().most_common_ordered();
    /// let expect = vec![('a', 5), ('b', 2), ('r', 2), ('c', 1), ('d', 1)];
    /// assert_eq!(mc, expect);
    /// ```
    ///
    /// # Time complexity
    ///
    /// *O*(*n* \* log *n*), where *n* is the number of items in the counter.  If all you want is
    /// the top *k* items and *k* < *n* then it can be more efficient to use
    /// [`k_most_common_ordered`].
    ///
    /// [`k_most_common_ordered`]: Counter::k_most_common_ordered
    pub fn most_common_ordered(&self) -> Vec<(T, N)> {
        self.most_common_tiebreaker(Ord::cmp)
    }

    /// Returns the `k` most common items in decreasing order of their counts.
    ///
    /// The returned vector is the same as would be obtained by calling `most_common_ordered` and
    /// then truncating the result to length `k`.  In particular, items with the same count are
    /// sorted in *increasing* order of their keys.  Further, if `k` is greater than the length of
    /// the counter then the returned vector will have length equal to that of the counter, not
    /// `k`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// # use counter::Counter;
    /// let counter: Counter<_> = "abracadabra".chars().collect();
    /// let top3 = counter.k_most_common_ordered(3);
    /// assert_eq!(top3, vec![('a', 5), ('b', 2), ('r', 2)]);
    /// ```
    ///
    /// # Time complexity
    ///
    /// This method can be much more efficient than [`most_common_ordered`] when *k* is much
    /// smaller than the length of the counter *n*.  When *k* = 1 the algorithm is equivalent
    /// to finding the minimum (or maximum) of *n* items, which requires *n* \- 1 comparisons.  For
    /// a fixed value of *k* > 1, the number of comparisons scales with *n* as *n* \+ *O*(log *n*)
    /// and the number of swaps scales as *O*(log *n*).  As *k* approaches *n*, this algorithm
    /// approaches a heapsort of the *n* items, which has complexity *O*(*n* \* log *n*).
    ///
    /// For values of *k* close to *n* the sorting algorithm used by [`most_common_ordered`] will
    /// generally be faster than the heapsort used by this method by a small constant factor.
    /// Exactly where the crossover point occurs will depend on several factors.  For small *k*
    /// choose this method.  If *k* is a substantial fraction of *n*, it may be that
    /// [`most_common_ordered`] is faster.  If performance matters in your application then it may
    /// be worth experimenting to see which of the two methods is faster.
    ///
    /// [`most_common_ordered`]: Counter::most_common_ordered
    #[allow(clippy::missing_panics_doc)] // current implementation does not panic
    pub fn k_most_common_ordered(&self, k: usize) -> Vec<(T, N)> {
        use std::cmp::Reverse;

        if k == 0 {
            return vec![];
        }

        // The quicksort implementation used by `most_common_ordered()` is generally faster than
        // the heapsort used below when sorting the entire counter.
        if k >= self.map.len() {
            return self.most_common_ordered();
        }

        // Clone the counts as we iterate over the map to eliminate an extra indirection when
        // comparing counts.  This will be an improvement in the typical case where `N: Copy`.
        // Defer cloning the keys until we have selected the top `k` items so that we clone only
        // `k` keys instead of all of them.
        let mut items = self.map.iter().map(|(t, n)| (Reverse(n.clone()), t));

        // Step 1. Make a heap out of the first `k` items; this makes O(k) comparisons.
        let mut heap: BinaryHeap<_> = items.by_ref().take(k).collect();

        // Step 2. Successively compare each of the remaining `n - k` items to the top of the heap,
        // replacing the root (and subsequently sifting down) whenever the item is less than the
        // root.  This takes at most n - k + k * (1 + log2(k)) * (H(n) - H(k)) comparisons, where
        // H(i) is the ith [harmonic number](https://en.wikipedia.org/wiki/Harmonic_number).  For
        // fixed `k`, this scales as *n* + *O*(log(*n*)).
        items.for_each(|item| {
            // If `items` is nonempty at this point then we know the heap contains `k > 0`
            // elements.
            let mut root = heap.peek_mut().expect("the heap is empty");
            if *root > item {
                *root = item;
            }
        });

        // Step 3. Sort the items in the heap with the second phases of heapsort.  The number of
        // comparisons is 2 * k * log2(k) + O(k).
        heap.into_sorted_vec()
            .into_iter()
            .map(|(Reverse(n), t)| (t.clone(), n))
            .collect()
    }
}

impl<T, N> Counter<T, N>
where
    T: Hash + Eq,
    N: PartialOrd + Zero,
{
    /// Test whether this counter is a superset of another counter.
    /// This is true if for all elements in this counter and the other,
    /// the count in this counter is greater than or equal to the count in the other.
    ///
    /// `c.is_superset(&d);` -> `c.iter().all(|(x, n)| n >= d[x]) && d.iter().all(|(x, n)| c[x] >= n)`
    ///
    /// ```rust
    /// # use counter::Counter;
    /// # use std::collections::HashMap;
    /// let c = "aaabbc".chars().collect::<Counter<_>>();
    /// let mut d = "abb".chars().collect::<Counter<_>>();
    ///
    /// assert!(c.is_superset(&d));
    /// d[&'e'] = 1;
    /// assert!(!c.is_superset(&d));
    /// ```
    pub fn is_superset(&self, other: &Self) -> bool {
        // need to test keys from both counters, because if N is signed, counts in `self`
        // could be < 0 for elements missing in `other`. For the unsigned case, only elements
        // from `other` would need to be tested.
        self.keys()
            .chain(other.keys())
            .all(|key| self[key] >= other[key])
    }

    /// Test whether this counter is a subset of another counter.
    /// This is true if for all elements in this counter and the other,
    /// the count in this counter is less than or equal to the count in the other.
    ///
    /// `c.is_subset(&d);` -> `c.iter().all(|(x, n)| n <= d[x]) && d.iter().all(|(x, n)| c[x] <= n)`
    ///
    /// ```rust
    /// # use counter::Counter;
    /// # use std::collections::HashMap;
    /// let mut c = "abb".chars().collect::<Counter<_>>();
    /// let mut d = "aaabbc".chars().collect::<Counter<_>>();
    ///
    /// assert!(c.is_subset(&d));
    /// c[&'e'] = 1;
    /// assert!(!c.is_subset(&d));
    /// ```
    pub fn is_subset(&self, other: &Self) -> bool {
        // need to test keys from both counters, because if N is signed, counts in `other`
        // could be < 0 for elements missing in `self`. For the unsigned case, only elements
        // from `self` would need to be tested.
        self.keys()
            .chain(other.keys())
            .all(|key| self[key] <= other[key])
    }
}
