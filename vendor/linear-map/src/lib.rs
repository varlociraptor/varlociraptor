//! A map implemented by searching linearly in a vector.
//!
//! See the [`LinearMap`](struct.LinearMap.html) type for details.

#![deny(missing_docs)]

// Optional Serde support
#[cfg(feature = "serde_impl")]
pub mod serde;
pub mod set;

use std::borrow::Borrow;
use std::fmt::{self, Debug};
use std::iter;
use std::mem;
use std::ops;
use std::slice;
use std::vec;

use self::Entry::{Occupied, Vacant};

/// A map implemented by searching linearly in a vector.
///
/// `LinearMap`'s keys are compared using the [`Eq`][eq] trait. All search operations
/// (`contains_key`, `get`, `get_mut`, `insert`, and `remove`) run in `O(n)` time, making this
/// implementation suitable only for small numbers of keys. The ordering of the keys in the
/// underlying vector is arbitrary.
///
/// It is a logic error for a key to be modified in such a way that the key's equality, as
/// determined by the [`Eq`][eq] trait, changes while it is in the map. This is normally only
/// possible through [`Cell`][cell], [`RefCell`][ref_cell], global state, I/O, or unsafe code.
///
/// [cell]: https://doc.rust-lang.org/nightly/std/cell/struct.Cell.html
/// [eq]: https://doc.rust-lang.org/nightly/std/cmp/trait.Eq.html
/// [ref_cell]: https://doc.rust-lang.org/nightly/std/cell/struct.RefCell.html
///
/// # Example
///
/// ```
/// use linear_map::LinearMap;
///
/// // type inference lets us omit an explicit type signature (which
/// // would be `LinearMap<&str, &str>` in this example).
/// let mut book_reviews = LinearMap::new();
///
/// // review some books.
/// book_reviews.insert("Adventures of Huckleberry Finn",    "My favorite book.");
/// book_reviews.insert("Grimms' Fairy Tales",               "Masterpiece.");
/// book_reviews.insert("Pride and Prejudice",               "Very enjoyable.");
/// book_reviews.insert("The Adventures of Sherlock Holmes", "Eye lyked it alot.");
///
/// // check for a specific one.
/// if !book_reviews.contains_key("Les Misérables") {
///     println!("We've got {} reviews, but Les Misérables ain't one.",
///              book_reviews.len());
/// }
///
/// // oops, this review has a lot of spelling mistakes. let's delete it.
/// book_reviews.remove("The Adventures of Sherlock Holmes");
///
/// // look up the values associated with some keys.
/// let to_find = ["Pride and Prejudice", "Alice's Adventure in Wonderland"];
/// for book in &to_find {
///     match book_reviews.get(book) {
///         Some(review) => println!("{}: {}", book, review),
///         None => println!("{} is unreviewed.", book)
///     }
/// }
///
/// // iterate over everything.
/// for (book, review) in &book_reviews {
///     println!("{}: \"{}\"", book, review);
/// }
/// ```
pub struct LinearMap<K, V> {
    storage: Vec<(K, V)>,
}

impl<K: Eq, V> LinearMap<K, V> {
    /// Creates an empty map. This method does not allocate.
    pub fn new() -> Self {
        LinearMap { storage: vec![] }
    }

    /// Creates an empty map with the given initial capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        LinearMap { storage: Vec::with_capacity(capacity) }
    }

    /// Returns the number of elements the map can hold without reallocating.
    pub fn capacity(&self) -> usize {
        self.storage.capacity()
    }

    /// Reserves capacity for at least `additional` more to be inserted in the
    /// map. The collection may reserve more space to avoid frequent
    /// reallocations.
    ///
    /// # Panics
    ///
    /// Panics if the new allocation size overflows `usize`.
    pub fn reserve(&mut self, additional: usize) {
        self.storage.reserve(additional);
    }

    /// Reserves the minimum capacity for exactly `additional` more elemnnts to
    /// be inserted in the map.
    ///
    /// Note that the allocator may give the collection more space than it
    /// requests. Therefore capacity cannot be relied upon to be precisely
    /// minimal. Prefer `reserve` if future insertions are expected.
    ///
    /// # Panics
    ///
    /// Panics if the new capacity overflows `usize`.
    pub fn reserve_exact(&mut self, additional: usize) {
        self.storage.reserve_exact(additional);
    }

    /// Shrinks the capacity of the map as much as possible.
    ///
    /// It will drop down as close as possible to the current length but the
    /// allocator may still inform the map that there is more space than
    /// necessary. Therefore capacity cannot be relid upon to be minimal.
    pub fn shrink_to_fit(&mut self) {
        self.storage.shrink_to_fit();
    }

    /// Returns the number of elements in the map.
    pub fn len(&self) -> usize {
        self.storage.len()
    }

    /// Returns true if the map contains no elements.
    pub fn is_empty(&self) -> bool {
        self.storage.is_empty()
    }

    /// Clears the map, removing all elements. Keeps the allocated memory for
    /// reuse.
    pub fn clear(&mut self) {
        self.storage.clear();
    }

    /// Scan through the map and keep those key-value pairs where the
    /// closure returns `true`.
    ///
    /// The order the elements are visited is not specified.
    pub fn retain<F>(&mut self, mut keep_fn: F)
    where F: FnMut(&K, &mut V) -> bool {
        let mut del = 0;
        {
            let v = &mut *self.storage;
            for i in 0..v.len() {
                if !keep_fn(&v[i].0, &mut v[i].1) {
                    del += 1;
                } else if del > 0 {
                    v.swap(i - del, i);
                }
            }
        }
        if del > 0 {
            let len = self.storage.len();
            self.storage.truncate(len - del);
        }
    }

    /// Removes all key-value pairs from the map and returns an iterator that yields them in
    /// arbitrary order.
    ///
    /// All key-value pairs are removed even if the iterator is not exhausted. However, the
    /// behavior of this method is unspecified if the iterator is leaked.
    ///
    /// The iterator's item type is `(K, V)`.
    pub fn drain(&mut self) -> Drain<K, V> {
        Drain { iter: self.storage.drain(..) }
    }

    /// Returns an iterator yielding references to the map's keys and their corresponding values in
    /// arbitrary order.
    ///
    /// The iterator's item type is `(&K, &V)`.
    pub fn iter(&self) -> Iter<K, V> {
        Iter { iter: self.storage.iter() }
    }

    /// Returns an iterator yielding references to the map's keys and mutable references to their
    /// corresponding values in arbitrary order.
    ///
    /// The iterator's item type is `(&K, &mut V)`.
    pub fn iter_mut(&mut self) -> IterMut<K, V> {
        IterMut { iter: self.storage.iter_mut() }
    }

    /// Returns an iterator yielding references to the map's keys in arbitrary order.
    ///
    /// The iterator's item type is `&K`.
    pub fn keys(&self) -> Keys<K, V> {
        Keys { iter: self.iter() }
    }

    /// Returns an iterator yielding references to the map's values in arbitrary order.
    ///
    /// The iterator's item type is `&V`.
    pub fn values(&self) -> Values<K, V> {
        Values { iter: self.iter() }
    }

    /// Returns a reference to the value in the map whose key is equal to the given key.
    ///
    /// Returns `None` if the map contains no such key.
    ///
    /// The given key may be any borrowed form of the map's key type, but `Eq` on the borrowed form
    /// *must* match that of the key type.
    pub fn get<Q: ?Sized + Eq>(&self, key: &Q) -> Option<&V> where K: Borrow<Q> {
        for (k, v) in self {
            if key == k.borrow() {
                return Some(v);
            }
        }
        None
    }

    /// Returns a mutable reference to the value in the map whose key is equal to the given key.
    ///
    /// Returns `None` if the map contains no such key.
    ///
    /// The given key may be any borrowed form of the map's key type, but `Eq` on the borrowed form
    /// *must* match that of the key type.
    pub fn get_mut<Q: ?Sized + Eq>(&mut self, key: &Q) -> Option<&mut V> where K: Borrow<Q> {
        for (k, v) in self {
            if key == k.borrow() {
                return Some(v);
            }
        }
        None
    }

    /// Checks if the map contains a key that is equal to the given key.
    ///
    /// The given key may be any borrowed form of the map's key type, but `Eq` on the borrowed form
    /// *must* match that of the key type.
    pub fn contains_key<Q: ?Sized + Eq>(&self, key: &Q) -> bool where K: Borrow<Q> {
        self.get(key).is_some()
    }

    /// Inserts a key-value pair into the map.
    ///
    /// Returns `None` if the map did not contain a key that is equal to the given key.
    ///
    /// If the map did contain such a key, its corresponding value is replaced with the given
    /// value, and the old value is returned. The key is not updated, though. This matters for
    /// values that can be `==` without being identical. See the [standard library's documentation]
    /// [std] for more details.
    ///
    /// [std]: https://doc.rust-lang.org/nightly/std/collections/index.html#insert-and-complex-keys
    pub fn insert(&mut self, key: K, value: V) -> Option<V> {
        match self.entry(key) {
            Occupied(mut e) => Some(e.insert(value)),
            Vacant(e) => { e.insert(value); None }
        }
    }

    /// Removes the key in the map that is equal to the given key and returns its corresponding
    /// value.
    ///
    /// Returns `None` if the map contained no such key.
    ///
    /// The given key may be any borrowed form of the map's key type, but `Eq` on the borrowed form
    /// *must* match that of the key type.
    pub fn remove<Q: ?Sized + Eq>(&mut self, key: &Q) -> Option<V> where K: Borrow<Q> {
        for i in 0..self.storage.len() {
            if self.storage[i].0.borrow() == key {
                return Some(self.storage.swap_remove(i).1);
            }
        }
        None
    }

    /// Returns the given key's corresponding entry in the map for in-place manipulation.
    pub fn entry(&mut self, key: K) -> Entry<K, V> {
        match self.storage.iter().position(|&(ref k, _)| key == *k) {
            None => Vacant(VacantEntry {
                map: self,
                key: key
            }),
            Some(index) => Occupied(OccupiedEntry {
                map: self,
                index: index
            })
        }
    }
}

impl<K: Clone, V: Clone> Clone for LinearMap<K, V> {
    fn clone(&self) -> Self {
        LinearMap { storage: self.storage.clone() }
    }

    fn clone_from(&mut self, other: &Self) {
        self.storage.clone_from(&other.storage);
    }
}

impl<K: Eq + Debug, V: Debug> Debug for LinearMap<K, V> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.debug_map().entries(self).finish()
    }
}

impl<K: Eq, V> Default for LinearMap<K, V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: Eq, V> Extend<(K, V)> for LinearMap<K, V> {
    fn extend<I: IntoIterator<Item = (K, V)>>(&mut self, key_values: I) {
        for (key, value) in key_values { self.insert(key, value); }
    }
}

impl<K: Eq, V> iter::FromIterator<(K, V)> for LinearMap<K, V> {
    fn from_iter<I: IntoIterator<Item = (K, V)>>(key_values: I) -> Self {
        let mut map = Self::new();
        map.extend(key_values);
        map
    }
}

impl<'a, K: Eq + Borrow<Q>, V, Q: ?Sized + Eq> ops::Index<&'a Q> for LinearMap<K, V> {
    type Output = V;

    fn index(&self, key: &'a Q) -> &V {
        self.get(key).expect("key not found")
    }
}

impl<K: Eq, V: PartialEq> PartialEq for LinearMap<K, V> {
    fn eq(&self, other: &Self) -> bool {
        if self.len() != other.len() {
            return false;
        }

        for (key, value) in self {
            if other.get(key) != Some(value) {
                return false;
            }
        }

        true
    }
}

impl<K: Eq, V: Eq> Eq for LinearMap<K, V> {}

impl<K: Eq, V> Into<Vec<(K, V)>> for LinearMap<K, V> {
    fn into(self) -> Vec<(K, V)> {
        self.storage
    }
}

/// Creates a `LinearMap` from a list of key-value pairs.
///
/// The created `LinearMap` has a capacity set to the number of entries provided.
///
/// # Example
///
/// ```
/// #[macro_use] extern crate linear_map;
/// # fn main() {
///
/// let map = linear_map!{
///     "a" => 1,
///     "b" => 2,
/// };
/// assert_eq!(map["a"], 1);
/// assert_eq!(map["b"], 2);
/// assert_eq!(map.get("c"), None);
/// # }
/// ```
#[macro_export]
macro_rules! linear_map {
    ($($key:expr => $value:expr,)+) => { linear_map!($($key => $value),+) };
    ($($key:expr => $value:expr),*) => {
        {
            let _cap = <[&str]>::len(&[$(stringify!($key)),*]);
            let mut _map = $crate::LinearMap::with_capacity(_cap);
            $(
                _map.insert($key, $value);
            )*
            _map
        }
    };
}

/// A view into a single occupied location in a `LinearMap`.
///
/// See [`LinearMap::entry`](struct.LinearMap.html#method.entry) for details.
pub struct OccupiedEntry<'a, K: 'a, V: 'a> {
    map: &'a mut LinearMap<K, V>,
    index: usize,
}

/// A view into a single vacant location in a `LinearMap`.
///
/// See [`LinearMap::entry`](struct.LinearMap.html#method.entry) for details.
pub struct VacantEntry<'a, K: 'a, V: 'a> {
    map: &'a mut LinearMap<K, V>,
    key: K,
}

/// A view into a single entry in a `LinearMap`.
///
/// See [`LinearMap::entry`](struct.LinearMap.html#method.entry) for details.
pub enum Entry<'a, K: 'a, V: 'a> {
    /// An occupied entry.
    Occupied(OccupiedEntry<'a, K, V>),

    /// A vacant entry.
    Vacant(VacantEntry<'a, K, V>)
}

impl<'a, K, V> Entry<'a, K, V> {
    /// Ensures that the entry is occupied by inserting the given value if it is vacant.
    ///
    /// Returns a mutable reference to the entry's value.
    pub fn or_insert(self, default: V) -> &'a mut V {
        match self {
            Occupied(entry) => entry.into_mut(),
            Vacant(entry) => entry.insert(default)
        }
    }

    /// Ensures that the entry is occupied by inserting the the result of the given function if it
    /// is vacant.
    ///
    /// Returns a mutable reference to the entry's value.
    pub fn or_insert_with<F: FnOnce() -> V>(self, default: F) -> &'a mut V {
        match self {
            Occupied(entry) => entry.into_mut(),
            Vacant(entry) => entry.insert(default())
        }
    }
}

impl<'a, K, V> OccupiedEntry<'a, K, V> {
    /// Returns a reference to the entry's value.
    pub fn get(&self) -> &V {
        &self.map.storage[self.index].1
    }

    /// Returns a mutable reference to the entry's value.
    pub fn get_mut(&mut self) -> &mut V {
        &mut self.map.storage[self.index].1
    }

    /// Returns a mutable reference to the entry's value with the same lifetime as the map.
    pub fn into_mut(self) -> &'a mut V {
        &mut self.map.storage[self.index].1
    }

    /// Replaces the entry's value with the given one and returns the previous value.
    pub fn insert(&mut self, value: V) -> V {
        mem::replace(self.get_mut(), value)
    }

    /// Removes the entry from the map and returns its value.
    pub fn remove(self) -> V {
        self.map.storage.swap_remove(self.index).1
    }
}

impl<'a, K, V> VacantEntry<'a, K, V> {
    /// Inserts the entry into the map with the given value.
    ///
    /// Returns a mutable reference to the entry's value with the same lifetime as the map.
    pub fn insert(self, value: V) -> &'a mut V {
        self.map.storage.push((self.key, value));
        &mut self.map.storage.last_mut().unwrap().1
    }
}

/// A consuming iterator over a `LinearMap`.
///
/// The iterator's order is arbitrary.
///
/// Acquire through [`IntoIterator`](struct.LinearMap.html#method.into_iter).
pub struct IntoIter<K, V> {
    iter: vec::IntoIter<(K, V)>,
}

impl<K, V> Iterator for IntoIter<K, V> {
    type Item = (K, V);

    fn next(&mut self) -> Option<(K, V)> {
        self.iter.next()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}

impl<K, V> DoubleEndedIterator for IntoIter<K, V> {
    fn next_back(&mut self) -> Option<(K, V)> {
        self.iter.next_back()
    }
}

impl<K, V> ExactSizeIterator for IntoIter<K, V> {
    fn len(&self) -> usize {
        self.iter.len()
    }
}

/// A draining iterator over a `LinearMap`.
///
/// See [`LinearMap::drain`](struct.LinearMap.html#method.drain) for details.
pub struct Drain<'a, K: 'a, V: 'a> {
    iter: vec::Drain<'a, (K, V)>,
}

/// An iterator yielding references to a `LinearMap`'s keys and their corresponding values.
///
/// See [`LinearMap::iter`](struct.LinearMap.html#method.iter) for details.
pub struct Iter<'a, K: 'a, V: 'a> {
    iter: slice::Iter<'a, (K, V)>,
}

/// An iterator yielding references to a `LinearMap`'s keys and mutable references to their
/// corresponding values.
///
/// See [`LinearMap::iter_mut`](struct.LinearMap.html#method.iter_mut) for details.
pub struct IterMut<'a, K: 'a, V: 'a> {
    iter: slice::IterMut<'a, (K, V)>,
}

/// An iterator yielding references to a `LinearMap`'s keys in arbitrary order.
///
/// See [`LinearMap::keys`](struct.LinearMap.html#method.keys) for details.
pub struct Keys<'a, K: 'a, V: 'a> {
    iter: Iter<'a, K, V>,
}

/// An iterator yielding references to a `LinearMap`'s values in arbitrary order.
///
/// See [`LinearMap::values`](struct.LinearMap.html#method.values) for details.
pub struct Values<'a, K: 'a, V: 'a> {
    iter: Iter<'a, K, V>,
}

macro_rules! impl_iter {($typ:ty, $item:ty, $map:expr) => {
    impl<'a, K, V> Iterator for $typ {
        type Item = $item;

        fn next(&mut self) -> Option<Self::Item> {
            self.iter.next().map($map)
        }

        fn size_hint(&self) -> (usize, Option<usize>) {
            self.iter.size_hint()
        }
    }

    impl<'a, K, V> DoubleEndedIterator for $typ {
        fn next_back(&mut self) -> Option<Self::Item> {
            self.iter.next_back().map($map)
        }
    }

    impl<'a, K, V> ExactSizeIterator for $typ {
        fn len(&self) -> usize {
            self.iter.len()
        }
    }
}}
impl_iter!{Drain<'a,K,V>,  (K,V),  |e| e }
impl_iter!{Iter<'a,K,V>,  (&'a K, &'a V),  |e| (&e.0, &e.1) }
impl_iter!{IterMut<'a,K,V>,  (&'a K, &'a mut V),  |e| (&e.0, &mut e.1) }
impl_iter!{Keys<'a,K,V>,  &'a K,  |e| e.0 }
impl_iter!{Values<'a,K,V>,  &'a V,  |e| e.1 }

impl<'a, K, V> Clone for Iter<'a, K, V> {
    fn clone(&self) -> Self {
        Iter { iter: self.iter.clone() }
    }
}

impl<'a, K, V> Clone for Keys<'a, K, V> {
    fn clone(&self) -> Self {
        Keys { iter: self.iter.clone() }
    }
}

impl<'a, K, V> Clone for Values<'a, K, V> {
    fn clone(&self) -> Self {
        Values { iter: self.iter.clone() }
    }
}

impl<K: Eq, V> IntoIterator for LinearMap<K, V> {
    type Item = (K, V);
    type IntoIter = IntoIter<K, V>;

    fn into_iter(self) -> IntoIter<K, V> {
        IntoIter { iter: self.storage.into_iter() }
    }
}

impl<'a, K: Eq, V> IntoIterator for &'a LinearMap<K, V> {
    type Item = (&'a K, &'a V);
    type IntoIter = Iter<'a, K, V>;

    fn into_iter(self) -> Iter<'a, K, V> {
        self.iter()
    }
}

impl<'a, K: Eq, V> IntoIterator for &'a mut LinearMap<K, V> {
    type Item = (&'a K, &'a mut V);
    type IntoIter = IterMut<'a, K, V>;

    fn into_iter(self) -> IterMut<'a, K, V> {
        self.iter_mut()
    }
}

#[allow(dead_code)]
fn assert_covariance() {
    fn a<'a, K, V>(x: LinearMap<&'static K, &'static V>) -> LinearMap<&'a K, &'a V> { x }

    fn b<'a, K, V>(x: IntoIter<&'static K, &'static V>) -> IntoIter<&'a K, &'a V> { x }

    fn c<'i, 'a, K, V>(x: Iter<'i, &'static K, &'static V>) -> Iter<'i, &'a K, &'a V> { x }

    fn d<'i, 'a, K, V>(x: Keys<'i, &'static K, &'static V>) -> Keys<'i, &'a K, &'a V> { x }

    fn e<'i, 'a, K, V>(x: Values<'i, &'static K, &'static V>) -> Values<'i, &'a K, &'a V> { x }
}
