// Copyright 2019 MaidSafe.net limited.
//
// This SAFE Network Software is licensed to you under the MIT license <LICENSE-MIT
// http://opensource.org/licenses/MIT> or the Modified BSD license <LICENSE-BSD
// https://opensource.org/licenses/BSD-3-Clause>, at your option. This file may not be copied,
// modified, or distributed except according to those terms. Please review the Licences for the
// specific language governing permissions and limitations relating to use of the SAFE Network
// Software.

//! Misc LRU cache iterators.

#[cfg(feature = "sn_fake_clock")]
use sn_fake_clock::FakeClock as Instant;
use std::collections::{BTreeMap, VecDeque};
use std::time::Duration;
#[cfg(not(feature = "sn_fake_clock"))]
use std::time::Instant;

/// An iterator over an `LruCache`'s entries that updates the timestamps as values are traversed.
/// Values are produced in the most recently used order.
pub struct Iter<'a, Key, Value> {
    /// Reference to the iterated cache.
    map: &'a mut BTreeMap<Key, (Value, Instant)>,
    /// Ordered cache entry keys where the least recently used items are first.
    list: &'a mut VecDeque<Key>,
    lru_cache_ttl: Option<Duration>,
    /// Index in `list` of the previously used item.
    item_index: usize,
}

impl<'a, Key, Value> Iter<'a, Key, Value>
where
    Key: Ord,
{
    #[doc(hidden)]
    pub fn new(
        map: &'a mut BTreeMap<Key, (Value, Instant)>,
        list: &'a mut VecDeque<Key>,
        lru_cache_ttl: Option<Duration>,
    ) -> Self {
        let item_index = list.len();
        Self {
            map,
            list,
            lru_cache_ttl,
            item_index,
        }
    }

    /// Returns next unexpired item in the cache or `None` if no such items.
    /// Expired items are removed from the cache.
    fn next_unexpired(&mut self, now: Instant) -> Option<Key> {
        loop {
            self.item_index = self.item_index.checked_sub(1)?;
            let key = self.list.remove(self.item_index)?;
            let value = self.map.get(&key)?;

            if let Some(ttl) = self.lru_cache_ttl {
                if value.1 + ttl > now {
                    return Some(key);
                } else {
                    let _ = self.map.remove(&key);
                }
            } else {
                return Some(key);
            }
        }
    }
}

impl<'a, Key, Value> Iterator for Iter<'a, Key, Value>
where
    Key: Ord + Clone,
{
    type Item = (&'a Key, &'a Value);

    /// Returns the next element in the cache and moves it to the top of the cache.
    /// The most recently used items are yield first.
    #[allow(unsafe_code)]
    fn next(&mut self) -> Option<(&'a Key, &'a Value)> {
        let now = Instant::now();
        let key = self.next_unexpired(now)?;
        self.list.push_back(key);
        let key = self.list.back()?;
        let mut value = self.map.get_mut(&key)?;
        value.1 = now;

        unsafe {
            let key = std::mem::transmute::<&Key, &'a Key>(key);
            let value = std::mem::transmute::<&Value, &'a Value>(&value.0);
            Some((key, value))
        }
    }
}

/// Entry produced by `NotifyIter` that might be still valid or expired.
pub enum TimedEntry<'a, Key: 'a, Value: 'a> {
    /// Entry has not yet expired.
    Valid(&'a Key, &'a Value),
    /// Entry got expired and was evicted from the cache.
    Expired(Key, Value),
}

/// Much like `Iter` except will produce expired entries too where `Iter` silently drops them.
pub struct NotifyIter<'a, Key, Value> {
    /// Reference to the iterated cache.
    map: &'a mut BTreeMap<Key, (Value, Instant)>,
    /// Ordered cache entry keys where the least recently used items are first.
    list: &'a mut VecDeque<Key>,
    lru_cache_ttl: Option<Duration>,
    /// Index in `list` of the previously used item.
    item_index: usize,
}

impl<'a, Key, Value> NotifyIter<'a, Key, Value>
where
    Key: Ord + Clone,
{
    #[doc(hidden)]
    pub fn new(
        map: &'a mut BTreeMap<Key, (Value, Instant)>,
        list: &'a mut VecDeque<Key>,
        lru_cache_ttl: Option<Duration>,
    ) -> Self {
        let item_index = list.len();
        Self {
            map,
            list,
            lru_cache_ttl,
            item_index,
        }
    }
}

impl<'a, Key, Value> Iterator for NotifyIter<'a, Key, Value>
where
    Key: Ord + Clone,
{
    type Item = TimedEntry<'a, Key, Value>;

    /// Returns the next element in the cache and moves it to the top of the cache.
    /// The most recently used items are yield first.
    #[allow(unsafe_code)]
    fn next(&mut self) -> Option<Self::Item> {
        self.item_index = self.item_index.checked_sub(1)?;
        let key = self.list.remove(self.item_index)?;
        let mut value = self.map.get_mut(&key)?;
        let now = Instant::now();

        if let Some(ttl) = self.lru_cache_ttl {
            if value.1 + ttl <= now {
                let value = self.map.remove(&key)?;
                return Some(TimedEntry::Expired(key, value.0));
            }
        }

        self.list.push_back(key);
        let key = self.list.back()?;
        value.1 = now;
        unsafe {
            let key = std::mem::transmute::<&Key, &'a Key>(key);
            let value = std::mem::transmute::<&Value, &'a Value>(&value.0);
            Some(TimedEntry::Valid(key, value))
        }
    }
}

/// An iterator over an `LruCache`'s entries that does not modify the timestamp.
pub struct PeekIter<'a, Key, Value> {
    /// Reference to the iterated cache.
    map: &'a BTreeMap<Key, (Value, Instant)>,
    /// Ordered cache entry keys where the least recently used items are first.
    list: &'a VecDeque<Key>,
    lru_cache_ttl: Option<Duration>,
    /// Index in `list` of the previously used item.
    item_index: usize,
}

impl<'a, Key, Value> PeekIter<'a, Key, Value>
where
    Key: Ord,
{
    #[doc(hidden)]
    pub fn new(
        map: &'a BTreeMap<Key, (Value, Instant)>,
        list: &'a VecDeque<Key>,
        lru_cache_ttl: Option<Duration>,
    ) -> Self {
        let item_index = list.len();
        Self {
            map,
            list,
            lru_cache_ttl,
            item_index,
        }
    }

    /// Returns next unexpired item in the cache or `None` if no such items.
    fn next_unexpired(&mut self, now: Instant) -> Option<()> {
        loop {
            self.item_index = self.item_index.checked_sub(1)?;
            let value = self.map.get(&self.list[self.item_index])?;

            if let Some(ttl) = self.lru_cache_ttl {
                if value.1 + ttl > now {
                    return Some(());
                }
            } else {
                return Some(());
            }
        }
    }
}

impl<'a, Key, Value> Iterator for PeekIter<'a, Key, Value>
where
    Key: Ord + Clone,
{
    type Item = (&'a Key, &'a Value);

    /// Returns the next element in the cache that has not expired yet.
    /// The most recently used items are yield first.
    #[allow(unsafe_code)]
    fn next(&mut self) -> Option<(&'a Key, &'a Value)> {
        let now = Instant::now();
        self.next_unexpired(now)?;
        let key = &self.list[self.item_index];
        let value = self.map.get(&key)?;

        unsafe {
            let key = std::mem::transmute::<&Key, &'a Key>(key);
            let value = std::mem::transmute::<&Value, &'a Value>(&value.0);
            Some((key, value))
        }
    }
}
