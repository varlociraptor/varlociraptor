use crate::Counter;

use std::collections::HashMap;
use std::hash::Hash;
use std::ops::{Deref, DerefMut};

type CounterMap<T, N> = HashMap<T, N>;

impl<T, N> Deref for Counter<T, N>
where
    T: Hash + Eq,
{
    type Target = CounterMap<T, N>;
    fn deref(&self) -> &CounterMap<T, N> {
        &self.map
    }
}

impl<T, N> DerefMut for Counter<T, N>
where
    T: Hash + Eq,
{
    fn deref_mut(&mut self) -> &mut CounterMap<T, N> {
        &mut self.map
    }
}
