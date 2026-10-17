#[macro_use]
extern crate linear_map;

use linear_map::LinearMap;
use linear_map::Entry::{Occupied, Vacant};

const TEST_CAPACITY: usize = 10;

#[test]
fn test_new() {
    let map: LinearMap<i32, i32> = LinearMap::new();
    assert_eq!(map.capacity(), 0);
    assert_eq!(map.len(), 0);
    assert!(map.is_empty());
}

#[test]
fn test_with_capacity() {
    let map: LinearMap<i32, i32> = LinearMap::with_capacity(TEST_CAPACITY);
    assert!(map.capacity() >= TEST_CAPACITY);
}

#[test]
fn test_capacity() {
    let mut map = LinearMap::new();
    map.insert(1, 2);
    assert!(map.capacity() >= 1);
    map.remove(&1);
    assert!(map.capacity() >= 1);
    map.reserve(TEST_CAPACITY);
    let capacity = map.capacity();
    assert!(capacity >= TEST_CAPACITY);
    for i in 0..TEST_CAPACITY as i32 {
        assert!(map.insert(i, i).is_none());
    }
    assert_eq!(capacity, map.capacity());
}

#[test]
fn test_reserve() {
    let mut map = LinearMap::new();
    map.reserve(TEST_CAPACITY);
    assert!(map.capacity() >= TEST_CAPACITY);
    for i in 0..TEST_CAPACITY as i32 {
        assert!(map.insert(i, i).is_none());
    }
    map.reserve(TEST_CAPACITY);
    assert!(map.capacity() >= 2 * TEST_CAPACITY);

    let mut map = LinearMap::new();
    map.reserve(TEST_CAPACITY);
    assert!(map.capacity() >= TEST_CAPACITY);
    for i in 0..TEST_CAPACITY as i32 {
        assert!(map.insert(i, i).is_none());
    }
    map.reserve(TEST_CAPACITY);
    assert!(map.capacity() >= 2 * TEST_CAPACITY);
}

#[test]
fn test_shrink_to_fit() {
    let mut map = LinearMap::new();
    map.shrink_to_fit();
    assert_eq!(map.capacity(), 0);
    map.reserve(TEST_CAPACITY);
    map.shrink_to_fit();
    assert_eq!(map.capacity(), 0);
    for i in 0..TEST_CAPACITY as i32 {
        assert!(map.insert(i, i).is_none());
    }
    map.shrink_to_fit();
    assert_eq!(map.len(), TEST_CAPACITY);
    assert!(map.capacity() >= TEST_CAPACITY);
}

#[test]
fn test_len_and_is_empty() {
    let mut map = LinearMap::new();
    assert_eq!(map.len(), 0);
    assert!(map.is_empty());
    map.insert(100, 100);
    assert_eq!(map.len(), 1);
    assert!(!map.is_empty());
    for i in 0..TEST_CAPACITY as i32 {
        assert!(map.insert(i, i).is_none());
    }
    assert_eq!(map.len(), 1 + TEST_CAPACITY);
    assert!(!map.is_empty());
    assert!(map.remove(&100).is_some());
    assert_eq!(map.len(), TEST_CAPACITY);
    assert!(!map.is_empty());
}

#[test]
fn test_clear() {
    let mut map = LinearMap::new();
    map.clear();
    assert_eq!(map.len(), 0);
    for i in 0..TEST_CAPACITY as i32 {
        assert!(map.insert(i, i).is_none());
    }
    map.clear();
    assert_eq!(map.len(), 0);
    assert!(map.capacity() > 0);
}

#[test]
fn test_iterators() {
    const ONE:   i32 = 0b0001;
    const TWO:   i32 = 0b0010;
    const THREE: i32 = 0b0100;
    const FOUR:  i32 = 0b1000;
    const ALL:   i32 = 0b1111;
    let mut map = LinearMap::new();
    assert!(map.insert(ONE, TWO).is_none());
    assert!(map.insert(TWO, THREE).is_none());
    assert!(map.insert(THREE, FOUR).is_none());
    assert!(map.insert(FOUR, ONE).is_none());

    {
        let mut result_k = 0;
        let mut result_v = 0;
        for (&k, &v) in map.iter() {
            result_k ^= k;
            result_v ^= v;
            assert_eq!(((k << 1) & ALL) | ((k >> 3) & ALL), v);
        }
        assert_eq!(result_k, ALL);
        assert_eq!(result_v, ALL);
    }
    {
        let mut result_k = 0;
        let mut result_v = 0;
        for (&k, &mut v) in map.iter_mut() {
            result_k ^= k;
            result_v ^= v;
            assert_eq!(((k << 1) & ALL) | ((k >> 3) & ALL), v);
        }
        assert_eq!(result_k, ALL);
        assert_eq!(result_v, ALL);
    }
    {
        let mut result = 0;
        for &k in map.keys() {
            result ^= k;
        }
        assert_eq!(result, ALL);
    }
    {
        let mut result = 0;
        for &v in map.values() {
            result ^= v;
        }
        assert_eq!(result, ALL);
    }
}

#[test]
fn test_insert_remove_get() {
    let mut map = LinearMap::new();
    assert!(map.insert(100, 101).is_none());
    assert!(map.contains_key(&100));
    assert_eq!(map.get(&100), Some(&101));
    assert_eq!(map.get_mut(&100), Some(&mut 101));
    for i in 0..TEST_CAPACITY as i32 {
        assert!(map.insert(i, i).is_none());
    }
    assert_eq!(map.insert(100, 102), Some(101));
    assert_eq!(map.remove(&100), Some(102));
    assert_eq!(map.remove(&100), None);
    assert_eq!(map.remove(&1000), None);
}

#[test]
fn test_entry() {
    let xs = [(1, 10), (2, 20), (3, 30), (4, 40), (5, 50), (6, 60)];

    let mut map = LinearMap::new();

    for &(k, v) in &xs {
        map.insert(k, v);
    }

    // Existing key (insert)
    match map.entry(1) {
        Vacant(_) => unreachable!(),
        Occupied(mut view) => {
            assert_eq!(view.get(), &10);
            assert_eq!(view.insert(100), 10);
        }
    }
    assert_eq!(map.get(&1).unwrap(), &100);
    assert_eq!(map.len(), 6);


    // Existing key (update)
    match map.entry(2) {
        Vacant(_) => unreachable!(),
        Occupied(mut view) => {
            let v = view.get_mut();
            let new_v = (*v) * 10;
            *v = new_v;
        }
    }
    assert_eq!(map.get(&2).unwrap(), &200);
    assert_eq!(map.len(), 6);

    // Existing key (take)
    match map.entry(3) {
        Vacant(_) => unreachable!(),
        Occupied(view) => {
            assert_eq!(view.remove(), 30);
        }
    }
    assert_eq!(map.get(&3), None);
    assert_eq!(map.len(), 5);


    // Inexistent key (insert)
    match map.entry(10) {
        Occupied(_) => unreachable!(),
        Vacant(view) => {
            assert_eq!(*view.insert(1000), 1000);
        }
    }
    assert_eq!(map.get(&10).unwrap(), &1000);
    assert_eq!(map.len(), 6);
}

#[test]
fn test_eq() {
    let kvs = vec![('a', 1), ('b', 2), ('c', 3)];

    let mut m1: LinearMap<_, _> = kvs.clone().into_iter().collect();
    let m2: LinearMap<_, _> = kvs.into_iter().rev().collect();
    assert_eq!(m1, m2);

    m1.insert('a', 11);
    assert!(m1 != m2);

    m1.insert('a', 1);
    assert_eq!(m1, m2);

    m1.remove(&'a');
    assert!(m1 != m2);
}

#[test]
fn test_macro() {
    let names = linear_map!{
        1 => "one",
        2 => "two",
    };
    assert_eq!(names.len(), 2);
    assert_eq!(names.capacity(), 2);
    assert_eq!(names[&1], "one");
    assert_eq!(names[&2], "two");
    assert_eq!(names.get(&3), None);

    let empty: LinearMap<i32, i32> = linear_map!{};
    assert_eq!(empty.len(), 0);
    assert_eq!(empty.capacity(), 0);

    let _nested_compiles = linear_map!{
        1 => linear_map!{0 => 1 + 2,},
        2 => linear_map!{1 => 1,},
    };
}

#[test]
fn test_retain() {
    let mut map: LinearMap<isize, isize> = (0..100).map(|x|(x, x*10)).collect();
    map.retain(|&k, _| k % 2 == 0);
    assert_eq!(map.len(), 50);
    assert_eq!(map[&2], 20);
    assert_eq!(map[&4], 40);
    assert_eq!(map[&6], 60);
}
