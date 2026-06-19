use editdistancek::{edit_distance, edit_distance_bounded, mismatch};
use levenshtein::levenshtein;
use rand::RngCore;
#[allow(unused_imports)]
use std::cmp::{max, min};

#[test]
fn test_equal_strings() {
    let tests = vec!["", "abacaba", "a"];
    for s in tests.iter() {
        assert_eq!(edit_distance(s.as_bytes(), s.as_bytes()), 0);
    }
}

#[test]
fn naive() {
    let tests = vec![
        ("kitten", "sitting"),
        ("book", "back"),
        ("table", "dinner"),
        ("person", "pardon"),
        ("person", "persons"),
        ("", "aba"),
        ("aba", ""),
    ];
    for (s, t) in tests.iter() {
        assert_eq!(
            edit_distance(s.as_bytes(), t.as_bytes()),
            levenshtein(&s, &t)
        );
    }
}

#[test]
fn test_basic_example() {
    let tests = vec![
        ("kitten", "sitting", 3),
        ("book", "back", 2),
        ("table", "dinner", 5),
        ("person", "pardon", 2),
        ("person", "persons", 1),
        ("", "aba", 3),
        ("aba", "", 3),
    ];
    for test in tests.iter() {
        assert_eq!(
            Some(test.2),
            edit_distance_bounded(test.0.as_bytes(), test.1.as_bytes(), test.2)
        );
        assert_eq!(test.2, edit_distance(test.0.as_bytes(), test.1.as_bytes()));
    }
}

#[test]
fn test_substitution() {
    assert_eq!(edit_distance("abacaba".as_bytes(), "abadabc".as_bytes()), 2);
    assert_eq!(
        edit_distance_bounded("abacaba".as_bytes(), "abadabc".as_bytes(), 0),
        None
    );
    assert_eq!(
        edit_distance_bounded("abacaba".as_bytes(), "abadabc".as_bytes(), 1),
        None
    );
    assert_eq!(
        edit_distance_bounded("abacaba".as_bytes(), "abadabc".as_bytes(), 2),
        Some(2)
    );
    assert_eq!(
        edit_distance_bounded("abacaba".as_bytes(), "abadabc".as_bytes(), 3),
        Some(2)
    );
}

#[test]
fn test_insert() {
    assert_eq!(edit_distance("abacab".as_bytes(), "abacaba".as_bytes()), 1);
}

#[test]
fn test_remove() {
    assert_eq!(edit_distance("abacaba".as_bytes(), "abacab".as_bytes()), 1);
}

#[test]
fn test_mismatch_128() {
    let mut s = [0u8; 128];
    rand::thread_rng().fill_bytes(&mut s);
    let mut t = s.clone();
    for i in 0..128 {
        t[i] = s[i] ^ 1;
        assert_eq!(mismatch(&s, &t), i);
        t[i] = s[i];
    }
}

#[test]
fn test_mismatch() {
    for l in 0..256 {
        let mut s = vec![0u8; l];
        rand::thread_rng().fill_bytes(&mut s);
        let mut t = s.clone();
        for i in 0..l {
            t[i] = s[i] ^ 1;
            assert_eq!(mismatch(&s, &t), i);
            t[i] = s[i];
        }
    }
}

extern crate quickcheck;

use quickcheck::quickcheck;

quickcheck! {
    fn equal(s : String) -> bool {
        edit_distance(s.as_bytes(), s.as_bytes()) == 0
    }

    fn symetry(s : String, t : String) -> bool {
        edit_distance(s.as_bytes(), t.as_bytes()) == edit_distance(t.as_bytes(), s.as_bytes())
    }

    fn triangle(a : String, b : String, c : String) -> bool {
        edit_distance(a.as_bytes(), b.as_bytes()) <=
        edit_distance(a.as_bytes(), c.as_bytes()) + edit_distance(c.as_bytes(), b.as_bytes())
    }

    fn size_difference_bound(s : String, t : String) -> bool {
        let lens = s.len();
        let lent = t.len();
        let diff = if lens > lent {lens - lent} else {lent - lens};
        edit_distance(s.as_bytes(), t.as_bytes()) >= diff
    }

    fn size_difference_bound_k(s : String, t : String) -> bool {
        let lens = s.chars().count();
        let lent = t.chars().count();
        let diff = if lens > lent {lens - lent} else {lent - lens};
        if let Some(x) = edit_distance_bounded(s.as_bytes(), t.as_bytes(), diff) {
           x >= diff
        } else {
            true
        }
    }

    fn upper_bound(s : String, t : String) -> bool {
        let lens = s.len();
        let lent = t.len();
        edit_distance(s.as_bytes(), t.as_bytes()) <= max(lens, lent)
    }

    fn upper_bound_k(s : String, t : String) -> bool {
        let lens = s.len();
        let lent = t.len();
        edit_distance_bounded(s.as_bytes(), t.as_bytes(), max(lens, lent)) <= Some(max(lens, lent))
    }

    fn bounded_max(s : String, t : String) -> bool {
        let lens = s.len();
        let lent = t.len();
        Some(edit_distance(s.as_bytes(), t.as_bytes())) == edit_distance_bounded(s.as_bytes(), t.as_bytes(), max(lens, lent))
    }

}
