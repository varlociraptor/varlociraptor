extern crate linear_map;

use linear_map::set::LinearSet;

#[test]
fn test_disjoint() {
    let mut xs = LinearSet::new();
    let mut ys = LinearSet::new();
    assert!(xs.is_disjoint(&ys));
    assert!(ys.is_disjoint(&xs));
    assert!(xs.insert(5));
    assert!(ys.insert(11));
    assert!(xs.is_disjoint(&ys));
    assert!(ys.is_disjoint(&xs));
    assert!(xs.insert(7));
    assert!(xs.insert(19));
    assert!(xs.insert(4));
    assert!(ys.insert(2));
    assert!(ys.insert(-11));
    assert!(xs.is_disjoint(&ys));
    assert!(ys.is_disjoint(&xs));
    assert!(ys.insert(7));
    assert!(!xs.is_disjoint(&ys));
    assert!(!ys.is_disjoint(&xs));
}

#[test]
fn test_subset_and_superset() {
    let mut a = LinearSet::new();
    assert!(a.insert(0));
    assert!(a.insert(5));
    assert!(a.insert(11));
    assert!(a.insert(7));

    let mut b = LinearSet::new();
    assert!(b.insert(0));
    assert!(b.insert(7));
    assert!(b.insert(19));
    assert!(b.insert(250));
    assert!(b.insert(11));
    assert!(b.insert(200));

    assert!(!a.is_subset(&b));
    assert!(!a.is_superset(&b));
    assert!(!b.is_subset(&a));
    assert!(!b.is_superset(&a));

    assert!(b.insert(5));

    assert!(a.is_subset(&b));
    assert!(!a.is_superset(&b));
    assert!(!b.is_subset(&a));
    assert!(b.is_superset(&a));
}

#[test]
fn test_iterate() {
    let mut a = LinearSet::new();
    for i in 0..32 {
        assert!(a.insert(i));
    }
    let mut observed: u32 = 0;
    for k in &a {
        observed |= 1 << *k;
    }
    assert_eq!(observed, 0xFFFF_FFFF);
}

#[test]
fn test_intersection() {
    let mut a = LinearSet::new();
    let mut b = LinearSet::new();

    assert!(a.insert(11));
    assert!(a.insert(1));
    assert!(a.insert(3));
    assert!(a.insert(77));
    assert!(a.insert(103));
    assert!(a.insert(5));
    assert!(a.insert(-5));

    assert!(b.insert(2));
    assert!(b.insert(11));
    assert!(b.insert(77));
    assert!(b.insert(-9));
    assert!(b.insert(-42));
    assert!(b.insert(5));
    assert!(b.insert(3));

    let mut i = 0;
    let expected = [3, 5, 11, 77];
    for x in a.intersection(&b) {
        assert!(expected.contains(x));
        i += 1
    }
    assert_eq!(i, expected.len());
}

#[test]
fn test_difference() {
    let mut a = LinearSet::new();
    let mut b = LinearSet::new();

    assert!(a.insert(1));
    assert!(a.insert(3));
    assert!(a.insert(5));
    assert!(a.insert(9));
    assert!(a.insert(11));

    assert!(b.insert(3));
    assert!(b.insert(9));

    let mut i = 0;
    let expected = [1, 5, 11];
    for x in a.difference(&b) {
        assert!(expected.contains(x));
        i += 1
    }
    assert_eq!(i, expected.len());
}

#[test]
fn test_symmetric_difference() {
    let mut a = LinearSet::new();
    let mut b = LinearSet::new();

    assert!(a.insert(1));
    assert!(a.insert(3));
    assert!(a.insert(5));
    assert!(a.insert(9));
    assert!(a.insert(11));

    assert!(b.insert(-2));
    assert!(b.insert(3));
    assert!(b.insert(9));
    assert!(b.insert(14));
    assert!(b.insert(22));

    let mut i = 0;
    let expected = [-2, 1, 5, 11, 14, 22];
    for x in a.symmetric_difference(&b) {
        assert!(expected.contains(x));
        i += 1
    }
    assert_eq!(i, expected.len());
}

#[test]
fn test_union() {
    let mut a = LinearSet::new();
    let mut b = LinearSet::new();

    assert!(a.insert(1));
    assert!(a.insert(3));
    assert!(a.insert(5));
    assert!(a.insert(9));
    assert!(a.insert(11));
    assert!(a.insert(16));
    assert!(a.insert(19));
    assert!(a.insert(24));

    assert!(b.insert(-2));
    assert!(b.insert(1));
    assert!(b.insert(5));
    assert!(b.insert(9));
    assert!(b.insert(13));
    assert!(b.insert(19));

    let mut i = 0;
    let expected = [-2, 1, 3, 5, 9, 11, 13, 16, 19, 24];
    for x in a.union(&b) {
        assert!(expected.contains(x));
        i += 1
    }
    assert_eq!(i, expected.len());
}

#[test]
fn test_from_iter() {
    let xs = [1, 2, 3, 4, 5, 6, 7, 8, 9];

    let set: LinearSet<_> = xs.iter().cloned().collect();

    for x in &xs {
        assert!(set.contains(x));
    }
}

#[test]
fn test_move_iter() {
    let hs = {
        let mut hs = LinearSet::new();

        hs.insert('a');
        hs.insert('b');

        hs
    };

    let v = hs.into_iter().collect::<Vec<char>>();
    assert!(v == ['a', 'b'] || v == ['b', 'a']);
}

#[test]
fn test_eq() {
    // These constants once happened to expose a bug in insert().
    // I'm keeping them around to prevent a regression.
    let mut s1 = LinearSet::new();

    s1.insert(1);
    s1.insert(2);
    s1.insert(3);

    let mut s2 = LinearSet::new();

    s2.insert(1);
    s2.insert(2);

    assert!(s1 != s2);

    s2.insert(3);

    assert_eq!(s1, s2);
}

#[test]
fn test_show() {
    let mut set = LinearSet::new();
    let empty = LinearSet::<i32>::new();

    set.insert(1);
    set.insert(2);

    let set_str = format!("{:?}", set);

    assert!(set_str == "{1, 2}" || set_str == "{2, 1}");
    assert_eq!(format!("{:?}", empty), "{}");
}

#[test]
fn test_trivial_drain() {
    let mut s = LinearSet::<i32>::new();
    for _ in s.drain() {}
    assert!(s.is_empty());
    drop(s);

    let mut s = LinearSet::<i32>::new();
    drop(s.drain());
    assert!(s.is_empty());
}

#[test]
fn test_drain() {
    let mut s: LinearSet<_> = (1..100).collect();

    // try this a bunch of times to make sure we don't screw up internal state.
    for _ in 0..20 {
        assert_eq!(s.len(), 99);

        {
            let mut last_i = 0;
            let mut d = s.drain();
            for (i, x) in d.by_ref().take(50).enumerate() {
                last_i = i;
                assert!(x != 0);
            }
            assert_eq!(last_i, 49);
        }

        for _ in &s { panic!("s should be empty!"); }

        // reset to try again.
        s.extend(1..100);
    }
}

#[test]
fn test_extend_ref() {
    let mut a = LinearSet::new();
    a.insert(1);

    a.extend(&[2, 3, 4]);

    assert_eq!(a.len(), 4);
    assert!(a.contains(&1));
    assert!(a.contains(&2));
    assert!(a.contains(&3));
    assert!(a.contains(&4));

    let mut b = LinearSet::new();
    b.insert(5);
    b.insert(6);

    a.extend(&b);

    assert_eq!(a.len(), 6);
    assert!(a.contains(&1));
    assert!(a.contains(&2));
    assert!(a.contains(&3));
    assert!(a.contains(&4));
    assert!(a.contains(&5));
    assert!(a.contains(&6));
}

#[test]
fn test_into_vec() {
    let mut a = LinearSet::new();
    assert!(a.insert(0));
    assert!(a.insert(5));
    assert!(a.insert(11));
    assert!(a.insert(7));
    let set_as_vec: Vec<_> = a.into();
    assert_eq!(set_as_vec, vec![0, 5, 11, 7]);
}

#[test]
fn test_retain() {
    let xs = [1,2,3,4,5,6];
    let mut set: LinearSet<isize> = xs.iter().cloned().collect();
    set.retain(|&k| k % 2 == 0);
    assert_eq!(set.len(), 3);
    assert!(set.contains(&2));
    assert!(set.contains(&4));
    assert!(set.contains(&6));
}
