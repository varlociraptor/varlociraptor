use crate::Counter;
use maplit::hashmap;
use std::collections::HashMap;
#[test]
fn test_creation() {
    let _: Counter<usize> = Counter::new();

    let counter = Counter::from_iter(&[1]);

    let mut expected = HashMap::new();
    static ONE: usize = 1;
    expected.insert(&ONE, 1);
    assert!(counter.map == expected);
}

#[test]
fn test_creation_with_capacity() {
    let counter: Counter<usize, usize> = Counter::with_capacity(3);
    assert_eq!(counter.map.capacity(), 3);
}

#[test]
fn test_update() {
    let mut counter = "abbccc".chars().collect::<Counter<_>>();
    let expected = hashmap! {
        'a' => 1,
        'b' => 2,
        'c' => 3,
    };
    assert!(counter.map == expected);

    counter.update("aeeeee".chars());
    let expected = hashmap! {
        'a' => 2,
        'b' => 2,
        'c' => 3,
        'e' => 5,
    };
    assert!(counter.map == expected);
}

#[test]
fn test_add_update_iterable() {
    let mut counter = "abbccc".chars().collect::<Counter<_>>();
    let expected = hashmap! {
        'a' => 1,
        'b' => 2,
        'c' => 3,
    };
    assert!(counter.map == expected);

    counter += "aeeeee".chars();
    let expected = hashmap! {
        'a' => 2,
        'b' => 2,
        'c' => 3,
        'e' => 5,
    };
    assert!(counter.map == expected);
}

#[test]
fn test_add_update_counter() {
    let mut counter = "abbccc".chars().collect::<Counter<_>>();
    let expected = hashmap! {
        'a' => 1,
        'b' => 2,
        'c' => 3,
    };
    assert!(counter.map == expected);

    let other = "aeeeee".chars().collect::<Counter<_>>();
    counter += other;
    let expected = hashmap! {
        'a' => 2,
        'b' => 2,
        'c' => 3,
        'e' => 5,
    };
    assert!(counter.map == expected);
}

#[test]
fn test_subtract() {
    let mut counter = "abbccc".chars().collect::<Counter<_>>();
    counter.subtract("bbccddd".chars());
    let expected = hashmap! {
        'a' => 1,
        'c' => 1,
    };
    assert!(counter.map == expected);
}

#[test]
fn test_sub_update_iterable() {
    let mut counter = "abbccc".chars().collect::<Counter<_>>();
    counter -= "bbccddd".chars();
    let expected = hashmap! {
        'a' => 1,
        'c' => 1,
    };
    assert!(counter.map == expected);
}

#[test]
fn test_sub_update_counter() {
    let mut counter = "abbccc".chars().collect::<Counter<_>>();
    let other = "bbccddd".chars().collect::<Counter<_>>();
    counter -= other;
    let expected = hashmap! {
        'a' => 1,
        'c' => 1,
    };
    assert!(counter.map == expected);
}

#[test]
fn test_from_iter_simple() {
    let counter = "abbccc".chars().collect::<Counter<_>>();
    let expected = hashmap! {
        'a' => 1,
        'b' => 2,
        'c' => 3,
    };
    assert!(counter.map == expected);
}

#[test]
fn test_from_iter_tuple() {
    let items = [('a', 1), ('b', 2), ('c', 3)];
    let counter = items.iter().cloned().collect::<Counter<_>>();
    let expected: HashMap<char, usize> = items.iter().cloned().collect();
    assert_eq!(counter.map, expected);
}

#[test]
fn test_from_iter_tuple_with_duplicates() {
    let items = [('a', 1), ('b', 2), ('c', 3)];
    let counter = items
        .iter()
        .cycle()
        .take(items.len() * 2)
        .cloned()
        .collect::<Counter<_>>();
    let expected: HashMap<char, usize> = items.iter().map(|(c, n)| (*c, n * 2)).collect();
    assert_eq!(counter.map, expected);
}

#[test]
fn test_extend_simple() {
    let mut counter = "abbccc".chars().collect::<Counter<_>>();
    counter.extend("bccddd".chars());
    let expected = hashmap! {
        'a' => 1,
        'b' => 3,
        'c' => 5,
        'd' => 3,
    };
    assert!(counter.map == expected);
}

#[test]
fn test_extend_tuple() {
    let mut counter = "bccddd".chars().collect::<Counter<_>>();
    let items = [('a', 1), ('b', 2), ('c', 3)];
    counter.extend(items.iter().cloned());
    let expected = hashmap! {
        'a' => 1,
        'b' => 3,
        'c' => 5,
        'd' => 3,
    };
    assert_eq!(counter.map, expected);
}

#[test]
fn test_extend_tuple_with_duplicates() {
    let mut counter = "ccc".chars().collect::<Counter<_>>();
    let items = [('a', 1), ('b', 2), ('c', 3)];
    counter.extend(items.iter().cycle().take(items.len() * 2 - 1).cloned());
    let expected: HashMap<char, usize> = items.iter().map(|(c, n)| (*c, n * 2)).collect();
    assert_eq!(counter.map, expected);
}

#[test]
fn test_count_minimal_type() {
    #[derive(Debug, Hash, PartialEq, Eq)]
    struct Inty {
        i: usize,
    }

    impl Inty {
        pub fn new(i: usize) -> Inty {
            Inty { i }
        }
    }

    // <https://en.wikipedia.org/wiki/867-5309/Jenny>
    let intys = vec![
        Inty::new(8),
        Inty::new(0),
        Inty::new(0),
        Inty::new(8),
        Inty::new(6),
        Inty::new(7),
        Inty::new(5),
        Inty::new(3),
        Inty::new(0),
        Inty::new(9),
    ];

    let inty_counts = intys.into_iter().collect::<Counter<_>>();
    // println!("{:?}", inty_counts.map); // test runner blanks this
    // {Inty { i: 8 }: 2, Inty { i: 0 }: 3, Inty { i: 9 }: 1, Inty { i: 3 }: 1,
    //  Inty { i: 7 }: 1, Inty { i: 6 }: 1, Inty { i: 5 }: 1}
    assert!(inty_counts.map.get(&Inty { i: 8 }) == Some(&2));
    assert!(inty_counts.map.get(&Inty { i: 0 }) == Some(&3));
    assert!(inty_counts.map.get(&Inty { i: 6 }) == Some(&1));
}

#[test]
fn test_collect() {
    let counter: Counter<_> = "abbccc".chars().collect();
    let expected = hashmap! {
        'a' => 1,
        'b' => 2,
        'c' => 3,
    };
    assert!(counter.map == expected);
}

#[test]
fn test_non_usize_count() {
    let counter: Counter<_, i8> = "abbccc".chars().collect();
    let expected = hashmap! {
        'a' => 1,
        'b' => 2,
        'c' => 3,
    };
    assert!(counter.map == expected);
}
