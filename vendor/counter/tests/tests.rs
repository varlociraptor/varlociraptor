#[cfg(test)]
mod tests {
    use counter::Counter;
    use rand::Rng;

    #[test]
    fn test_composite_add_sub() {
        let mut counts = "able babble table babble rabble table able fable scrabble"
            .split_whitespace()
            .collect::<Counter<_>>();
        // add or subtract an iterable of the same type
        counts += "cain and abel fable table cable".split_whitespace();
        // or add or subtract from another Counter of the same type
        let other_counts = "scrabble cabbie fable babble"
            .split_whitespace()
            .collect::<Counter<_>>();
        let _diff = counts - other_counts;
    }

    #[test]
    fn test_most_common() {
        let counter = "abbccc".chars().collect::<Counter<_>>();
        let by_common = counter.most_common();
        let expected = vec![('c', 3), ('b', 2), ('a', 1)];
        assert!(by_common == expected);
    }

    #[test]
    fn test_most_common_tiebreaker() {
        let counter = "eaddbbccc".chars().collect::<Counter<_>>();
        let by_common = counter.most_common_tiebreaker(|&a, &b| a.cmp(&b));
        let expected = vec![('c', 3), ('b', 2), ('d', 2), ('a', 1), ('e', 1)];
        assert!(by_common == expected);
    }

    #[test]
    fn test_most_common_tiebreaker_reversed() {
        let counter = "eaddbbccc".chars().collect::<Counter<_>>();
        let by_common = counter.most_common_tiebreaker(|&a, &b| b.cmp(&a));
        let expected = vec![('c', 3), ('d', 2), ('b', 2), ('e', 1), ('a', 1)];
        assert!(by_common == expected);
    }

    // The main purpose of this test is to see that we can call `Counter::most_common_tiebreaker()`
    // with a closure that is `FnMut` but not `Fn`.
    #[test]
    fn test_most_common_tiebreaker_fn_mut() {
        let counter: Counter<_> = "abracadabra".chars().collect::<Counter<_>>();
        // Count how many times the tiebreaker closure is called.
        let mut num_ties = 0;
        let sorted = counter.most_common_tiebreaker(|a, b| {
            num_ties += 1;
            a.cmp(b)
        });
        let expected = vec![('a', 5), ('b', 2), ('r', 2), ('c', 1), ('d', 1)];
        assert_eq!(sorted, expected);
        // We should have called the tiebreaker twice: once to resolve the tie between `'b'` and
        // `'r'` and once to resolve the tie between `'c'` and `'d'`.
        assert_eq!(num_ties, 2);
    }

    #[test]
    fn test_most_common_ordered() {
        let counter = "eaddbbccc".chars().collect::<Counter<_>>();
        let by_common = counter.most_common_ordered();
        let expected = vec![('c', 3), ('b', 2), ('d', 2), ('a', 1), ('e', 1)];
        assert!(by_common == expected);
    }

    #[test]
    fn test_k_most_common_ordered() {
        let counter: Counter<_> = "abracadabra".chars().collect();
        let all = counter.most_common_ordered();
        for k in 0..=counter.len() {
            let topk = counter.k_most_common_ordered(k);
            assert_eq!(&topk, &all[..k]);
        }
    }

    /// This test is fundamentally the same as `test_k_most_common_ordered`, but it operates on
    /// a wider variety of data. In particular, it tests both longer, narrower, and wider
    /// distributions of data than the other test does.
    #[test]
    fn test_k_most_common_ordered_heavy() {
        let mut rng = rand::thread_rng();

        for container_size in [5, 10, 25, 100, 256] {
            for max_value_factor in [0.25, 0.5, 1.0, 1.25, 2.0, 10.0, 100.0] {
                let max_value = ((container_size as f64) * max_value_factor) as u32;
                let mut values = vec![0; container_size];
                for value in values.iter_mut() {
                    *value = rng.gen_range(0..=max_value);
                }

                let counter: Counter<_> = values.into_iter().collect();
                let all = counter.most_common_ordered();
                for k in 0..=counter.len() {
                    let topk = counter.k_most_common_ordered(k);
                    assert_eq!(&topk, &all[..k]);
                }
            }
        }
    }

    #[test]
    fn test_total() {
        let counter = "".chars().collect::<Counter<_>>();
        let total: usize = counter.total();
        assert_eq!(total, 0);

        let counter = "eaddbbccc".chars().collect::<Counter<_>>();
        let total: usize = counter.total();
        assert_eq!(total, 9);
    }

    #[test]
    fn test_add() {
        let d = "abbccc".chars().collect::<Counter<_>>();
        let e = "bccddd".chars().collect::<Counter<_>>();

        let out = d + e;
        let expected = "abbbcccccddd".chars().collect::<Counter<_>>();
        assert!(out == expected);
    }

    #[test]
    fn test_sub() {
        let d = "abbccc".chars().collect::<Counter<_>>();
        let e = "bccddd".chars().collect::<Counter<_>>();

        let out = d - e;
        let expected = "abc".chars().collect::<Counter<_>>();
        assert!(out == expected);
    }

    #[test]
    fn test_intersection() {
        let d = "abbccc".chars().collect::<Counter<_>>();
        let e = "bccddd".chars().collect::<Counter<_>>();

        let out = d & e;
        let expected = "bcc".chars().collect::<Counter<_>>();
        assert!(out == expected);
    }

    #[test]
    fn test_union() {
        let d = "abbccc".chars().collect::<Counter<_>>();
        let e = "bccddd".chars().collect::<Counter<_>>();

        let out = d | e;
        let expected = "abbcccddd".chars().collect::<Counter<_>>();
        assert!(out == expected);
    }

    #[test]
    fn test_delete_key_from_backing_map() {
        let mut counter = "aa-bb-cc".chars().collect::<Counter<_>>();
        counter.remove(&'-');
        assert!(counter == "aabbcc".chars().collect::<Counter<_>>());
    }

    #[test]
    fn test_superset_non_usize_count() {
        let mut a: Counter<_, i8> = "abbcccc".chars().collect();
        let mut b: Counter<_, i8> = "abb".chars().collect();
        assert!(a.is_superset(&b));
        // Negative values are possible, a is no longer a superset
        a[&'e'] = -1;
        assert!(!a.is_superset(&b));
        // Adjust b to make a a superset again
        b[&'e'] = -2;
        assert!(a.is_superset(&b));
    }

    #[test]
    fn test_subset_non_usize_count() {
        let mut a: Counter<_, i8> = "abb".chars().collect();
        let mut b: Counter<_, i8> = "abbcccc".chars().collect();
        assert!(a.is_subset(&b));
        // Negative values are possible; a is no longer a subset
        b[&'e'] = -1;
        assert!(!a.is_subset(&b));
        // Adjust a to make it a subset again
        a[&'e'] = -2;
        assert!(a.is_subset(&b));
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_serialize_deserialize() {
        let a = "abbccc".chars().collect::<Counter<_>>();
        let serialized = serde_json::to_string(&a).unwrap();
        let b: Counter<char> = serde_json::from_str(&serialized).unwrap();
        assert!(a == b)
    }
}
