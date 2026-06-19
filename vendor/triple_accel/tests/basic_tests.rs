use triple_accel::*;
use triple_accel::hamming::*;
use triple_accel::levenshtein::*;

#[test]
fn test_basic_hamming_naive() {
    let a1 = b"abc";
    let b1 = b"abd";
    let dist = hamming_naive(a1, b1);
    assert!(dist == 1);

    let a2 = b"";
    let b2 = b"";
    let dist = hamming_naive(a2, b2);
    assert!(dist == 0);
}

#[test]
fn test_basic_hamming_search_naive() {
    let a1 = b"abc";
    let b1 = b"  abc  abb";
    let mut res: Vec<Match> = hamming_search_naive_with_opts(a1, b1, 1, SearchType::All).collect();
    assert!(res == vec![Match{start: 2, end: 5, k: 0}, Match{start: 7, end: 10, k: 1}]);

    let a2 = b"abc";
    let b2 = b"  abc  abb";
    res = hamming_search_naive(a2, b2).collect();
    assert!(res == vec![Match{start: 2, end: 5, k: 0}]);
}

#[test]
fn test_basic_hamming_search_simd() {
    let a1 = b"abc";
    let b1 = b"  abc  abb aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
    let mut res: Vec<Match> = hamming_search_simd_with_opts(a1, b1, 1, SearchType::All).collect();
    assert!(res == vec![Match{start: 2, end: 5, k: 0}, Match{start: 7, end: 10, k: 1}]);

    let a2 = b"abc";
    let b2 = b"  abc  abb aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
    res = hamming_search_simd(a2, b2).collect();
    assert!(res == vec![Match{start: 2, end: 5, k: 0}]);
}

#[test]
fn test_basic_hamming_words_64() {
    let a_str = b"abc";
    let b_str = b"abd";

    let mut a = alloc_str(a_str.len());
    fill_str(&mut a, a_str);

    let mut b = alloc_str(b_str.len());
    fill_str(&mut b, b_str);

    let dist = hamming_words_64(&a, &b);
    assert!(dist == 1);
}

#[test]
fn test_basic_hamming_words_128() {
    let a_str = b"abc";
    let b_str = b"abd";

    let mut a = alloc_str(a_str.len());
    fill_str(&mut a, a_str);

    let mut b = alloc_str(b_str.len());
    fill_str(&mut b, b_str);

    let dist = hamming_words_128(&a, &b);
    assert!(dist == 1);
}

#[test]
fn test_basic_hamming_simd_movemask() {
    let a1 = b"abcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
    let b1 = b"abdaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
    let dist = hamming_simd_movemask(a1, b1);
    assert!(dist == 1);

    let a2 = b"";
    let b2 = b"";
    let dist = hamming_simd_movemask(a2, b2);
    assert!(dist == 0);
}

#[test]
fn test_basic_hamming_simd_parallel() {
    let a1 = b"abcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
    let b1 = b"abdaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
    let dist = hamming_simd_parallel(a1, b1);
    assert!(dist == 1);

    let a2 = b"";
    let b2 = b"";
    let dist = hamming_simd_parallel(a2, b2);
    assert!(dist == 0);
}

#[test]
fn test_basic_levenshtein_naive() {
    let a1 = b"abcde";
    let b1 = b" ab cde";
    let mut res = levenshtein_naive(a1, b1);
    assert!(res == 2);

    let a2 = b"abcde";
    let b2 = b"";
    res = levenshtein_naive(a2, b2);
    assert!(res == 5);

    let a3 = b"abcde";
    let b3 = b"abcdee";
    res = levenshtein_naive(a3, b3);
    assert!(res == 1);

    let a4 = b"abcde";
    let b4 = b"acde";
    res = levenshtein_naive(a4, b4);
    assert!(res == 1);

    let a5 = b"abcde";
    let b5 = b"abbde";
    res = levenshtein_naive(a5, b5);
    assert!(res == 1);

    let a6 = b"abcde";
    let b6 = b"acbde";
    res = levenshtein_naive_with_opts(a6, b6, false, EditCosts::new(1, 1, 0, Some(1))).0;
    assert!(res == 1);

    let a7 = b"ab";
    let b7 = b"ba";
    res = levenshtein_naive_with_opts(a7, b7, false, EditCosts::new(1, 1, 0, Some(1))).0;
    assert!(res == 1);

    let a8 = b"abc";
    let b8 = b"aac";
    res = levenshtein_naive_with_opts(a8, b8, false, EditCosts::new(2, 3, 0, None)).0;
    assert!(res == 2);

    let a9 = b"abc";
    let b9 = b"aac";
    res = levenshtein_naive_with_opts(a9, b9, false, EditCosts::new(3, 1, 0, None)).0;
    assert!(res == 2);

    let a10 = b"abc";
    let b10 = b"ac";
    res = levenshtein_naive_with_opts(a10, b10, false, EditCosts::new(1, 1, 2, None)).0;
    assert!(res == 3);

    let a11 = b"acde";
    let b11 = b"abce";
    res = levenshtein_naive_with_opts(a11, b11, false, EditCosts::new(2, 1, 2, None)).0;
    assert!(res == 4);

    let a12 = b"abcde";
    let b12 = b"abe";
    res = levenshtein_naive_with_opts(a12, b12, false, EditCosts::new(1, 1, 2, None)).0;
    assert!(res == 4);
}

#[test]
fn test_trace_on_levenshtein_naive() {
    let a1 = b"abcde";
    let b1 = b" ab cde";
    let mut res = levenshtein_naive_with_opts(a1, b1, true, LEVENSHTEIN_COSTS);
    assert!(res.0 == 2);
    assert!(res.1.unwrap() == vec![Edit{edit: EditType::AGap, count: 1},
                                   Edit{edit: EditType::Match, count: 2},
                                   Edit{edit: EditType::AGap, count: 1},
                                   Edit{edit: EditType::Match, count: 3}]);

    let a2 = b"abcde";
    let b2 = b"";
    res = levenshtein_naive_with_opts(a2, b2, true, LEVENSHTEIN_COSTS);
    assert!(res.0 == 5);
    assert!(res.1.unwrap() == vec![Edit{edit: EditType::BGap, count: 5}]);

    let a3 = b"abcde";
    let b3 = b"abcce";
    res = levenshtein_naive_with_opts(a3, b3, true, LEVENSHTEIN_COSTS);
    assert!(res.0 == 1);
    assert!(res.1.unwrap() == vec![Edit{edit: EditType::Match, count: 3},
                                   Edit{edit: EditType::Mismatch, count: 1},
                                   Edit{edit: EditType::Match, count: 1}]);

    let a4 = b"abcde";
    let b4 = b"acbde";
    res = levenshtein_naive_with_opts(a4, b4, true, EditCosts::new(1, 1, 0, Some(1)));
    assert!(res.0 == 1);
    assert!(res.1.unwrap() == vec![Edit{edit: EditType::Match, count: 1},
                                   Edit{edit: EditType::Transpose, count: 1},
                                   Edit{edit: EditType::Match, count: 2}]);
}

#[test]
fn test_basic_levenshtein() {
    let a1 = b"abcde";
    let b1 = b" ab cde";
    let mut res = levenshtein(a1, b1);
    assert!(res == 2);

    let a2 = b"abcde";
    let b2 = b"";
    res = levenshtein(a2, b2);
    assert!(res == 5);

    let a3 = b"abcde";
    let b3 = b"abcdee";
    res = levenshtein(a3, b3);
    assert!(res == 1);

    let a4 = b"abcde";
    let b4 = b"acde";
    res = levenshtein(a4, b4);
    assert!(res == 1);

    let a5 = b"abcde";
    let b5 = b"abbde";
    res = levenshtein(a5, b5);
    assert!(res == 1);
}

#[test]
fn test_basic_levenshtein_exp() {
    let a1 = b"abcde";
    let b1 = b" ab cde";
    let mut res = levenshtein_exp(a1, b1);
    assert!(res == 2);

    let a2 = b"abcde";
    let b2 = b"";
    res = levenshtein_exp(a2, b2);
    assert!(res == 5);

    let a3 = b"abcde";
    let b3 = b"abcdee";
    res = levenshtein_exp(a3, b3);
    assert!(res == 1);

    let a4 = b"abcde";
    let b4 = b"acde";
    res = levenshtein_exp(a4, b4);
    assert!(res == 1);

    let a5 = b"abcde";
    let b5 = b"abbde";
    res = levenshtein_exp(a5, b5);
    assert!(res == 1);
}

#[test]
fn test_basic_rdamerau() {
    let a1 = b"abcde";
    let b1 = b" ab dce";
    let mut res = rdamerau(a1, b1);
    assert!(res == 3);

    let a2 = b"abcde";
    let b2 = b"";
    res = rdamerau(a2, b2);
    assert!(res == 5);

    let a3 = b"abcde";
    let b3 = b"bacdee";
    res = rdamerau(a3, b3);
    assert!(res == 2);

    let a4 = b"abcde";
    let b4 = b"acde";
    res = rdamerau(a4, b4);
    assert!(res == 1);

    let a5 = b"abcde";
    let b5 = b"abbde";
    res = rdamerau(a5, b5);
    assert!(res == 1);
}

#[test]
fn test_basic_rdamerau_exp() {
    let a1 = b"abcde";
    let b1 = b" ab dce";
    let mut res = rdamerau_exp(a1, b1);
    assert!(res == 3);

    let a2 = b"abcde";
    let b2 = b"";
    res = rdamerau_exp(a2, b2);
    assert!(res == 5);

    let a3 = b"abcde";
    let b3 = b"bacdee";
    res = rdamerau_exp(a3, b3);
    assert!(res == 2);

    let a4 = b"abcde";
    let b4 = b"acde";
    res = rdamerau_exp(a4, b4);
    assert!(res == 1);

    let a5 = b"abcde";
    let b5 = b"abbde";
    res = rdamerau_exp(a5, b5);
    assert!(res == 1);
}

#[test]
fn test_basic_levenshtein_naive_k_with_opts() {
    let a1 = b"abcde";
    let b1 = b" ab cde";
    let mut res = levenshtein_naive_k_with_opts(a1, b1, 2, false, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 2);
    assert!(res.1.is_none());

    let a2 = b"abcde";
    let b2 = b"";
    res = levenshtein_naive_k_with_opts(a2, b2, 10, false, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 5);
    assert!(res.1.is_none());

    let a3 = b"abcde";
    let b3 = b"abcdee";
    res = levenshtein_naive_k_with_opts(a3, b3, 2, false, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.is_none());

    let a4 = b"abcde";
    let b4 = b"acde";
    res = levenshtein_naive_k_with_opts(a4, b4, 2, false, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.is_none());

    let a5 = b"abcde";
    let b5 = b"abbde";
    res = levenshtein_naive_k_with_opts(a5, b5, 2, false, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.is_none());

    let a6 = b"abcde";
    let b6 = b"abbde";
    res = levenshtein_naive_k_with_opts(a6, b6, 1, false, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.is_none());

    let a7 = b"abcde";
    let b7 = b"acbde";
    res = levenshtein_naive_k_with_opts(a7, b7, 1, false, EditCosts::new(1, 1, 0, Some(1))).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.is_none());

    let a8 = b"ab";
    let b8 = b"ba";
    res = levenshtein_naive_k_with_opts(a8, b8, 1, false, EditCosts::new(1, 1, 0, Some(1))).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.is_none());

    let a9 = b"abc";
    let b9 = b"aac";
    res = levenshtein_naive_k_with_opts(a9, b9, 5, false, EditCosts::new(2, 3, 0, None)).unwrap();
    assert!(res.0 == 2);
    assert!(res.1.is_none());

    let a10 = b"abc";
    let b10 = b"aac";
    res = levenshtein_naive_k_with_opts(a10, b10, 5, false, EditCosts::new(3, 1, 0, None)).unwrap();
    assert!(res.0 == 2);
    assert!(res.1.is_none());

    let a11 = b"abc";
    let b11 = b"ac";
    res = levenshtein_naive_k_with_opts(a11, b11, 5, false, EditCosts::new(1, 1, 2, None)).unwrap();
    assert!(res.0 == 3);
    assert!(res.1.is_none());

    let a12 = b"acde";
    let b12 = b"abce";
    res = levenshtein_naive_k_with_opts(a12, b12, 5, false, EditCosts::new(2, 1, 2, None)).unwrap();
    assert!(res.0 == 4);
    assert!(res.1.is_none());

    let a13 = b"abcde";
    let b13 = b"abe";
    res = levenshtein_naive_k_with_opts(a13, b13, 5, false, EditCosts::new(1, 1, 2, None)).unwrap();
    assert!(res.0 == 4);
    assert!(res.1.is_none());

    let a14 = b"abcde";
    let b14 = b"hello";
    let res1 = levenshtein_naive_k_with_opts(a14, b14, 1, false, RDAMERAU_COSTS);
    assert!(res1.is_none());
}

#[test]
fn test_trace_on_levenshtein_naive_k_with_opts() {
    let a1 = b"abcde";
    let b1 = b" ab cde";
    let mut res = levenshtein_naive_k_with_opts(a1, b1, 2, true, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 2);
    assert!(res.1.unwrap() == vec![Edit{edit: EditType::AGap, count: 1},
                                   Edit{edit: EditType::Match, count: 2},
                                   Edit{edit: EditType::AGap, count: 1},
                                   Edit{edit: EditType::Match, count: 3}]);

    let a2 = b"abcde";
    let b2 = b"";
    res = levenshtein_naive_k_with_opts(a2, b2, 10, true, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 5);
    assert!(res.1.unwrap() == vec![Edit{edit: EditType::BGap, count: 5}]);

    let a3 = b"abcde";
    let b3 = b"abcce";
    res = levenshtein_naive_k_with_opts(a3, b3, 2, true, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.unwrap() == vec![Edit{edit: EditType::Match, count: 3},
                                   Edit{edit: EditType::Mismatch, count: 1},
                                   Edit{edit: EditType::Match, count: 1}]);

    let a4 = b"abcde";
    let b4 = b"acbde";
    res = levenshtein_naive_k_with_opts(a4, b4, 2, true, EditCosts::new(1, 1, 0, Some(1))).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.unwrap() == vec![Edit{edit: EditType::Match, count: 1},
                                   Edit{edit: EditType::Transpose, count: 1},
                                   Edit{edit: EditType::Match, count: 2}]);
}

#[test]
fn test_basic_levenshtein_simd_k_with_opts() {
    let a1 = b"abcde";
    let b1 = b" ab cde";
    let mut res = levenshtein_simd_k_with_opts(a1, b1, 2, false, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 2);
    assert!(res.1.is_none());

    let a2 = b"abcde";
    let b2 = b"";
    res = levenshtein_simd_k_with_opts(a2, b2, 30, false, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 5);
    assert!(res.1.is_none());

    let a3 = b"abcde";
    let b3 = b"abcdee";
    res = levenshtein_simd_k_with_opts(a3, b3, 20, false, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.is_none());

    let a4 = b"abcde";
    let b4 = b"acde";
    res = levenshtein_simd_k_with_opts(a4, b4, 1, false, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.is_none());

    let a5 = b"abcde";
    let b5 = b"abbde";
    res = levenshtein_simd_k_with_opts(a5, b5, 2, false, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.is_none());

    let a6 = b"abcde";
    let b6 = b"acbde";
    res = levenshtein_simd_k_with_opts(a6, b6, 2, false, EditCosts::new(1, 1, 0, Some(1))).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.is_none());

    let a7 = b"ab";
    let b7 = b"ba";
    res = levenshtein_simd_k_with_opts(a7, b7, 2, false, EditCosts::new(1, 1, 0, Some(1))).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.is_none());

    let a8 = b"abc";
    let b8 = b"aac";
    res = levenshtein_simd_k_with_opts(a8, b8, 5, false, EditCosts::new(2, 3, 0, None)).unwrap();
    assert!(res.0 == 2);
    assert!(res.1.is_none());

    let a9 = b"abc";
    let b9 = b"aac";
    res = levenshtein_simd_k_with_opts(a9, b9, 5, false, EditCosts::new(3, 1, 0, None)).unwrap();
    assert!(res.0 == 2);
    assert!(res.1.is_none());

    let a10 = b"abc";
    let b10 = b"ac";
    res = levenshtein_simd_k_with_opts(a10, b10, 5, false, EditCosts::new(1, 1, 2, None)).unwrap();
    assert!(res.0 == 3);
    assert!(res.1.is_none());

    let a11 = b"acde";
    let b11 = b"abce";
    res = levenshtein_simd_k_with_opts(a11, b11, 5, false, EditCosts::new(2, 1, 2, None)).unwrap();
    assert!(res.0 == 4);
    assert!(res.1.is_none());

    let a12 = b"abcde";
    let b12 = b"abe";
    res = levenshtein_simd_k_with_opts(a12, b12, 5, false, EditCosts::new(1, 1, 2, None)).unwrap();
    assert!(res.0 == 4);
    assert!(res.1.is_none());

    let a13 = b"\0";
    let b13 = b"";
    res = levenshtein_simd_k_with_opts(a13, b13, 2, false, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.is_none());

    let a14 = b"ab\0de";
    let b14 = b"a\0bde";
    res = levenshtein_simd_k_with_opts(a14, b14, 2, false, EditCosts::new(1, 1, 0, Some(1))).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.is_none());

    let a15 = b"\0b";
    let b15 = b"b\0";
    res = levenshtein_simd_k_with_opts(a15, b15, 2, false, EditCosts::new(1, 1, 0, Some(1))).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.is_none());

    let a16 = b"\0";
    let b16 = b"\0\0";
    res = levenshtein_simd_k_with_opts(a16, b16, 2, false, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.is_none());

    let a17 = b"\0";
    let b17 = b"\0";
    res = levenshtein_simd_k_with_opts(a17, b17, 2, false, EditCosts::new(1, 1, 0, Some(1))).unwrap();
    assert!(res.0 == 0);
    assert!(res.1.is_none());

    let a18 = b"\0\0b\0";
    let b18 = b"\0b\0\0";
    res = levenshtein_simd_k_with_opts(a18, b18, 2, false, EditCosts::new(1, 1, 0, Some(1))).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.is_none());

    let a19 = b"abcde";
    let b19 = b"hello";
    let res1 = levenshtein_simd_k_with_opts(a19, b19, 1, false, RDAMERAU_COSTS);
    assert!(res1.is_none());
}

#[test]
fn test_trace_on_levenshtein_simd_k_with_opts() {
    let a1 = b"abcde";
    let b1 = b" ab cde";
    let mut res = levenshtein_simd_k_with_opts(a1, b1, 30, true, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 2);
    assert!(res.1.unwrap() == vec![Edit{edit: EditType::AGap, count: 1},
                                   Edit{edit: EditType::Match, count: 2},
                                   Edit{edit: EditType::AGap, count: 1},
                                   Edit{edit: EditType::Match, count: 3}]);

    let a2 = b"abcde";
    let b2 = b"";
    res = levenshtein_simd_k_with_opts(a2, b2, 5, true, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 5);
    assert!(res.1.unwrap() == vec![Edit{edit: EditType::BGap, count: 5}]);

    let a3 = b"abcde";
    let b3 = b"abcce";
    res = levenshtein_simd_k_with_opts(a3, b3, 1, true, LEVENSHTEIN_COSTS).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.unwrap() == vec![Edit{edit: EditType::Match, count: 3},
                                   Edit{edit: EditType::Mismatch, count: 1},
                                   Edit{edit: EditType::Match, count: 1}]);

    let a4 = b"abcde";
    let b4 = b"acbde";
    res = levenshtein_simd_k_with_opts(a4, b4, 2, true, EditCosts::new(1, 1, 0, Some(1))).unwrap();
    assert!(res.0 == 1);
    assert!(res.1.unwrap() == vec![Edit{edit: EditType::Match, count: 1},
                                   Edit{edit: EditType::Transpose, count: 1},
                                   Edit{edit: EditType::Match, count: 2}]);
}

#[test]
fn test_basic_levenshtein_search_naive() {
    let a1 = b"bcc";
    let b1 = b"abcde";
    let k1 = 1;
    let mut res: Vec<Match> = levenshtein_search_naive_with_opts(a1, b1, k1, SearchType::All, LEVENSHTEIN_COSTS, false).collect();
    assert!(res == vec![Match{start: 1, end: 3, k: 1}, Match{start: 1, end: 4, k: 1}]);

    let a2 = b"";
    let b2 = b"";
    let k2 = 1;
    res = levenshtein_search_naive_with_opts(a2, b2, k2, SearchType::All, LEVENSHTEIN_COSTS, false).collect();
    assert!(res == vec![]);

    let a3 = b"tast";
    let b3 = b"testing 123 tating!";
    let k3 = 1;
    res = levenshtein_search_naive_with_opts(a3, b3, k3, SearchType::All, LEVENSHTEIN_COSTS, false).collect();
    assert!(res == vec![Match{start: 0, end: 4, k: 1}, Match{start: 12, end: 15, k: 1}]);

    let a4 = b"tst";
    let b4 = b"testing 123 tasting!";
    let k4 = 1;
    res = levenshtein_search_naive_with_opts(a4, b4, k4, SearchType::All, LEVENSHTEIN_COSTS, false).collect();
    assert!(res == vec![Match{start: 0, end: 4, k: 1}, Match{start: 12, end: 16, k: 1}]);

    let a5 = b"tst";
    let b5 = b"testing 123 tasting!";
    res = levenshtein_search_naive(a5, b5).collect();
    assert!(res == vec![Match{start: 0, end: 4, k: 1}, Match{start: 12, end: 16, k: 1}]);

    let a6 = b"ab";
    let b6 = b"ba";
    let k6 = 1;
    res = levenshtein_search_naive_with_opts(a6, b6, k6, SearchType::All, EditCosts::new(1, 1, 0, Some(1)), false).collect();
    assert!(res == vec![Match{start: 0, end: 1, k: 1}, Match{start: 0, end: 2, k: 1}]);

    let a7 = b"test";
    let b7 = b"...tseting!";
    let k7 = 1;
    res = levenshtein_search_naive_with_opts(a7, b7, k7, SearchType::All, EditCosts::new(1, 1, 0, Some(1)), false).collect();
    assert!(res == vec![Match{start: 3, end: 7, k: 1}]);

    let a8 = b"test";
    let b8 = b"...tssting!";
    let k8 = 2;
    res = levenshtein_search_naive_with_opts(a8, b8, k8, SearchType::All, EditCosts::new(3, 1, 0, None), false).collect();
    assert!(res == vec![Match{start: 3, end: 5, k: 2}, Match{start: 3, end: 7, k: 2}]);

    let a9 = b"tst";
    let b9 = b"testing 123 tasting";
    let k9 = 1;
    let res1 = levenshtein_search_naive_with_opts(a9, b9, k9, SearchType::All, LEVENSHTEIN_COSTS, false).next().unwrap();
    assert!(res1 == Match{start: 0, end: 4, k: 1});

    let a10 = b"test";
    let b10 = b" testing 123 tasting";
    let k10 = 1;
    res = levenshtein_search_naive_with_opts(a10, b10, k10, SearchType::All, LEVENSHTEIN_COSTS, true).collect();
    assert!(res == vec![Match{start: 1, end: 5, k: 1}]);

    let a11 = b"test";
    let b11 = b" etsting 123 tasting";
    let k11 = 2;
    res = levenshtein_search_naive_with_opts(a11, b11, k11, SearchType::All, RDAMERAU_COSTS, true).collect();
    assert!(res == vec![Match{start: 0, end: 3, k: 2}, Match{start: 0, end: 4, k: 2}, Match{start: 1, end: 5, k: 2}]);

    let a12 = b"test";
    let b12 = b"etsting";
    let k12 = 1;
    res = levenshtein_search_naive_with_opts(a12, b12, k12, SearchType::All, RDAMERAU_COSTS, true).collect();
    assert!(res == vec![Match{start: 0, end: 4, k: 1}]);

    let a13 = b"test";
    let b13 = b"est";
    let k13 = 3;
    res = levenshtein_search_naive_with_opts(a13, b13, k13, SearchType::All, EditCosts::new(1, 1, 2, None), true).collect();
    assert!(res == vec![Match{start: 0, end: 3, k: 3}]);

    let a14 = b"testing";
    let b14 = b"   teing";
    let k14 = 4;
    res = levenshtein_search_naive_with_opts(a14, b14, k14, SearchType::All, EditCosts::new(1, 1, 2, None), false).collect();
    assert!(res == vec![Match{start: 1, end: 8, k: 4}]);

    let a15 = b"testing";
    let b15 = b"   teing";
    let k15 = 4;
    res = levenshtein_search_naive_with_opts(a15, b15, k15, SearchType::All, EditCosts::new(2, 1, 2, None), false).collect();
    assert!(res == vec![Match{start: 3, end: 8, k: 4}]);

    let a16 = b"abc";
    let b16 = b"";
    let k16 = 5;
    res = levenshtein_search_naive_with_opts(a16, b16, k16, SearchType::All, LEVENSHTEIN_COSTS, false).collect();
    assert!(res == vec![Match{start: 0, end: 0, k: 3}]);

    let a17 = b"";
    let b17 = b"abc";
    let k17 = 2;
    res = levenshtein_search_naive_with_opts(a17, b17, k17, SearchType::All, LEVENSHTEIN_COSTS, true).collect();
    assert!(res == vec![Match{start: 0, end: 0, k: 0}, Match{start: 0, end: 1, k: 1}, Match{start: 0, end: 2, k: 2}]);
}

#[test]
fn test_basic_levenshtein_search_simd() {
    let a1 = b"bcc";
    let b1 = b"abcde";
    let k1 = 1;
    let mut res: Vec<Match> = levenshtein_search_simd_with_opts(a1, b1, k1, SearchType::All, LEVENSHTEIN_COSTS, false).collect();
    assert!(res == vec![Match{start: 1, end: 3, k: 1}, Match{start: 1, end: 4, k: 1}]);

    let a2 = b"";
    let b2 = b"";
    let k2 = 1;
    res = levenshtein_search_simd_with_opts(a2, b2, k2, SearchType::All, LEVENSHTEIN_COSTS, false).collect();
    assert!(res == vec![]);

    let a3 = b"tast";
    let b3 = b"testing 123 tating!";
    let k3 = 1;
    res = levenshtein_search_simd_with_opts(a3, b3, k3, SearchType::All, LEVENSHTEIN_COSTS, false).collect();
    assert!(res == vec![Match{start: 0, end: 4, k: 1}, Match{start: 12, end: 15, k: 1}]);

    let a4 = b"tst";
    let b4 = b"testing 123 tasting!";
    let k4 = 1;
    res = levenshtein_search_simd_with_opts(a4, b4, k4, SearchType::All, LEVENSHTEIN_COSTS, false).collect();
    assert!(res == vec![Match{start: 0, end: 4, k: 1}, Match{start: 12, end: 16, k: 1}]);

    let a5 = b"tst";
    let b5 = b"testing 123 tasting!";
    res = levenshtein_search_simd(a5, b5).collect();
    assert!(res == vec![Match{start: 0, end: 4, k: 1}, Match{start: 12, end: 16, k: 1}]);

    let a6 = b"ab";
    let b6 = b"ba";
    let k6 = 1;
    res = levenshtein_search_simd_with_opts(a6, b6, k6, SearchType::All, EditCosts::new(1, 1, 0, Some(1)), false).collect();
    assert!(res == vec![Match{start: 0, end: 1, k: 1}, Match{start: 0, end: 2, k: 1}]);

    let a7 = b"test";
    let b7 = b"...tseting!";
    let k7 = 1;
    res = levenshtein_search_simd_with_opts(a7, b7, k7, SearchType::All, EditCosts::new(1, 1, 0, Some(1)), false).collect();
    assert!(res == vec![Match{start: 3, end: 7, k: 1}]);

    let a8 = b"test";
    let b8 = b"...tssting!";
    let k8 = 2;
    res = levenshtein_search_simd_with_opts(a8, b8, k8, SearchType::All, EditCosts::new(3, 1, 0, None), false).collect();
    assert!(res == vec![Match{start: 3, end: 5, k: 2}, Match{start: 3, end: 7, k: 2}]);

    let a9 = b"tst";
    let b9 = b"testing 123 tasting";
    let k9 = 1;
    let res1 = levenshtein_search_simd_with_opts(a9, b9, k9, SearchType::All, LEVENSHTEIN_COSTS, false).next().unwrap();
    assert!(res1 == Match{start: 0, end: 4, k: 1});

    let a10 = b"test";
    let b10 = b" testing 123 tasting";
    let k10 = 1;
    res = levenshtein_search_simd_with_opts(a10, b10, k10, SearchType::All, LEVENSHTEIN_COSTS, true).collect();
    assert!(res == vec![Match{start: 1, end: 5, k: 1}]);

    let a11 = b"test";
    let b11 = b" etsting 123 tasting";
    let k11 = 2;
    res = levenshtein_search_simd_with_opts(a11, b11, k11, SearchType::All, RDAMERAU_COSTS, true).collect();
    assert!(res == vec![Match{start: 0, end: 3, k: 2}, Match{start: 0, end: 4, k: 2}, Match{start: 1, end: 5, k: 2}]);

    let a12 = b"test";
    let b12 = b"etsting";
    let k12 = 1;
    res = levenshtein_search_simd_with_opts(a12, b12, k12, SearchType::All, RDAMERAU_COSTS, true).collect();
    assert!(res == vec![Match{start: 0, end: 4, k: 1}]);

    let a13 = b"test";
    let b13 = b"est";
    let k13 = 3;
    res = levenshtein_search_simd_with_opts(a13, b13, k13, SearchType::All, EditCosts::new(1, 1, 2, None), true).collect();
    assert!(res == vec![Match{start: 0, end: 3, k: 3}]);

    let a14 = b"testing";
    let b14 = b"   teing";
    let k14 = 4;
    res = levenshtein_search_simd_with_opts(a14, b14, k14, SearchType::All, EditCosts::new(1, 1, 2, None), false).collect();
    assert!(res == vec![Match{start: 1, end: 8, k: 4}]);

    let a15 = b"testing";
    let b15 = b"   teing";
    let k15 = 4;
    res = levenshtein_search_simd_with_opts(a15, b15, k15, SearchType::All, EditCosts::new(2, 1, 2, None), false).collect();
    assert!(res == vec![Match{start: 3, end: 8, k: 4}]);

    let a16 = b"\0b";
    let b16 = b"b\0";
    let k16 = 1;
    res = levenshtein_search_simd_with_opts(a16, b16, k16, SearchType::All, RDAMERAU_COSTS, true).collect();
    assert!(res == vec![Match{start: 0, end: 1, k: 1}, Match{start: 0, end: 2, k: 1}]);

    let a17 = b"\0\0";
    let b17 = b"\0\0";
    let k17 = 0;
    res = levenshtein_search_simd_with_opts(a17, b17, k17, SearchType::All, RDAMERAU_COSTS, true).collect();
    assert!(res == vec![Match{start: 0, end: 2, k: 0}]);

    let a18 = b"testing";
    let b18 = b"   \0esting";
    let k18 = 1;
    res = levenshtein_search_simd_with_opts(a18, b18, k18, SearchType::All, LEVENSHTEIN_COSTS, false).collect();
    assert!(res == vec![Match{start: 3, end: 10, k: 1}]);

    let a19 = b"\0\0\0";
    let b19 = b"\0\0";
    let k19 = 1;
    res = levenshtein_search_simd_with_opts(a19, b19, k19, SearchType::All, LEVENSHTEIN_COSTS, true).collect();
    assert!(res == vec![Match{start: 0, end: 2, k: 1}]);

    let a20 = b"\0\0";
    let b20 = b"   \0\0";
    let k20 = 0;
    res = levenshtein_search_simd_with_opts(a20, b20, k20, SearchType::All, RDAMERAU_COSTS, false).collect();
    assert!(res == vec![Match{start: 3, end: 5, k: 0}]);

    let a21 = b"abc";
    let b21 = b"";
    let k21 = 5;
    res = levenshtein_search_simd_with_opts(a21, b21, k21, SearchType::All, LEVENSHTEIN_COSTS, false).collect();
    assert!(res == vec![Match{start: 0, end: 0, k: 3}]);

    let a22 = b"";
    let b22 = b"abc";
    let k22 = 2;
    res = levenshtein_search_simd_with_opts(a22, b22, k22, SearchType::All, LEVENSHTEIN_COSTS, true).collect();
    assert!(res == vec![Match{start: 0, end: 0, k: 0}, Match{start: 0, end: 1, k: 1}, Match{start: 0, end: 2, k: 2}]);
}

