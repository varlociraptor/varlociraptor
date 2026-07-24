extern crate itertools_num;
extern crate itertools;

use itertools_num::linspace;
use itertools::assert_equal;

#[test]
fn test_linspace() {
    let iter = linspace::<f32>(0., 2., 3);
    assert_equal(iter, vec![0., 1., 2.]);

    let iter = linspace::<f32>(0., 1.0, 11);
    for (a, b) in iter.zip(vec![ 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]) {
        assert!((a - b).abs() < 1.0e-6);
    }

    let iter = linspace::<f32>(0., 100.0, 1001);
    assert_eq!(iter.last(), Some(100.0));

    let iter = linspace::<f32>(0., -2., 4);
    assert_equal(iter, vec![0., -0.666666666667, -1.333333333333, -2.]);

    let iter = linspace::<f32>(0., -2., 4);
    assert_equal(iter.rev(), vec![-2., -1.333333333333, -0.666666666667, 0.]);

    let iter = linspace::<f32>(0., 1., 1);
    assert_equal(iter, vec![0.]);

    let iter = linspace::<f32>(0., 1., 1);
    assert_equal(iter.rev(), vec![0.]);

    let mut iter = linspace::<f32>(0., 1., 0);
    assert_eq!(iter.next(), None);
}

