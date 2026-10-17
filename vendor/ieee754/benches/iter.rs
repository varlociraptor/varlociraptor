#[macro_use] extern crate criterion;
extern crate ieee754;

use criterion::{Criterion, Fun, black_box};

use ieee754::Ieee754;

fn f32_iter(c: &mut Criterion) {
    type B = <f32 as Ieee754>::Bits;
    type S = <f32 as Ieee754>::Significand;

    let positive = Fun::new("positive", |b, &i: &usize| {
        let end = f32::recompose(false, 0, i as S);
        b.iter(|| {
               assert_eq!(black_box(1_f32).upto(black_box(end))
                          .map(black_box).count(),
                          i + 1);
        })
    });
    let over_zero = Fun::new("over_zero", |b, &i: &usize| {
        let x = f32::recompose(false, -f32::exponent_bias(), (i / 2) as S);
        b.iter(|| {
               assert_eq!(black_box(-x).upto(black_box(x))
                          .map(black_box).count(),
                          i + 1);
        })
    });
    let baseline = Fun::new("baseline", |b, &i: &usize| {
        b.iter(|| {
            assert_eq!((black_box(0 as B)..=black_box(i as B))
                       .map(black_box).count(),
                       i + 1);
        })
    });

    c.bench_functions("f32_iter", vec![positive, over_zero, baseline], 40);
}

fn f64_iter(c: &mut Criterion) {
    type B = <f64 as Ieee754>::Bits;
    type S = <f64 as Ieee754>::Significand;

    let positive = Fun::new("positive", |b, &i: &usize| {
        let end = f64::recompose(false, 0, i as S);
        b.iter(|| {
               assert_eq!(black_box(1_f64).upto(black_box(end))
                          .map(black_box).count(),
                          i + 1);
        })
    });
    let over_zero = Fun::new("over_zero", |b, &i: &usize| {
        let x = f64::recompose(false, -f64::exponent_bias(), (i / 2) as S);
        b.iter(|| {
               assert_eq!(black_box(-x).upto(black_box(x))
                          .map(black_box).count(),
                          i + 1);
        })
    });
    let baseline = Fun::new("baseline", |b, &i: &usize| {
        b.iter(|| {
            assert_eq!((black_box(0 as B)..=black_box(i as B))
                       .map(black_box).count(),
                       i + 1);
        })
    });

    c.bench_functions("f64_iter", vec![positive, over_zero, baseline], 40);
}

criterion_group!(benches, f32_iter, f64_iter);
criterion_main!(benches);
