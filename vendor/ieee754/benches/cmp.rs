#[macro_use] extern crate criterion;
extern crate ieee754;

use criterion::{Criterion, Fun, black_box};

use std::{f32, f64};
use ieee754::Ieee754;
use std::cmp::Ordering;

fn f32_cmp(c: &mut Criterion) {
    let std = Fun::new("std", |b, &data: &&[f32]| {
        b.iter(|| {
            for x in data {
                for y in data {
                    if black_box(x).partial_cmp(y).unwrap_or(Ordering::Less) == Ordering::Less {
                        black_box(*x);
                    }
                }
            }
        })
    });
    let total = Fun::new("total", |b, &data: &&[f32]| {
        b.iter(|| {
            for x in data {
                for y in data {
                    if black_box(x).total_cmp(y) == Ordering::Less {
                        black_box(*x);
                    }
                }
            }
        })
    });

    const DATA: &[f32] = &[
        -f32::INFINITY,-0.97,0.66,-0.00,-0.49,1.68,f32::NAN,0.13,-0.41,0.59,
        -0.47,0.00,1.03,0.89,1.97,f32::INFINITY,0.17,-0.30,-0.16,-f32::NAN,
    ];
    c.bench_functions("f32_cmp", vec![std, total], DATA);
}

fn f64_cmp(c: &mut Criterion) {
    let std = Fun::new("std", |b, &data: &&[f64]| {
        b.iter(|| {
            for x in data {
                for y in data {
                    if black_box(x).partial_cmp(y).unwrap_or(Ordering::Less) == Ordering::Less {
                        black_box(*x);
                    }
                }
            }
        })
    });
    let total = Fun::new("total", |b, &data: &&[f64]| {
        b.iter(|| {
            for x in data {
                for y in data {
                    if black_box(x).total_cmp(y) == Ordering::Less {
                        black_box(*x);
                    }
                }
            }
        })
    });

    const DATA: &[f64] = &[
        -f64::INFINITY,-0.97,0.66,-0.00,-0.49,1.68,f64::NAN,0.13,-0.41,0.59,
        -0.47,0.00,1.03,0.89,1.97,f64::INFINITY,0.17,-0.30,-0.16,-f64::NAN,
    ];
    c.bench_functions("f64_cmp", vec![std, total], DATA);
}

fn f32_sort(c: &mut Criterion) {
    let baseline = Fun::new("baseline", |b, &data: &&[f32]| {
        let v = data.iter().cloned().cycle().take(1000).collect::<Vec<_>>();

        b.iter(|| v.clone());
    });
    let std = Fun::new("std", |b, &data: &&[f32]| {
        let v = data.iter().cloned().cycle().take(1000).collect::<Vec<_>>();

        b.iter(|| {
            v.clone().sort_by(
                |a, b| a.partial_cmp(b).unwrap_or_else(|| {
                    if a.is_nan() { Ordering::Less } else { Ordering::Greater }
                }))
        })
    });
    let total = Fun::new("total", |b, &data: &&[f32]| {
        let v = data.iter().cloned().cycle().take(1000).collect::<Vec<_>>();

        b.iter(|| {
            v.clone().sort_by(|a, b| a.total_cmp(b))
        })
    });

    const DATA: &[f32] = &[
        -f32::INFINITY,-0.97,0.66,-0.00,-0.49,1.68,f32::NAN,0.13,-0.41,0.59,
        -0.47,0.00,1.03,0.89,1.97,f32::INFINITY,0.17,-0.30,-0.16,-f32::NAN,
    ];
    c.bench_functions("f32_sort", vec![baseline, std, total], DATA);
}

fn f64_sort(c: &mut Criterion) {
    let baseline = Fun::new("baseline", |b, &data: &&[f64]| {
        let v = data.iter().cloned().cycle().take(1000).collect::<Vec<_>>();

        b.iter(|| v.clone());
    });
    let std = Fun::new("std", |b, &data: &&[f64]| {
        let v = data.iter().cloned().cycle().take(1000).collect::<Vec<_>>();

        b.iter(|| {
            v.clone().sort_by(
                |a, b| a.partial_cmp(b).unwrap_or_else(|| {
                    if a.is_nan() { Ordering::Less } else { Ordering::Greater }
                }))
        })
    });
    let total = Fun::new("total", |b, &data: &&[f64]| {
        let v = data.iter().cloned().cycle().take(1000).collect::<Vec<_>>();

        b.iter(|| {
            v.clone().sort_by(|a, b| a.total_cmp(b))
        })
    });

    const DATA: &[f64] = &[
        -f64::INFINITY,-0.97,0.66,-0.00,-0.49,1.68,f64::NAN,0.13,-0.41,0.59,
        -0.47,0.00,1.03,0.89,1.97,f64::INFINITY,0.17,-0.30,-0.16,-f64::NAN,
    ];
    c.bench_functions("f64_sort", vec![baseline, std, total], DATA);
}

criterion_group!(benches, f32_cmp, f64_cmp, f32_sort, f64_sort);
criterion_main!(benches);
