extern crate criterion;
extern crate ieee754;
use ieee754::Ieee754;
use std::f32 as f;

fn main() {
    let count = f::NEG_INFINITY.upto(f::INFINITY).map(criterion::black_box).count();
    let expected =
        /* bits */ (1u64 << 32)
        - /* NaNs */ (1 << 24)
        - /* -0.0 */ 1
        + /* infinities */ 2;

    assert_eq!(count, expected as usize);
    println!("there are {} non-NaN floats", count);
}
