extern crate criterion;

// base line comparison for the all example.
fn main() {
    let expected =
        /* bits */ (1u64 << 32)
        - /* NaNs */ (1 << 24)
        - /* -0.0 */ 1
        + /* infinities */ 2;

    println!("count {}", (0..expected as usize).map(criterion::black_box).count());
}
