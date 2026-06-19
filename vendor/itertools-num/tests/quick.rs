#[macro_use]
extern crate quickcheck;
extern crate itertools_num;

fn exact_size<I: ExactSizeIterator>(mut it: I) -> bool {
    // check every iteration
    let (mut low, mut hi) = it.size_hint();
    if Some(low) != hi { return false; }
    while let Some(_) = it.next() {
        let (xlow, xhi) = it.size_hint();
        if low != xlow + 1 { return false; }
        low = xlow;
        hi = xhi;
        if Some(low) != hi { return false; }
    }
    let (low, hi) = it.size_hint();
    low == 0 && hi == Some(0)
}


quickcheck! {
    fn size_linspace(a: f32, b: f32, n: usize) -> bool {
        let it = itertools_num::linspace(a, b, n);
        it.len() == n &&
            exact_size(it)
    }
}

