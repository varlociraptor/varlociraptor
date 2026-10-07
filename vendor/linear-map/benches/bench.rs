#![feature(test)]

extern crate linear_map;
extern crate test;

use linear_map::LinearMap;

const SMALL:  u32 =   10;
const MEDIUM: u32 =  100;
const BIG:    u32 = 1000;

fn insert(b: &mut test::Bencher, num: u32) {
    b.iter(|| {
        let mut map = LinearMap::new();
        for i in 0..num {
            map.insert(i, i);
        }
    })
}

fn remove_insert(b: &mut test::Bencher, num: u32) {
    b.iter(|| {
        let mut map = LinearMap::new();
        for i in 0..num {
            map.insert(i, i);
        }
        for i in 0..num {
            map.remove(&i);
        }
    })
}

fn remove_rev_insert(b: &mut test::Bencher, num: u32) {
    b.iter(|| {
        let mut map = LinearMap::new();
        for i in 0..num {
            map.insert(i, i);
        }
        for i in 0..num {
            map.remove(&(num - i - 1));
        }
    })
}

fn get_middle(b: &mut test::Bencher, num: u32) {
    let mut map = LinearMap::new();
    for i in 0..num {
        map.insert(i, i);
    }
    let middle = num / 2;
    b.iter(|| {
        test::black_box(map.get(&middle));
        test::black_box(map.get_mut(&middle));
    })
}

fn get_none(b: &mut test::Bencher, num: u32) {
    let mut map = LinearMap::new();
    for i in 0..num {
        map.insert(i, i);
    }
    let none = num + 1;
    b.iter(|| {
        test::black_box(map.get(&none));
        test::black_box(map.get_mut(&none));
    })
}

#[bench] fn bench_insert_small (b: &mut test::Bencher) { insert(b, SMALL);  }
#[bench] fn bench_insert_medium(b: &mut test::Bencher) { insert(b, MEDIUM); }
#[bench] fn bench_insert_big   (b: &mut test::Bencher) { insert(b, BIG);    }

#[bench] fn bench_remove_insert_small (b: &mut test::Bencher) { remove_insert(b, SMALL);  }
#[bench] fn bench_remove_insert_medium(b: &mut test::Bencher) { remove_insert(b, MEDIUM); }
#[bench] fn bench_remove_insert_big   (b: &mut test::Bencher) { remove_insert(b, BIG);    }

#[bench] fn bench_remove_rev_insert_small (b: &mut test::Bencher) { remove_rev_insert(b, SMALL);  }
#[bench] fn bench_remove_rev_insert_medium(b: &mut test::Bencher) { remove_rev_insert(b, MEDIUM); }
#[bench] fn bench_remove_rev_insert_big   (b: &mut test::Bencher) { remove_rev_insert(b, BIG);    }

#[bench] fn bench_get_middle_small (b: &mut test::Bencher) { get_middle(b, SMALL);  }
#[bench] fn bench_get_middle_medium(b: &mut test::Bencher) { get_middle(b, MEDIUM); }
#[bench] fn bench_get_middle_big   (b: &mut test::Bencher) { get_middle(b, BIG);    }

#[bench] fn bench_get_none_small (b: &mut test::Bencher) { get_none(b, SMALL);  }
#[bench] fn bench_get_none_medium(b: &mut test::Bencher) { get_none(b, MEDIUM); }
#[bench] fn bench_get_none_big   (b: &mut test::Bencher) { get_none(b, BIG);    }
