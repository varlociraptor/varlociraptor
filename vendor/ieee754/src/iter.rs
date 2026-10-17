use core::{mem, usize};
use {Bits, Ieee754};

/// An iterator over floating point numbers, created by `Ieee754::upto`.
pub struct Iter<T: Ieee754> {
    from: T,
    to: T,
    done: bool
}
pub fn new_iter<T: Ieee754>(from: T, to: T) -> Iter<T> {
    Iter { from, to, done: false }
}

impl<T: Ieee754> Iterator for Iter<T> {
    type Item = T;
    fn next(&mut self) -> Option<T> {
        if self.done { return None }

        let x = self.from;
        let y = x.next();
        // we've canonicalised negative zero to positive zero, and
        // we're guaranteed that neither is NaN, so comparing bitwise
        // is valid (and 20% faster for the `all` example).
        if x.bits() == self.to.bits() {
            self.done = true;
        }
        self.from = y;
        return Some(x)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        if self.done {
            return (0, Some(0))
        }

        let high_pos = 8 * mem::size_of::<T>() - 1;
        let high_mask = 1 << high_pos;

        let from_ = self.from.bits().as_u64();
        let (from, from_sign) = (from_ & !high_mask,
                                 from_ & high_mask != 0);
        let to_ = self.to.bits().as_u64();
        let (to, to_sign) = (to_ & !high_mask,
                             to_ & high_mask != 0);
        let from = if from_sign { -(from as i64) } else { from as i64 };
        let to = if to_sign { -(to as i64) } else { to as i64 };

        let distance = (to - from + 1) as u64;
        if distance <= usize::MAX as u64 {
            let d = distance as usize;
            (d, Some(d))
        } else {
            (usize::MAX, None)
        }
    }
}
impl<T: Ieee754> DoubleEndedIterator for Iter<T> {
    fn next_back(&mut self) -> Option<T> {
        if self.done { return None }

        let x = self.to;
        let y = x.prev();
        if x == self.from {
            self.done = true;
        }
        self.to = y;
        return Some(x)
    }
}
