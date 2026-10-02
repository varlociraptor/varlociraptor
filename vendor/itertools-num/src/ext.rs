
use std::ops::Add;
use num_traits::{Zero};

/// Extension trait for iterators: extra adaptors and methods
/// for numerical iterators
pub trait ItertoolsNum : Iterator {
    /// Return an iterator that produces the sequence of cumulative sums
    /// of the base iterator. The type of the sum is `S`.
    fn cumsum<S>(self) -> Cumsum<Self, S>
        where Self: Sized,
              S: Add<Self::Item, Output=S>,
              S: Zero,
    {
        cumsum(self)
    }
}

impl<I: ?Sized> ItertoolsNum for I where I: Iterator { }

pub struct Cumsum<I, S> {
    sum: S,
    iter: I,
}

fn cumsum<I, S>(iter: I) -> Cumsum<I, S>
    where I: Iterator,
          S: Add<I::Item, Output=S>,
          S: Zero,
{
    Cumsum {
        sum: Zero::zero(),
        iter: iter,
    }
}

impl<I, S> Iterator for Cumsum<I, S>
    where I: Iterator,
          S: Add<I::Item, Output=S>,
          S: Zero + Clone,
{
    type Item = S;
    fn next(&mut self) -> Option<Self::Item> {
        let z = &mut self.sum;
        self.iter.next().map(|x| {
            *z = z.clone() + x;
            z.clone()
        })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}

impl<I, S> ExactSizeIterator for Cumsum<I, S>
    where I: Iterator,
          S: Add<I::Item, Output=S>,
          S: Zero + Clone,
{ }

#[test]
fn test_cumsum() {
    let data = [1., 2., 3.];
    let mut iter = data.iter().cumsum::<f64>();
    assert_eq!(iter.next(), Some(1.));
    assert_eq!(iter.next(), Some(3.));
    assert_eq!(iter.next(), Some(6.));

}
