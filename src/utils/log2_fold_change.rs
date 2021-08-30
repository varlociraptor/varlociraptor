use std::ops::Not;

use ordered_float::NotNan;

use crate::utils::comparison::ComparisonOperator;
use std::fmt::{Display, Formatter};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Derefable)]
pub(crate) struct Log2FoldChange(#[deref] NotNan<f64>);

impl Log2FoldChange {
    /// Calculate log2 fold change of a by b (i.e. (a / b).log2()).
    /// For a == b == 0.0, returns 0.0.
    /// Returns None if fold change cannot be calculated (if both a and b are zero).
    pub(crate) fn new(a: NotNan<f64>, b: NotNan<f64>) -> Self {
        if *a == 0.0 && *b == 0.0 {
            Log2FoldChange(NotNan::new(0.0).unwrap())
        } else {
            let value = a.log2() - b.log2();
            Log2FoldChange(NotNan::new(value).expect("bug: NaN when calculating Log2FoldChange"))
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct Log2FoldChangePredicate {
    pub(crate) comparison: ComparisonOperator,
    pub(crate) value: NotNan<f64>,
}

impl Display for Log2FoldChangePredicate {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        self.comparison.fmt(f)?;
        f.write_fmt(format_args!(" {}", self.value))
    }
}

impl Log2FoldChangePredicate {
    pub(crate) fn is_true(&self, lfc: &Log2FoldChange) -> bool {
        let v = *self.value;
        let lfc = *lfc.0;
        match self.comparison {
            ComparisonOperator::Equal => relative_eq!(lfc, v),
            ComparisonOperator::Greater => lfc > v,
            ComparisonOperator::GreaterEqual => lfc >= v,
            ComparisonOperator::Less => lfc < v,
            ComparisonOperator::LessEqual => lfc <= v,
            ComparisonOperator::NotEqual => relative_ne!(lfc, v),
        }
    }
}

impl Not for Log2FoldChangePredicate {
    type Output = Self;

    fn not(self) -> Self::Output {
        Self {
            comparison: !self.comparison,
            value: self.value,
        }
    }
}
