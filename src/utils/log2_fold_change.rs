use std::ops::Not;

use num_traits::ops;
use ordered_float::NotNan;

use crate::{
    grammar::VAFRange, utils::comparison::ComparisonOperator, variants::model::AlleleFreq,
};
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

    /// For the given VAF, infer the possible VAF range based on the log2 fold change predicate,
    /// given the vaf of the left operand.
    pub(crate) fn infer_vaf_bounds(&self, vaf: AlleleFreq) -> VAFRange {
        let projection_right_operand = vaf / self.value.exp2();
        if *projection_right_operand < 0.0 || *projection_right_operand > 1.0 {
            return VAFRange::empty();
        }
        match self.comparison {
            ComparisonOperator::Equal => VAFRange::builder()
                .inner(projection_right_operand..projection_right_operand)
                .left_exclusive(false)
                .right_exclusive(false)
                .build(),
            ComparisonOperator::Greater => VAFRange::builder()
                .inner(AlleleFreq(0.0)..projection_right_operand)
                .left_exclusive(false)
                .right_exclusive(true)
                .build(),
            ComparisonOperator::GreaterEqual => VAFRange::builder()
                .inner(AlleleFreq(0.0)..projection_right_operand)
                .left_exclusive(false)
                .right_exclusive(false)
                .build(),
            ComparisonOperator::Less => VAFRange::builder()
                .inner(projection_right_operand..AlleleFreq(1.0))
                .left_exclusive(true)
                .right_exclusive(false)
                .build(),
            ComparisonOperator::LessEqual => VAFRange::builder()
                .inner(projection_right_operand..AlleleFreq(1.0))
                .left_exclusive(false)
                .right_exclusive(false)
                .build(),
            ComparisonOperator::NotEqual => VAFRange::builder()
                .inner(AlleleFreq(0.0)..AlleleFreq(1.0))
                .left_exclusive(false)
                .right_exclusive(false)
                .build(), // not much we can do here
        }
    }

    pub(crate) fn invert(&self) -> Self {
        match self.comparison {
            ComparisonOperator::Equal => Self {
                comparison: ComparisonOperator::Equal,
                value: self.value,
            },
            ComparisonOperator::Greater => Self {
                comparison: ComparisonOperator::LessEqual,
                value: -self.value,
            },
            ComparisonOperator::GreaterEqual => Self {
                comparison: ComparisonOperator::Less,
                value: -self.value,
            },
            ComparisonOperator::Less => Self {
                comparison: ComparisonOperator::GreaterEqual,
                value: -self.value,
            },
            ComparisonOperator::LessEqual => Self {
                comparison: ComparisonOperator::Greater,
                value: -self.value,
            },
            ComparisonOperator::NotEqual => Self {
                comparison: ComparisonOperator::NotEqual,
                value: self.value,
            },
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
