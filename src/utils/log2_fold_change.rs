use std::ops::Not;

use ordered_float::NotNan;

use crate::utils::comparison::ComparisonOperator;
use crate::variants::model::AlleleFreq;
use std::fmt::{Display, Formatter};

#[derive(Debug, Clone)]
pub(crate) struct Log2FoldChange {
    vaf_a: AlleleFreq,
    vaf_b: AlleleFreq,
    value: f64,
}

impl Log2FoldChange {
    pub(crate) fn new(vaf_a: AlleleFreq, vaf_b: AlleleFreq) -> Self {
        Self {
            vaf_a,
            vaf_b,
            value: vaf_a.log2() - vaf_b.log2(),
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
        f.write_fmt(format_args!("{}", self.value))
    }
}

impl Log2FoldChangePredicate {
    pub(crate) fn is_true(&self, lfc: &Log2FoldChange) -> bool {
        let v = *self.value;
        match self.comparison {
            ComparisonOperator::Equal => relative_eq!(lfc.value, v),
            ComparisonOperator::Greater => lfc.value > v,
            ComparisonOperator::GreaterEqual => lfc.value >= v,
            ComparisonOperator::Less => lfc.value < v,
            ComparisonOperator::LessEqual => lfc.value <= v,
            ComparisonOperator::NotEqual => relative_ne!(lfc.value, v),
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
