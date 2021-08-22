use crate::variants::model::AlleleFreq;
use core::ops::Not;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ComparisonOperator {
    Equal,
    Greater,
    GreaterEqual,
    Less,
    LessEqual,
    NotEqual,
}

impl ComparisonOperator {
    pub(crate) fn is_true(&self, vaf_a: AlleleFreq, vaf_b: AlleleFreq) -> bool {
        match self {
            ComparisonOperator::Equal => relative_eq!(*vaf_a, *vaf_b),
            ComparisonOperator::Greater => vaf_a > vaf_b,
            ComparisonOperator::GreaterEqual => vaf_a >= vaf_b,
            ComparisonOperator::Less => vaf_a < vaf_b,
            ComparisonOperator::LessEqual => vaf_a <= vaf_b,
            ComparisonOperator::NotEqual => relative_ne!(*vaf_a, *vaf_b),
        }
    }
}

impl Not for ComparisonOperator {
    type Output = Self;

    fn not(self) -> Self::Output {
        match self {
            ComparisonOperator::Equal => ComparisonOperator::NotEqual,
            ComparisonOperator::Greater => ComparisonOperator::LessEqual,
            ComparisonOperator::GreaterEqual => ComparisonOperator::Less,
            ComparisonOperator::Less => ComparisonOperator::GreaterEqual,
            ComparisonOperator::LessEqual => ComparisonOperator::Greater,
            ComparisonOperator::NotEqual => ComparisonOperator::Equal,
        }
    }
}
