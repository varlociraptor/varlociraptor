use core::ops::Not;
use std::fmt::{Display, Formatter};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ComparisonOperator {
    Equal,
    Greater,
    GreaterEqual,
    Less,
    LessEqual,
    NotEqual,
}

impl Display for ComparisonOperator {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            ComparisonOperator::Equal => "==",
            ComparisonOperator::Greater => ">",
            ComparisonOperator::GreaterEqual => ">=",
            ComparisonOperator::Less => "<",
            ComparisonOperator::LessEqual => "<=",
            ComparisonOperator::NotEqual => "!=",
        };
        f.write_str(s)
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
