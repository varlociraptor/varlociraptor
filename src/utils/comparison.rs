use core::ops::Not;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Comparison {
    Equal,
    Greater,
    GreaterEqual,
    Less,
    LessEqual,
    NotEqual,
}

impl Not for Comparison {
    type Output = Self;

    fn not(self) -> Self::Output {
        match self {
            Comparison::Equal => Comparison::NotEqual,
            Comparison::Greater => Comparison::LessEqual,
            Comparison::GreaterEqual => Comparison::Less,
            Comparison::Less => Comparison::GreaterEqual,
            Comparison::LessEqual => Comparison::Greater,
            Comparison::NotEqual => Comparison::Equal,
        }
    }
}
