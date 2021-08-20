use crate::variants::model::AlleleFreq;
use ordered_float::NotNan;
use std::ops::Not;

#[derive(Debug, Clone)]
pub(crate) struct Log2FoldChange {
    vaf_a: AlleleFreq,
    vaf_b: AlleleFreq,
    l2fc: f64,
}

impl Log2FoldChange {
    pub(crate) fn new(vaf_a: AlleleFreq, vaf_b: AlleleFreq) -> Self {
        Self {
            vaf_a,
            vaf_b,
            l2fc: vaf_a.log2() - vaf_b.log2(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum Comparison {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct Log2FoldChangePredicate(pub(crate) Comparison, pub(crate) NotNan<f64>);

impl Log2FoldChangePredicate {
    pub(crate) fn is_true(&self, lfc: &Log2FoldChange) -> bool {
        // FIXME: equality checks should probably allow for some error
        let v = *self.1;
        match self.0 {
            Comparison::Equal => lfc.l2fc == v,
            Comparison::Greater => lfc.l2fc > v,
            Comparison::GreaterEqual => lfc.l2fc >= v,
            Comparison::Less => lfc.l2fc < v,
            Comparison::LessEqual => lfc.l2fc <= v,
            Comparison::NotEqual => lfc.l2fc != v,
        }
    }
}

impl Not for Log2FoldChangePredicate {
    type Output = Self;

    fn not(self) -> Self::Output {
        Self(!self.0, self.1)
    }
}
