use std::ops::Not;

use ordered_float::NotNan;

use crate::utils::comparison::Comparison;
use crate::variants::model::AlleleFreq;

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
    pub(crate) comparison: Comparison,
    pub(crate) value: NotNan<f64>,
}

impl Log2FoldChangePredicate {
    pub(crate) fn is_true(&self, lfc: &Log2FoldChange) -> bool {
        // FIXME: equality checks should probably allow for some error
        let v = *self.value;
        match self.comparison {
            Comparison::Equal => lfc.value == v,
            Comparison::Greater => lfc.value > v,
            Comparison::GreaterEqual => lfc.value >= v,
            Comparison::Less => lfc.value < v,
            Comparison::LessEqual => lfc.value <= v,
            Comparison::NotEqual => lfc.value != v,
        }
    }
}

impl Not for Log2FoldChangePredicate {
    type Output = Self;

    fn not(self) -> Self::Output {
        Self(!self.comparison, self.value)
    }
}
