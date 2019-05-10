use pest::iterators::{Pair, Pairs};

use crate::grammar::vafrange::VAFRange;

pub enum Formula {
    Conjunction {
        operands: Vec<Box<Formula>>
    },
    Disjunction {
        operands: Vec<Box<Formula>>
    },
    Negation {
        operand: Box<Formula>
    },
    Atom {
        sample: T,
        vafs: VAFRange
    }
}




fn parse_formula(pair: Pair) -> Result<>
