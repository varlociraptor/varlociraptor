use bio::stats::bayesian::model::Prior;
use bio::stats::LogProb;

use crate::model::AlleleFreq;

#[derive(Default, Clone)]
pub struct FlatPrior {}

impl FlatPrior {
    pub fn new() -> Self {
        FlatPrior {}
    }
}

impl Prior for FlatPrior {
    type Event = Vec<AlleleFreq>;

    fn compute(&self, _event: &Self::Event) -> LogProb {
        LogProb::ln_one()
    }
}
