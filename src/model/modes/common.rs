use bio::stats::bayesian::model::Prior;
use bio::stats::LogProb;

use crate::model::AlleleFreq;

pub struct FlatPairPrior {}

impl Prior for FlatPairPrior {
    type Event = (AlleleFreq, AlleleFreq);

    fn compute(&self, _event: &(AlleleFreq, AlleleFreq)) -> LogProb {
        LogProb::ln_one()
    }
}
