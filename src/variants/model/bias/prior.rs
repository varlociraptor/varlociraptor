use bio::stats::{LogProb, bayesian::model::Likelihood};

use crate::variants::{model::AlleleFreq, sample::Pileup};

use super::Bias;

pub(crate) struct BiasPrior<B> {
    bias: B
}

impl<B> Likelihood for BiasPrior<B> 
where
    B: Bias
{
    type Event = AlleleFreq;
    type Data = Pileup;

    fn compute(
        &self, 
        event: &Self::Event, 
        data: &Self::Data, 
        payload: &mut ()
    ) -> LogProb {
        for obs in data {
            self.bias.
        }
    }
}