use std::cmp;

use bio::stats::bayesian::model::{Likelihood, Posterior};
use bio::stats::LogProb;

use crate::model;
use crate::model::likelihood;
use crate::model::sample::Pileup;
use crate::model::{AlleleFreq, ContinuousAlleleFreqs};

pub struct GenericPosterior {
    resolutions: Vec<usize>,
    // for each sample, this vector points to the contaminating one (if any)
    contaminations: Vec<Option<usize>>
}

impl GenericPosterior {
    fn grid_points(&self, pileups: &Self::Data) -> Vec<usize> {
        pileups.iter().zip(self.resolutions.iter()).map(|(pileup, res)| {
            let n_obs = pileup.len();
            let mut n = cmp::min(cmp::max(n_obs + 1, 5), res);
            if n % 2 == 0 {
                n += 1;
            }
            n
        }).collect()
    }
}

impl Posterior for GenericPosterior {
    type BaseEvent = Vec<likelihood::Event>;
    type Event = Vec<model::Event<ContinuousAlleleFreqs>>;
    type Data = Vec<Pileup>;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        event: &Self::Event,
        pileups: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        let grid_points = self.grid_points(pileups);


    }
}
