// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use bio::stats::bayesian::model::Prior;
use bio::stats::LogProb;

use crate::model::likelihood::Event;

#[derive(Default, Clone)]
pub struct FlatPrior {}

impl FlatPrior {
    pub fn new() -> Self {
        FlatPrior {}
    }
}

impl Prior for FlatPrior {
    type Event = Vec<Event>;

    fn compute(&self, _event: &Self::Event) -> LogProb {
        LogProb::ln_one()
    }
}
