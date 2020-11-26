// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

pub(crate) mod generic;
pub(crate) mod tumor;

use crate::grammar;

pub(crate) trait UpdatablePrior {
    fn set_universe(&mut self, universe: grammar::SampleInfo<grammar::VAFUniverse>);
    fn set_ploidies(&mut self, ploidies: grammar::SampleInfo<Option<u32>>);
}
