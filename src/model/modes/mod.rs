// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

pub mod generic;
pub mod tumor;

use crate::grammar;

pub trait UniverseDrivenPrior {
    fn set_universe(&mut self, universe: grammar::SampleInfo<grammar::VAFUniverse>);
}