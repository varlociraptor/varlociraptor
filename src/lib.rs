// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

#[macro_use]
extern crate log;
#[macro_use]
extern crate approx;
#[macro_use]
extern crate serde_derive;
#[macro_use]
extern crate lazy_static;
extern crate askama;
#[macro_use]
extern crate derive_new;
#[macro_use]
extern crate pest_derive;
#[macro_use]
extern crate getset;
#[macro_use]
extern crate snafu;

pub mod calling;
pub mod cli;
pub mod constants;
pub mod conversion;
pub mod errors;
pub mod estimation;
pub mod filtration;
pub mod grammar;
pub mod model;
pub mod testcase;
pub mod utils;

pub use crate::estimation::alignment_properties::{AlignmentProperties, InsertSize};
pub use crate::model::likelihood;
pub use crate::model::modes;
pub use crate::model::sample::Sample;

/// Event to call.
pub trait Event {
    fn name(&self) -> &str;

    fn tag_name(&self, prefix: &str) -> String {
        format!("{}_{}", prefix, self.name().to_ascii_uppercase())
    }

    fn header_entry(&self, prefix: &str, desc: &str) -> String {
        format!(
            "##INFO=<ID={tag_name},Number=A,Type=Float,\
             Description=\"{desc} {name} variant\">",
            name = self.name().replace("_", "-"),
            desc = desc,
            tag_name = &self.tag_name(prefix)
        )
    }
}

/// Complement of other given events (i.e. 1 - Pr(other events)).
pub struct ComplementEvent {
    /// event name
    pub name: String,
}

impl Event for ComplementEvent {
    fn name(&self) -> &str {
        &self.name
    }
}

/// A simple event that just has a name.
#[derive(Debug)]
pub struct SimpleEvent {
    /// event name
    pub name: String,
}

impl Event for SimpleEvent {
    fn name(&self) -> &str {
        &self.name
    }
}
