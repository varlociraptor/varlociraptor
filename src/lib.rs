// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

#![warn(clippy::all)]

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
extern crate strum_macros;
#[macro_use]
extern crate derive_builder;
#[macro_use]
extern crate shrinkwraprs;
#[macro_use]
extern crate derefable;
#[macro_use]
extern crate typed_builder;
extern crate paste;

pub mod calling;
pub mod candidates;
pub mod cli;
pub(crate) mod conversion;
pub(crate) mod errors;
pub(crate) mod estimation;
pub mod filtration;
pub(crate) mod grammar;
pub(crate) mod reference;
pub mod testcase;
pub mod utils;
pub mod variants;

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
            name = self.name().replace('_', "-"),
            desc = desc,
            tag_name = &self.tag_name(prefix)
        )
    }
}

/// A simple event that just has a name.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct SimpleEvent {
    /// event name
    pub name: String,
}

impl SimpleEvent {
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
        }
    }
}

impl Event for SimpleEvent {
    fn name(&self) -> &str {
        &self.name
    }
}
