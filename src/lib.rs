#![cfg_attr(feature = "flame_it", feature(plugin))]
#![cfg_attr(feature = "flame_it", plugin(flamer))]
// activate flame for the whole crate

extern crate bio;
extern crate rust_htslib;
#[macro_use]
extern crate log;
extern crate itertools;
extern crate itertools_num;
extern crate rgsl;
#[macro_use]
extern crate approx;
extern crate ndarray;
extern crate ordered_float;
extern crate rusty_machine;
extern crate serde;
#[macro_use]
extern crate serde_derive;
#[macro_use]
extern crate quick_error;
extern crate csv;
#[macro_use]
extern crate lazy_static;
extern crate regex;
extern crate statrs;
extern crate vec_map;

#[cfg(feature = "flame_it")]
extern crate flame;

pub mod call;
pub mod constants;
pub mod estimation;
pub mod model;
pub mod utils;

pub use estimation::alignment_properties::{AlignmentProperties, InsertSize};
pub use model::likelihood;
pub use model::priors;
pub use model::sample::Sample;

quick_error! {
    #[derive(Debug)]
    pub enum BCFError {
        MissingTag(name: String) {
            description("unexpected missing tag")
            display("expected tag {} missing from BCF record", name)
        }
        InvalidRecord(msg: String) {
            description("invalid record")
            display("{}", msg)
        }
    }
}

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
            name = self.name(),
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
pub struct SimpleEvent {
    /// event name
    pub name: String,
}

impl Event for SimpleEvent {
    fn name(&self) -> &str {
        &self.name
    }
}
