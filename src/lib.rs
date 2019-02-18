#[macro_use]
extern crate log;
#[macro_use]
extern crate approx;
#[macro_use]
extern crate serde_derive;
#[macro_use]
extern crate quick_error;
#[macro_use]
extern crate lazy_static;

pub mod call;
pub mod constants;
pub mod estimation;
pub mod filtration;
pub mod model;
pub mod utils;

pub use crate::estimation::alignment_properties::{AlignmentProperties, InsertSize};
pub use crate::model::likelihood;
pub use crate::model::modes;
pub use crate::model::sample::Sample;

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
