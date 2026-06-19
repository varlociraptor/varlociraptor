
extern crate num_traits;

mod linspace;
mod ext;

pub use linspace::linspace;
#[doc(no_inline)]
pub use structs::Linspace;
pub use ext::ItertoolsNum;

/// The concrete iterator structs.
pub mod structs {
    pub use ext::Cumsum;
    pub use linspace::Linspace;
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
    }
}
