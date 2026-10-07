//! Showing how to use shrinkwrap for different types of structs.

#![cfg(feature = "std")]
#![allow(dead_code)]

#[macro_use] extern crate shrinkwraprs;

#[derive(Shrinkwrap)]
struct Foo(i32);

#[derive(Shrinkwrap)]
struct Bar(i32, #[shrinkwrap(main_field)] String);

#[derive(Shrinkwrap)]
#[shrinkwrap(mutable)]
struct Baz {
  field1: String
}

#[derive(Shrinkwrap)]
struct Quux {
  field1: u32,
  #[shrinkwrap(main_field)] field2: String
}

#[derive(Shrinkwrap)]
#[shrinkwrap(mutable, unsafe_ignore_visibility)]
pub struct MyStruct {
  field1: u32
}

fn is_commercial(b: &Baz) -> bool {
  (**b).contains(".co")
}

fn main() {
  let mut email = Baz { field1: "chiya+snacks@natsumeya.jp".into() };

  println!("is_commercial: {}", is_commercial(&email));
  (*email).push_str(".com");
  println!("is_commercial: {}", is_commercial(&email));
}
