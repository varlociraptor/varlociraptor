#![cfg(feature = "std")]

#[macro_use] extern crate shrinkwraprs;

#[derive(Shrinkwrap)]
pub struct Email<'a>(&'a str);

#[derive(Shrinkwrap)]
pub struct GenericType<T>(T);

fn main() {
  let email = Email("satsuki.kiryuuin@honnouji.edu");

  println!("email len: {}", email.len());

  let generic = GenericType(20_i32);

  println!("inside value: {}", *generic);
}
