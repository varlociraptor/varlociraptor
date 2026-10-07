//! Example usages of the mapping functions.

#![cfg(feature = "std")]
#![allow(unused_variables)]

#[macro_use] extern crate shrinkwraprs;

#[derive(Shrinkwrap)]
#[shrinkwrap(transformers)]
pub struct Email(String);
#[derive(Debug)]
pub struct Sanitized(String);

fn main() {
  let mut email = Email("aoi.miyamori@musashino".into());

  email
    .transform(|s| s.push_str(".co"))
    .transform(|s| s.push_str(".jp"));

  println!("{}", *email);

  // email has been moved out of after siphoning
  let output = email.siphon(Sanitized);

  println!("{:?}", output);
}
