#![allow(unused_variables)]

#[macro_use] extern crate shrinkwraprs;
extern crate core;

#[derive(Shrinkwrap)]
#[shrinkwrap(transformers)]
pub struct Email(String);

#[test]
fn test_map_mut() {
  let mut email = Email("aoi.miyamori@musashino".into());
  let orig_len = email.len();

  email
    .transform(|s| s.push_str(".co"))
    .transform(|s| s.push_str(".jp"));

  let new_len = email.len();
  let siphoned_len = email.siphon(|s| s.len());

  assert_eq!(new_len, orig_len + 6);
  assert_eq!(new_len, siphoned_len);
}
