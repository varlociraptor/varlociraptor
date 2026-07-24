#![allow(unused_variables, dead_code)]

#[macro_use] extern crate shrinkwraprs;
extern crate core;

#[derive(Shrinkwrap)]
#[shrinkwrap(mutable)]
struct Email(String);

#[derive(Shrinkwrap)]
#[shrinkwrap(mutable)]
struct CodeSpan(u64, u64, #[shrinkwrap(main_field)] String);

#[derive(Shrinkwrap)]
#[shrinkwrap(mutable)]
struct PhoneNumber {
  number: String
}

#[derive(Shrinkwrap)]
#[shrinkwrap(mutable)]
struct FileContents {
  #[shrinkwrap(main_field)] contents: String,
  linked_inodes: u64
}

#[test]
fn test_tuple_can_deref_mut() {
  let mut email = Email("chiya+snacks@natsumeya.jp".into());

  email.push_str(".co");
}

#[test]
fn test_nary_tuple_can_deref_mut() {
  let mut span = CodeSpan(0, 24, "  impl  ".into());

  span.push_str("!Sync for MyCell");
}

#[test]
fn test_single_can_deref_mut() {
  let mut number = PhoneNumber {
    number: "+1 (800) 273-8255".into()
  };

  number.push_str(" (20)");
}

#[test]
fn test_multi_can_deref_mut() {
  let mut contents = FileContents {
    contents: "fjkfdlsjfkdlsjflks".into(),
    linked_inodes: 3
  };

  contents.push_str("fdjskl");
}
