#![allow(unused_variables, dead_code)]

#[macro_use] extern crate shrinkwraprs;
extern crate core;

#[derive(Shrinkwrap)]
struct Email(String);

#[derive(Shrinkwrap)]
struct CodeSpan(u64, u64, #[shrinkwrap(main_field)] String);

#[derive(Shrinkwrap)]
struct PhoneNumber {
  number: String
}

#[derive(Shrinkwrap)]
struct FileContents {
  #[shrinkwrap(main_field)] contents: String,
  linked_inodes: u64
}

#[test]
fn test_tuple_can_deref() {
  let email = Email("chiya+snacks@natsumeya.jp".into());

  assert!(email.contains("+"));
}

#[test]
fn test_nary_tuple_can_deref() {
  let span = CodeSpan(0, 24, "  impl  ".into());

  assert_eq!(span.trim(), "impl");
}

#[test]
fn test_single_can_deref() {
  let number = PhoneNumber {
    number: "+1 (800) 273-8255".into()
  };
  let is_collect_call = number.contains("(800)");

  assert!(is_collect_call);
}

#[test]
fn test_multi_can_deref() {
  let contents = FileContents {
    contents: "fjkfdlsjfkdlsjflks".into(),
    linked_inodes: 3
  };

  assert!(contents.len() > 0);
}
