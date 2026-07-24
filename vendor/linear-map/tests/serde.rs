#![cfg(feature = "serde_impl")]

extern crate linear_map;
use linear_map::LinearMap;

extern crate serde_test;
use serde_test::{Token, assert_tokens};

#[test]
fn test_ser_de_empty() {
    let map = LinearMap::<char, u32>::new();

    assert_tokens(&map, &[
        Token::Map { len: Some(0) },
        Token::MapEnd,
    ]);
}

#[test]
fn test_ser_de() {
    let mut map = LinearMap::new();
    map.insert('b', 20);
    map.insert('a', 10);
    map.insert('c', 30);

    assert_tokens(&map, &[
        Token::Map {len: Some(3) },
            Token::Char('b'),
            Token::I32(20),

            Token::Char('a'),
            Token::I32(10),

            Token::Char('c'),
            Token::I32(30),
        Token::MapEnd,
    ]);
}

mod set {
    use serde_test::{Token, assert_tokens};
    use linear_map::set::LinearSet;

    #[test]
    fn test_ser_de_empty() {
        let set = LinearSet::<char>::new();
        assert_tokens(&set, &[
            Token::Seq { len: Some(0) },
            Token::SeqEnd,
        ]);
    }

    #[test]
    fn test_ser_de() {
        let mut set = LinearSet::new();
        set.insert('b');
        set.insert('a');
        set.insert('c');

        assert_tokens(&set, &[
            Token::Seq { len: Some(3) },
            Token::Char('b'),
            Token::Char('a'),
            Token::Char('c'),
            Token::SeqEnd,
        ]);
    }
}
