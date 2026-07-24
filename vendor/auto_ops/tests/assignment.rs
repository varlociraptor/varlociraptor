use auto_ops::{impl_op, impl_op_ex};

mod kong {
    #[derive(Clone, Copy, Debug, Default, PartialEq)]
    pub struct Barrel<T> {
        pub bananas: T,
    }

    impl<T> Barrel<T> {
        pub fn new(bananas: T) -> Barrel<T> {
            Barrel { bananas }
        }
    }

    #[derive(Clone, Copy, Debug, Default, PartialEq)]
    pub struct Donkey {
        pub bananas: i32,
    }

    impl Donkey {
        pub fn new(bananas: i32) -> Donkey {
            Donkey { bananas }
        }
    }

    #[derive(Clone, Copy, Debug, Default, PartialEq)]
    pub struct Diddy {
        pub bananas: i32,
    }

    impl Diddy {
        pub fn new(bananas: i32) -> Diddy {
            Diddy { bananas }
        }
    }

    #[derive(Clone, Copy, Debug, Default, PartialEq)]
    pub struct Dixie {
        pub bananas: i32,
    }

    impl Dixie {
        pub fn new(bananas: i32) -> Dixie {
            Dixie { bananas }
        }
    }
}

mod impl_op_operators {
    use super::*;

    impl_op!(+= |a: &mut kong::Donkey, b: kong::Diddy| { a.bananas += b.bananas; });
    #[test]
    fn add_assign() {
        let mut dk = kong::Donkey::new(3);
        dk += kong::Diddy::new(1);
        assert_eq!(kong::Donkey::new(3 + 1), dk);
    }

    impl_op!(-= |a: &mut kong::Donkey, b: kong::Diddy| { a.bananas -= b.bananas; });
    #[test]
    fn sub_assign() {
        let mut dk = kong::Donkey::new(3);
        dk -= kong::Diddy::new(1);
        assert_eq!(kong::Donkey::new(3 - 1), dk);
    }

    impl_op!(*= |a: &mut kong::Donkey, b: kong::Diddy| { a.bananas *= b.bananas; });
    #[test]
    fn mul_assign() {
        let mut dk = kong::Donkey::new(3);
        dk *= kong::Diddy::new(1);
        assert_eq!(kong::Donkey::new(3 * 1), dk);
    }

    impl_op!(/= |a: &mut kong::Donkey, b: kong::Diddy| { a.bananas /= b.bananas; });
    #[test]
    fn div_assign() {
        let mut dk = kong::Donkey::new(3);
        dk /= kong::Diddy::new(1);
        assert_eq!(kong::Donkey::new(3 / 1), dk);
    }

    impl_op!(%= |a: &mut kong::Donkey, b: kong::Diddy| { a.bananas %= b.bananas; });
    #[test]
    fn rem_assign() {
        let mut dk = kong::Donkey::new(3);
        dk %= kong::Diddy::new(1);
        assert_eq!(kong::Donkey::new(3 % 1), dk);
    }

    impl_op!(&= |a: &mut kong::Donkey, b: kong::Diddy| { a.bananas &= b.bananas; });
    #[test]
    fn bitand_assign() {
        let mut dk = kong::Donkey::new(3);
        dk &= kong::Diddy::new(1);
        assert_eq!(kong::Donkey::new(3 & 1), dk);
    }

    impl_op!(|= |a: &mut kong::Donkey, b: kong::Diddy| { a.bananas |= b.bananas; });
    #[test]
    fn bitor_assign() {
        let mut dk = kong::Donkey::new(3);
        dk |= kong::Diddy::new(1);
        assert_eq!(kong::Donkey::new(3 | 1), dk);
    }

    impl_op!(^= |a: &mut kong::Donkey, b: kong::Diddy| { a.bananas ^= b.bananas; });
    #[test]
    fn bitxor_assign() {
        let mut dk = kong::Donkey::new(3);
        dk ^= kong::Diddy::new(1);
        assert_eq!(kong::Donkey::new(3 ^ 1), dk);
    }

    impl_op!(<<= |a: &mut kong::Donkey, b: kong::Diddy| { a.bananas <<= b.bananas; });
    #[test]
    fn shl_assign() {
        let mut dk = kong::Donkey::new(3);
        dk <<= kong::Diddy::new(1);
        assert_eq!(kong::Donkey::new(3 << 1), dk);
    }

    impl_op!(>>= |a: &mut kong::Donkey, b: kong::Diddy| { a.bananas >>= b.bananas; });
    #[test]
    fn shr_assign() {
        let mut dk = kong::Donkey::new(3);
        dk >>= kong::Diddy::new(1);
        assert_eq!(kong::Donkey::new(3 >> 1), dk);
    }
}

mod impl_op_variants {
    use super::*;

    impl_op!(+= |a: &mut kong::Donkey, b: kong::Dixie| { a.bananas += b.bananas; });
    #[test]
    fn owned() {
        let mut dk = kong::Donkey::new(3);
        dk += kong::Dixie::new(1);
        assert_eq!(kong::Donkey::new(3 + 1), dk);
    }

    impl_op!(+= |a: &mut kong::Donkey, b: &kong::Dixie| { a.bananas += b.bananas; });
    #[test]
    fn borrowed() {
        let mut dk = kong::Donkey::new(3);
        dk += &kong::Dixie::new(1);
        assert_eq!(kong::Donkey::new(3 + 1), dk);
    }
}

mod impl_op_ex_variants {
    use super::*;

    impl_op_ex!(-= |a: &mut kong::Donkey, b: kong::Dixie| { a.bananas -= b.bananas; });
    #[test]
    fn owned() {
        let mut dk = kong::Donkey::new(3);
        dk -= kong::Dixie::new(1);
        assert_eq!(kong::Donkey::new(3 - 1), dk);
    }

    impl_op_ex!(*= |a: &mut kong::Donkey, b: &kong::Dixie| { a.bananas *= b.bananas; });
    #[test]
    fn borrowed() {
        let mut dk = kong::Donkey::new(3);
        dk *= &kong::Dixie::new(2);
        assert_eq!(kong::Donkey::new(3 * 2), dk);
        dk *= kong::Dixie::new(2);
        assert_eq!(kong::Donkey::new(3 * 2 * 2), dk);
    }
}

mod generics {
    use super::*;

    impl_op!(+= |a: &mut kong::Barrel<i32>, b: kong::Barrel<u32>| { a.bananas += b.bananas as i32; });
    #[test]
    fn impl_op() {
        let mut barrel = kong::Barrel::new(3);
        barrel += kong::Barrel::new(2u32);
        assert_eq!(kong::Barrel::new(3 + 2), barrel);
    }

    impl_op_ex!(-= |a: &mut kong::Barrel<i32>, b: &kong::Barrel<u32>| { a.bananas -= b.bananas as i32; });
    #[test]
    fn impl_op_ex() {
        let mut barrel = kong::Barrel::new(3);
        barrel -= kong::Barrel::new(2u32);
        assert_eq!(kong::Barrel::new(3 - 2), barrel);
        barrel -= &kong::Barrel::new(1u32);
        assert_eq!(kong::Barrel::new(0), barrel);
    }
}

mod multiline {
    use super::*;

    impl_op!(+= |a: &mut kong::Donkey, b: kong::Barrel<i32>| {
        a.bananas += 0;
        a.bananas += b.bananas;
    });
    #[test]
    fn impl_op() {
        let mut dk = kong::Donkey::new(3);
        dk += kong::Barrel::new(2);
        assert_eq!(kong::Donkey::new(3 + 2), dk);
    }

    impl_op_ex!(-= |a: &mut kong::Donkey, b: &kong::Barrel<i32>| {
        a.bananas += 0;
        a.bananas -= b.bananas;
    });
    #[test]
    fn impl_op_ex() {
        let mut dk = kong::Donkey::new(3);
        dk -= kong::Barrel::new(2);
        assert_eq!(kong::Donkey::new(3 - 2), dk);
        dk -= &kong::Barrel::new(1);
        assert_eq!(kong::Donkey::new(0), dk);
    }
}

mod impl_op_sequence_types {
    use super::*;

    impl_op!(+= |a: &mut kong::Donkey, b: (i32, i32)| { a.bananas += b.0; });
    #[test]
    fn tuple() {
        let mut dk = kong::Donkey::new(3);
        dk += (1, 2);
        assert_eq!(kong::Donkey::new(3 + 1), dk);
    }

    impl_op!(+= |a: &mut kong::Donkey, b: [i32; 2]| { a.bananas += b[0]; });
    #[test]
    fn array() {
        let mut dk = kong::Donkey::new(3);
        dk += [1, 2];
        assert_eq!(kong::Donkey::new(3 + 1), dk);
    }

    impl_op!(+= |a: &mut kong::Donkey, b: &[i32]| { a.bananas += b[0]; });
    #[test]
    fn slice() {
        let mut dk = kong::Donkey::new(3);
        dk += vec![1, 2].as_slice();
        assert_eq!(kong::Donkey::new(3 + 1), dk);
    }
}
