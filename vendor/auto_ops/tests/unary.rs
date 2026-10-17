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

    impl_op!(!|a: kong::Donkey| -> kong::Diddy { kong::Diddy::new(a.bananas) });
    #[test]
    fn not() {
        assert_eq!(kong::Diddy::new(3), !kong::Donkey::new(3));
    }

    impl_op!(-|a: kong::Donkey| -> kong::Diddy { kong::Diddy::new(-a.bananas) });
    #[test]
    fn neg() {
        assert_eq!(kong::Diddy::new(-3), -kong::Donkey::new(3));
    }
}

mod impl_op_variants {
    use super::*;

    impl_op!(!|a: kong::Diddy| -> kong::Dixie { kong::Dixie::new(a.bananas) });
    #[test]
    fn owned() {
        assert_eq!(kong::Dixie::new(4), !kong::Diddy::new(4));
    }

    impl_op!(!|a: &kong::Diddy| -> kong::Dixie { kong::Dixie::new(a.bananas) });
    #[test]
    fn borrowed() {
        assert_eq!(kong::Dixie::new(4), !&kong::Diddy::new(4));
    }
}

mod impl_op_ex_variants {
    use super::*;

    impl_op_ex!(!|a: kong::Dixie| -> kong::Donkey { kong::Donkey::new(a.bananas) });
    #[test]
    fn owned() {
        assert_eq!(kong::Donkey::new(4), !kong::Dixie::new(4));
    }

    impl_op_ex!(-|a: &kong::Dixie| -> kong::Donkey { kong::Donkey::new(-a.bananas) });
    #[test]
    fn borrowed() {
        assert_eq!(kong::Donkey::new(-4), -&kong::Dixie::new(4));
        assert_eq!(kong::Donkey::new(-4), -kong::Dixie::new(4));
    }
}

mod generics {
    use super::*;

    impl_op!(!|a: kong::Barrel<u32>| -> kong::Barrel<i32> { kong::Barrel::new(a.bananas as i32) });
    #[test]
    fn impl_op() {
        assert_eq!(kong::Barrel::new(3), !kong::Barrel::new(3u32));
    }

    impl_op_ex!(-|a: &kong::Barrel<u32>| -> kong::Barrel<i32> {
        kong::Barrel::new(-(a.bananas as i32))
    });
    #[test]
    fn impl_op_ex() {
        assert_eq!(kong::Barrel::new(-3), -&kong::Barrel::new(3u32));
        assert_eq!(kong::Barrel::new(-3), -kong::Barrel::new(3u32));
    }
}

mod multiline {
    use super::*;

    impl_op!(!|a: kong::Barrel<i32>| -> kong::Donkey {
        let bananas = a.bananas;
        kong::Donkey::new(bananas)
    });
    #[test]
    fn impl_op() {
        assert_eq!(kong::Donkey::new(3), !kong::Barrel::new(3));
    }

    impl_op_ex!(-|a: &kong::Barrel<i32>| -> kong::Donkey {
        let bananas = a.bananas;
        kong::Donkey::new(-bananas)
    });
    #[test]
    fn impl_op_ex() {
        assert_eq!(kong::Donkey::new(-3), -&kong::Barrel::new(3));
        assert_eq!(kong::Donkey::new(-3), -kong::Barrel::new(3));
    }
}
