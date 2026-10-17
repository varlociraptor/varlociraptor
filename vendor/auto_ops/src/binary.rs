#[doc(hidden)]
#[macro_export]
macro_rules! _parse_binary_op {
    (+, $($t:tt)+) => ($crate::_impl_binary_op_internal!(Add, add, $($t)+););
    (-, $($t:tt)+) => ($crate::_impl_binary_op_internal!(Sub, sub, $($t)+););
    (*, $($t:tt)+) => ($crate::_impl_binary_op_internal!(Mul, mul, $($t)+););
    (/, $($t:tt)+) => ($crate::_impl_binary_op_internal!(Div, div, $($t)+););
    (%, $($t:tt)+) => ($crate::_impl_binary_op_internal!(Rem, rem, $($t)+););
    (&, $($t:tt)+) => ($crate::_impl_binary_op_internal!(BitAnd, bitand, $($t)+););
    (|, $($t:tt)+) => ($crate::_impl_binary_op_internal!(BitOr, bitor, $($t)+););
    (^, $($t:tt)+) => ($crate::_impl_binary_op_internal!(BitXor, bitxor, $($t)+););
    (<<, $($t:tt)+) => ($crate::_impl_binary_op_internal!(Shl, shl, $($t)+););
    (>>, $($t:tt)+) => ($crate::_impl_binary_op_internal!(Shr, shr, $($t)+););
}

#[doc(hidden)]
#[macro_export]
macro_rules! _impl_binary_op_internal {
    ($ops_trait:ident, $ops_fn:ident, &$lhs:ty, &$rhs:ty, $out:ty, $lhs_i:ident, $rhs_i:ident, $body:block) => {
        $crate::_impl_binary_op_borrowed_borrowed!(
            $ops_trait, $ops_fn, $lhs, $rhs, $out, $lhs_i, $rhs_i, $body
        );
    };
    ($ops_trait:ident, $ops_fn:ident, &$lhs:ty, $rhs:ty, $out:ty, $lhs_i:ident, $rhs_i:ident, $body:block) => {
        $crate::_impl_binary_op_borrowed_owned!(
            $ops_trait, $ops_fn, $lhs, $rhs, $out, $lhs_i, $rhs_i, $body
        );
    };
    ($ops_trait:ident, $ops_fn:ident, $lhs:ty, &$rhs:ty, $out:ty, $lhs_i:ident, $rhs_i:ident, $body:block) => {
        $crate::_impl_binary_op_owned_borrowed!(
            $ops_trait, $ops_fn, $lhs, $rhs, $out, $lhs_i, $rhs_i, $body
        );
    };
    ($ops_trait:ident, $ops_fn:ident, $lhs:ty, $rhs:ty, $out:ty, $lhs_i:ident, $rhs_i:ident, $body:block) => {
        $crate::_impl_binary_op_owned_owned!(
            $ops_trait, $ops_fn, $lhs, $rhs, $out, $lhs_i, $rhs_i, $body
        );
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! _impl_binary_op_owned_owned {
    ($ops_trait:ident, $ops_fn:ident, $lhs:ty, $rhs:ty, $out:ty, $lhs_i:ident, $rhs_i:ident, $body:block) => {
        impl ::std::ops::$ops_trait<$rhs> for $lhs {
            type Output = $out;

            fn $ops_fn(self, $rhs_i: $rhs) -> Self::Output {
                let $lhs_i = self;
                $body
            }
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! _impl_binary_op_owned_borrowed {
    ($ops_trait:ident, $ops_fn:ident, $lhs:ty, $rhs:ty, $out:ty, $lhs_i:ident, $rhs_i:ident, $body:block) => {
        impl ::std::ops::$ops_trait<&$rhs> for $lhs {
            type Output = $out;

            fn $ops_fn(self, $rhs_i: &$rhs) -> Self::Output {
                let $lhs_i = self;
                $body
            }
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! _impl_binary_op_borrowed_owned {
    ($ops_trait:ident, $ops_fn:ident, $lhs:ty, $rhs:ty, $out:ty, $lhs_i:ident, $rhs_i:ident, $body:block) => {
        impl ::std::ops::$ops_trait<$rhs> for &$lhs {
            type Output = $out;

            fn $ops_fn(self, $rhs_i: $rhs) -> Self::Output {
                let $lhs_i = self;
                $body
            }
        }
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! _impl_binary_op_borrowed_borrowed {
    ($ops_trait:ident, $ops_fn:ident, $lhs:ty, $rhs:ty, $out:ty, $lhs_i:ident, $rhs_i:ident, $body:block) => {
        impl ::std::ops::$ops_trait<&$rhs> for &$lhs {
            type Output = $out;

            fn $ops_fn(self, $rhs_i: &$rhs) -> Self::Output {
                let $lhs_i = self;
                $body
            }
        }
    };
}
