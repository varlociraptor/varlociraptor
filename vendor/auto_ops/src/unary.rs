#[doc(hidden)]
#[macro_export]
macro_rules! _parse_unary_op {
    (-, $($t:tt)+) => ($crate::_impl_unary_op_internal!(Neg, neg, $($t)+););
    (!, $($t:tt)+) => ($crate::_impl_unary_op_internal!(Not, not, $($t)+););
}

#[doc(hidden)]
#[macro_export]
macro_rules! _impl_unary_op_internal {
    ($ops_trait:ident, $ops_fn:ident, &$lhs:ty, $out:ty, $lhs_i:ident, $body:block) => {
        impl ::std::ops::$ops_trait for &$lhs {
            type Output = $out;

            fn $ops_fn(self) -> Self::Output {
                let $lhs_i = self;
                $body
            }
        }
    };
    ($ops_trait:ident, $ops_fn:ident, $lhs:ty, $out:ty, $lhs_i:ident, $body:block) => {
        impl ::std::ops::$ops_trait for $lhs {
            type Output = $out;

            fn $ops_fn(self) -> Self::Output {
                let $lhs_i = self;
                $body
            }
        }
    };
}
