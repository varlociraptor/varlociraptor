#[doc(hidden)]
#[macro_export]
macro_rules! _parse_assignment_op {
    (+=, $($t:tt)+) => ($crate::_impl_assignment_op_internal!(AddAssign, add_assign, $($t)+););
    (-=, $($t:tt)+) => ($crate::_impl_assignment_op_internal!(SubAssign, sub_assign, $($t)+););
    (*=, $($t:tt)+) => ($crate::_impl_assignment_op_internal!(MulAssign, mul_assign, $($t)+););
    (/=, $($t:tt)+) => ($crate::_impl_assignment_op_internal!(DivAssign, div_assign, $($t)+););
    (%=, $($t:tt)+) => ($crate::_impl_assignment_op_internal!(RemAssign, rem_assign, $($t)+););
    (&=, $($t:tt)+) => ($crate::_impl_assignment_op_internal!(BitAndAssign, bitand_assign, $($t)+););
    (|=, $($t:tt)+) => ($crate::_impl_assignment_op_internal!(BitOrAssign, bitor_assign, $($t)+););
    (^=, $($t:tt)+) => ($crate::_impl_assignment_op_internal!(BitXorAssign, bitxor_assign, $($t)+););
    (<<=, $($t:tt)+) => ($crate::_impl_assignment_op_internal!(ShlAssign, shl_assign, $($t)+););
    (>>=, $($t:tt)+) => ($crate::_impl_assignment_op_internal!(ShrAssign, shr_assign, $($t)+););
}

#[doc(hidden)]
#[macro_export]
macro_rules! _impl_assignment_op_internal {
    ($ops_trait:ident, $ops_fn:ident, $lhs:ty, &$rhs:ty, $lhs_i:ident, $rhs_i:ident, $body:block) => {
        impl ::std::ops::$ops_trait<&$rhs> for $lhs {
            fn $ops_fn(&mut self, $rhs_i: &$rhs) {
                let mut $lhs_i = self;
                $body
            }
        }
    };
    ($ops_trait:ident, $ops_fn:ident, $lhs:ty, $rhs:ty, $lhs_i:ident, $rhs_i:ident, $body:block) => {
        impl ::std::ops::$ops_trait<$rhs> for $lhs {
            fn $ops_fn(&mut self, $rhs_i: $rhs) {
                let mut $lhs_i = self;
                $body
            }
        }
    };
}
