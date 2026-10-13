//! Unit systems and units shipped with the library (API-05, API-06).

pub mod natural;
pub mod si;

/// Implements unit arithmetic (`a * b`, `a / b`, `a.pow::<K>()`) and `f64 * unit`
/// for a unit type. Used by [`unit!`](crate::unit).
#[doc(hidden)]
#[macro_export]
macro_rules! unit_ops {
    ($t:ident) => {
        impl<Y: $crate::engine::dimension::Unit> core::ops::Mul<Y> for $t {
            type Output = $crate::engine::dimension::UnitProd<$t, Y>;
            fn mul(self, _: Y) -> Self::Output {
                Default::default()
            }
        }
        impl<Y: $crate::engine::dimension::Unit> core::ops::Div<Y> for $t {
            type Output = $crate::engine::dimension::UnitQuot<$t, Y>;
            fn div(self, _: Y) -> Self::Output {
                Default::default()
            }
        }
        impl $t {
            /// This unit raised to the integer power `E` (`-20..=20`).
            pub fn pow<const E: i32>(
                self,
            ) -> $crate::engine::dimension::UnitPow<$t, <$crate::utility::typelevel::ConstInt<{ E }> as $crate::utility::typelevel::ToInt>::Out>
            where
                $crate::utility::typelevel::ConstInt<{ E }>: $crate::utility::typelevel::ToInt,
            {
                Default::default()
            }
        }
        impl core::ops::Mul<$t> for f64 {
            type Output = $crate::Quantity<<$t as $crate::engine::dimension::Unit>::Dim>;
            fn mul(self, _: $t) -> Self::Output {
                $crate::Quantity::of(self, $t)
            }
        }
    };
}

/// Declares a named unit with its dimension and its factor to the base unit
/// (EXT-02), with unit arithmetic and `f64 * unit`.
///
/// ```
/// use tensorcheck::units::si::Length;
/// tensorcheck::unit!(/** Light year. */ pub light_year: Length = 9.460_730_472_580_8e15);
/// let d = 2.0 * light_year;
/// ```
#[macro_export]
macro_rules! unit {
    ($(#[$doc:meta])* $vis:vis $t:ident : $dim:ty = $factor:expr) => {
        $(#[$doc])*
        #[allow(non_camel_case_types)]
        #[derive(Clone, Copy, Debug, Default)]
        $vis struct $t;
        impl $crate::engine::dimension::Unit for $t {
            type Dim = $dim;
            const FACTOR: f64 = $factor;
        }
        $crate::unit_ops!($t);
    };
}

// Composite units built with `*`, `/` and `pow` support the same operations.
mod composite {
    use crate::Quantity;
    use crate::engine::dimension::{Unit, UnitPow, UnitProd, UnitQuot};
    use crate::utility::typelevel::{ConstInt, ToInt};
    use core::ops::{Div, Mul};

    macro_rules! composite_ops {
        ($t:ident) => {
            impl<X, Z, Y: Unit> Mul<Y> for $t<X, Z> {
                type Output = UnitProd<$t<X, Z>, Y>;
                fn mul(self, _: Y) -> Self::Output {
                    Default::default()
                }
            }
            impl<X, Z, Y: Unit> Div<Y> for $t<X, Z> {
                type Output = UnitQuot<$t<X, Z>, Y>;
                fn div(self, _: Y) -> Self::Output {
                    Default::default()
                }
            }
            impl<X, Z> $t<X, Z> {
                /// This unit raised to the integer power `E` (`-20..=20`).
                pub fn pow<const E: i32>(self) -> UnitPow<$t<X, Z>, <ConstInt<{ E }> as ToInt>::Out>
                where
                    ConstInt<{ E }>: ToInt,
                {
                    Default::default()
                }
            }
            impl<X, Z> Mul<$t<X, Z>> for f64
            where
                $t<X, Z>: Unit,
            {
                type Output = Quantity<<$t<X, Z> as Unit>::Dim>;
                fn mul(self, _: $t<X, Z>) -> Self::Output {
                    Quantity::from_store(crate::engine::dimension::to_base::<$t<X, Z>>(self))
                }
            }
        };
    }
    composite_ops!(UnitProd);
    composite_ops!(UnitQuot);
    composite_ops!(UnitPow);
}
