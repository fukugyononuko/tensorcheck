//! Macros that define spaces, charts and labels (DD-18).

/// Defines a space with `n` dimensions and the reference dimension of its
/// coordinates (EXT-03).
///
/// ```
/// use tensorcheck::units::si::Dimensionless;
/// tensorcheck::space!(pub FieldSpace, 2, Dimensionless);
/// ```
#[macro_export]
macro_rules! space {
    ($(#[$m:meta])* $vis:vis $name:ident, $n:literal, $refdim:ty) => {
        $(#[$m])*
        #[derive(Clone, Copy, Debug)]
        $vis struct $name;
        impl $crate::access::component::Extent for $name {
            const N: usize = $n;
            type Arr<T: Copy> = [T; $n];
            type SymArr<T: Copy> = [T; $n * ($n + 1) / 2];
            type AntiArr<T: Copy> = [T; $n * ($n - 1) / 2];
            fn from_fn<T: Copy>(f: impl FnMut(usize) -> T) -> [T; $n] {
                core::array::from_fn(f)
            }
            fn sym_from_fn<T: Copy>(f: impl FnMut(usize) -> T) -> [T; $n * ($n + 1) / 2] {
                core::array::from_fn(f)
            }
            fn anti_from_fn<T: Copy>(f: impl FnMut(usize) -> T) -> [T; $n * ($n - 1) / 2] {
                core::array::from_fn(f)
            }
            fn flatten<T: Copy>(xs: &[[T; $n]]) -> &[T] {
                xs.as_flattened()
            }
        }
        impl $crate::engine::index::Space for $name {
            type RefDim = $refdim;
            const NAME: &'static str = stringify!($name);
        }
    };
}

/// Defines a chart on a space, with one marker type per coordinate (EXT-04).
///
/// ```
/// use tensorcheck::units::si::{Dimensionless, Length};
/// tensorcheck::space!(pub Plane, 2, Length);
/// tensorcheck::chart!(pub Polar on Plane { R: Length, Theta: Dimensionless });
/// ```
///
/// Add `; non_coordinate` after the braces for a basis that is not a coordinate basis.
#[macro_export]
macro_rules! chart {
    ($(#[$m:meta])* $vis:vis $name:ident on $space:ty { $($x:ident : $d:ty),+ $(,)? } $(; $nc:ident)?) => {
        $(#[$m])*
        #[derive(Clone, Copy, Debug)]
        $vis struct $name;
        $(
            #[doc = concat!("Coordinate `", stringify!($x), "` of [`", stringify!($name), "`].")]
            #[derive(Clone, Copy, Debug)]
            $vis struct $x;
        )+
        impl $crate::engine::index::Chart for $name {
            type Space = $space;
            type Coords = $crate::chart!(@list $($x),+);
            type CoordDims = $crate::chart!(@list $($d),+);
            type Basis = $crate::chart!(@basis $($nc)?);
            const NAME: &'static str = stringify!($name);
        }
    };
    (@list $h:ty $(, $t:ty)*) => { $crate::utility::typelevel::Cons<$h, $crate::chart!(@list $($t),*)> };
    (@list) => { $crate::utility::typelevel::Nil };
    (@basis) => { $crate::engine::index::CoordinateBasis };
    (@basis non_coordinate) => { $crate::engine::index::NonCoordinateBasis };
}

/// Defines an index label with an explicit number from 40 to 63 (the library's
/// labels use 0 to 39). Two labels with the same number are treated as the same
/// label, so choose distinct numbers.
///
/// ```
/// tensorcheck::label!(pub kappa = 40);
/// ```
#[macro_export]
macro_rules! label {
    ($vis:vis $name:ident = $id:literal) => {
        #[doc = concat!("Index label `", stringify!($name), "`.")]
        #[allow(non_camel_case_types)]
        #[derive(Clone, Copy, Debug)]
        $vis struct $name;
        impl $crate::engine::index::Label for $name {
            type Id = <$crate::utility::typelevel::ConstNat<$id> as $crate::utility::typelevel::ToNat>::Out;
            const NAME: &'static str = stringify!($name);
        }
    };
}
