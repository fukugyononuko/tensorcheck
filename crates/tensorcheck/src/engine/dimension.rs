//! Dimensions and units (CMP-E01, DD-03).
//!
//! A dimension is a unit system plus a list of integer exponents, one per base
//! dimension. Products and quotients add and subtract the exponents at compile time.

use core::marker::PhantomData;

use crate::utility::diagnostics::{UnitDimensionIs, UnitMatches};
use crate::utility::typelevel::{Bool, IntList, IntListEq, Len, MapMul, ZipAdd, ZipSub};

/// A system of base dimensions (EXT-01), for example SI with seven base dimensions.
pub trait UnitSystem: 'static {
    /// The exponent list of a dimensionless quantity in this system.
    type Dimensionless;
    /// Symbols of the base units, used when printing dimensions.
    const SYMBOLS: &'static [&'static str];
}

/// A physical dimension: unit system `Sys` and exponent list `Exps`.
pub struct Dim<Sys, Exps>(PhantomData<(Sys, Exps)>);

/// Product of two dimensions of the same unit system.
pub trait DimMul<B> {
    /// The product.
    type Out;
}
impl<Sys, A: ZipAdd<B>, B> DimMul<Dim<Sys, B>> for Dim<Sys, A> {
    type Out = Dim<Sys, <A as ZipAdd<B>>::Out>;
}

/// Quotient of two dimensions of the same unit system.
pub trait DimDiv<B> {
    /// The quotient.
    type Out;
}
impl<Sys, A: ZipSub<B>, B> DimDiv<Dim<Sys, B>> for Dim<Sys, A> {
    type Out = Dim<Sys, <A as ZipSub<B>>::Out>;
}

/// Integer power of a dimension.
pub trait DimPow<K> {
    /// `Self^K`.
    type Out;
}
impl<Sys, A: MapMul<K>, K> DimPow<K> for Dim<Sys, A> {
    type Out = Dim<Sys, <A as MapMul<K>>::Out>;
}

/// `A · B`.
pub type Prod<A, B> = <A as DimMul<B>>::Out;
/// `A / B`.
pub type Quot<A, B> = <A as DimDiv<B>>::Out;
/// `A^K` for a type-level integer `K`.
pub type Pow<A, K> = <A as DimPow<K>>::Out;

/// Reads a dimension's exponents and symbols at run time (for printing).
pub trait DimInfo {
    /// Number of base dimensions.
    const COUNT: usize;
    /// Base unit symbols.
    const SYMBOLS: &'static [&'static str];
    /// Writes the exponents into `out`.
    fn exponents(out: &mut [i32]);
}
impl<Sys: UnitSystem, Exps: IntList + Len> DimInfo for Dim<Sys, Exps> {
    const COUNT: usize = Exps::LEN;
    const SYMBOLS: &'static [&'static str] = Sys::SYMBOLS;
    fn exponents(out: &mut [i32]) {
        Exps::write(out)
    }
}

/// A unit (EXT-02): a dimension and the factor to the base unit of its system.
pub trait Unit: Copy + 'static {
    /// The dimension.
    type Dim: DimInfo;
    /// One of this unit equals `FACTOR` base units.
    const FACTOR: f64;
}

/// Equality of two dimensions of the same unit system, as a type-level boolean.
pub trait DimEq<B> {
    /// [`True`](crate::utility::typelevel::True) when equal.
    type Out: Bool;
}
impl<Sys, A: IntListEq<B>, B> DimEq<Dim<Sys, B>> for Dim<Sys, A> {
    type Out = <A as IntListEq<B>>::Out;
}

impl<U: Unit, D> UnitMatches<D> for U
where
    U::Dim: DimEq<D>,
    <U::Dim as DimEq<D>>::Out: UnitDimensionIs<U>,
{
}

/// The product of two units, written `a * b`.
pub struct UnitProd<A, B>(PhantomData<(A, B)>);
/// The quotient of two units, written `a / b`.
pub struct UnitQuot<A, B>(PhantomData<(A, B)>);
/// The integer power of a unit, written `u.pow::<K>()`.
pub struct UnitPow<A, K>(PhantomData<(A, K)>);

macro_rules! marker_impls {
    ($t:ident) => {
        impl<X, Y> Clone for $t<X, Y> {
            fn clone(&self) -> Self {
                *self
            }
        }
        impl<X, Y> Copy for $t<X, Y> {}
        impl<X, Y> Default for $t<X, Y> {
            fn default() -> Self {
                $t(PhantomData)
            }
        }
        impl<X, Y> core::fmt::Debug for $t<X, Y> {
            fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
                f.write_str(stringify!($t))
            }
        }
    };
}
marker_impls!(UnitProd);
marker_impls!(UnitQuot);
marker_impls!(UnitPow);

impl<A: Unit, B: Unit> Unit for UnitProd<A, B>
where
    A::Dim: DimMul<B::Dim>,
    Prod<A::Dim, B::Dim>: DimInfo + 'static,
{
    type Dim = Prod<A::Dim, B::Dim>;
    const FACTOR: f64 = A::FACTOR * B::FACTOR;
}
impl<A: Unit, B: Unit> Unit for UnitQuot<A, B>
where
    A::Dim: DimDiv<B::Dim>,
    Quot<A::Dim, B::Dim>: DimInfo + 'static,
{
    type Dim = Quot<A::Dim, B::Dim>;
    const FACTOR: f64 = A::FACTOR / B::FACTOR;
}
impl<A: Unit, K: crate::utility::typelevel::Int> Unit for UnitPow<A, K>
where
    A::Dim: DimPow<K>,
    Pow<A::Dim, K>: DimInfo + 'static,
{
    type Dim = Pow<A::Dim, K>;
    const FACTOR: f64 = powi_const(A::FACTOR, K::VALUE);
}

/// `x^k` evaluated at compile time.
pub const fn powi_const(x: f64, k: i32) -> f64 {
    let mut r = 1.0;
    let mut i = 0;
    let n = if k < 0 { -k } else { k };
    while i < n {
        r *= x;
        i += 1;
    }
    if k < 0 { 1.0 / r } else { r }
}

/// Converts a value given in unit `U` to the base unit.
pub fn to_base<U: Unit>(x: f64) -> f64 {
    x * U::FACTOR
}

/// Converts a value in the base unit to unit `U`.
pub fn from_base<U: Unit>(x: f64) -> f64 {
    x / U::FACTOR
}
