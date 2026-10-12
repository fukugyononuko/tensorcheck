//! Tensors and quantities, and the operations that assemble expressions (CMP-M01, DD-05).

use core::fmt;
use core::marker::PhantomData;
use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::access::component::Slots;
use crate::engine::dimension::{DimDiv, DimInfo, DimMul, Prod, Quot, Unit, from_base, to_base};
use crate::engine::index::Tensorial;
use crate::engine::symmetry::NoSym;
use crate::utility::diagnostics::{SameDimension, UnitMatches};
use crate::utility::scalar::Scalar;
use crate::utility::typelevel::{Cons, Nil};

/// A tensor with dimension `D`, index slots `S`, symmetry `Sym`, tensoriality `K`
/// and component type `N` (API-01).
///
/// Components are stored in the base unit of the unit system. A quantity is a
/// tensor of rank 0 ([`Quantity`]).
pub struct Tensor<D, S: Slots, Sym = NoSym, K = Tensorial, N: Scalar = f64> {
    pub(crate) c: S::Store<N>,
    _m: PhantomData<(D, Sym, K)>,
}

/// A rank-0 tensor: a number with a physical dimension (API-02).
pub type Quantity<D, N = f64> = Tensor<D, Nil, NoSym, Tensorial, N>;

impl<D, S: Slots, Sym, K, N: Scalar> Clone for Tensor<D, S, Sym, K, N> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<D, S: Slots, Sym, K, N: Scalar> Copy for Tensor<D, S, Sym, K, N> {}

impl<D, S: Slots, Sym, K, N: Scalar> Tensor<D, S, Sym, K, N> {
    /// Wraps a component store (crate-internal: the store must already be in base units
    /// and satisfy the symmetry).
    #[inline(always)]
    pub(crate) fn from_store(c: S::Store<N>) -> Self {
        Tensor { c, _m: PhantomData }
    }

    /// Builds a tensor from a function of the flat (row-major) component index.
    #[inline(always)]
    pub(crate) fn from_flat_fn(f: impl FnMut(usize) -> N) -> Self {
        Self::from_store(S::from_fn(f))
    }

    /// Component number `flat` in row-major order, in base units.
    #[inline(always)]
    pub(crate) fn flat(&self, flat: usize) -> N {
        S::get(&self.c, flat)
    }

    /// Number of components.
    pub const LEN: usize = S::LEN;

    /// Applies `f` to every component.
    #[inline(always)]
    pub(crate) fn map(self, mut f: impl FnMut(N) -> N) -> Self {
        Self::from_flat_fn(|i| f(self.flat(i)))
    }

    /// The tensor with every component equal to zero.
    pub fn zero() -> Self {
        Self::from_flat_fn(|_| N::zero())
    }
}

impl<D, N: Scalar> Quantity<D, N> {
    /// The value in base units. Prefer [`Quantity::value_in`] to state the unit.
    pub fn base_value(self) -> N {
        self.c
    }
}

impl<D> Quantity<D> {
    /// Creates a quantity from a value in `unit` (API-20).
    ///
    /// ```
    /// use tensorcheck::Quantity;
    /// use tensorcheck::units::si::{km, Length};
    /// let d: Quantity<Length> = Quantity::new(1.5, km);
    /// assert_eq!(d.value_in(tensorcheck::units::si::m), 1500.0);
    /// ```
    pub fn new<U: Unit + UnitMatches<D>>(value: f64, _unit: U) -> Self {
        Self::from_store(to_base::<U>(value))
    }

    /// The value in `unit` (API-23).
    pub fn value_in<U: Unit + UnitMatches<D>>(self, _unit: U) -> f64 {
        from_base::<U>(self.c)
    }
}

// ---------------------------------------------------------------------------
// Addition and subtraction: same dimension (ERR-01). Slots, symmetry and
// tensoriality must be identical types; their dedicated messages come in later phases.

impl<D1, D2, S: Slots, Sym, K, N: Scalar> Add<Tensor<D2, S, Sym, K, N>> for Tensor<D1, S, Sym, K, N>
where
    D1: SameDimension<D2>,
{
    type Output = Tensor<D1, S, Sym, K, N>;
    #[inline(always)]
    fn add(self, rhs: Tensor<D2, S, Sym, K, N>) -> Self::Output {
        Tensor::from_flat_fn(|i| self.flat(i) + rhs.flat(i))
    }
}

impl<D1, D2, S: Slots, Sym, K, N: Scalar> Sub<Tensor<D2, S, Sym, K, N>> for Tensor<D1, S, Sym, K, N>
where
    D1: SameDimension<D2>,
{
    type Output = Tensor<D1, S, Sym, K, N>;
    #[inline(always)]
    fn sub(self, rhs: Tensor<D2, S, Sym, K, N>) -> Self::Output {
        Tensor::from_flat_fn(|i| self.flat(i) - rhs.flat(i))
    }
}

impl<D, S: Slots, Sym, K, N: Scalar> Neg for Tensor<D, S, Sym, K, N> {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        self.map(|x| -x)
    }
}

// ---------------------------------------------------------------------------
// Multiplication by a scalar quantity. At least one operand has rank 0; products of
// two tensors of positive rank are written with labels (`idx`) or `outer`.

/// Slots of the product of a rank-0 tensor and a tensor (either order).
pub trait ScalarProductSlots<S2> {
    /// The slots of the product.
    type Out: Slots;
    /// Whether the left operand is the scalar.
    const LEFT_IS_SCALAR: bool;
}
impl<S2: Slots> ScalarProductSlots<S2> for Nil {
    type Out = S2;
    const LEFT_IS_SCALAR: bool = true;
}
impl<H, T> ScalarProductSlots<Nil> for Cons<H, T>
where
    Cons<H, T>: Slots,
{
    type Out = Cons<H, T>;
    const LEFT_IS_SCALAR: bool = false;
}

impl<D1, D2, S1, S2, Sym1, Sym2, K1, K2, N> Mul<Tensor<D2, S2, Sym2, K2, N>>
    for Tensor<D1, S1, Sym1, K1, N>
where
    S1: Slots + ScalarProductSlots<S2>,
    S2: Slots,
    N: Scalar,
    D1: DimMul<D2>,
    (S1, Sym1, K1): PickTraits<S2, Sym2, K2>,
{
    type Output = Tensor<
        Prod<D1, D2>,
        <S1 as ScalarProductSlots<S2>>::Out,
        <(S1, Sym1, K1) as PickTraits<S2, Sym2, K2>>::Sym,
        <(S1, Sym1, K1) as PickTraits<S2, Sym2, K2>>::K,
        N,
    >;
    #[inline(always)]
    fn mul(self, rhs: Tensor<D2, S2, Sym2, K2, N>) -> Self::Output {
        if <S1 as ScalarProductSlots<S2>>::LEFT_IS_SCALAR {
            let k = self.flat(0);
            Tensor::from_flat_fn(|i| k * rhs.flat(i))
        } else {
            let k = rhs.flat(0);
            Tensor::from_flat_fn(|i| self.flat(i) * k)
        }
    }
}

impl<D1, D2, S1, Sym1, K1, N> Div<Quantity<D2, N>> for Tensor<D1, S1, Sym1, K1, N>
where
    S1: Slots,
    N: Scalar,
    D1: DimDiv<D2>,
{
    type Output = Tensor<Quot<D1, D2>, S1, Sym1, K1, N>;
    #[inline(always)]
    fn div(self, rhs: Quantity<D2, N>) -> Self::Output {
        let k = rhs.c;
        Tensor::from_flat_fn(|i| self.flat(i) / k)
    }
}

/// Chooses symmetry and tensoriality of a scalar product from the non-scalar operand.
pub trait PickTraits<S2, Sym2, K2> {
    /// Symmetry of the product.
    type Sym;
    /// Tensoriality of the product.
    type K;
}
impl<S2, Sym1, Sym2, K1, K2> PickTraits<S2, Sym2, K2> for (Nil, Sym1, K1) {
    type Sym = Sym2;
    type K = K2;
}
impl<H, T, Sym1, Sym2, K1, K2> PickTraits<Nil, Sym2, K2> for (Cons<H, T>, Sym1, K1) {
    type Sym = Sym1;
    type K = K1;
}

// Multiplication and division by a plain number.
impl<D, S: Slots, Sym, K> Mul<f64> for Tensor<D, S, Sym, K, f64> {
    type Output = Self;
    #[inline(always)]
    fn mul(self, k: f64) -> Self {
        self.map(|x| x * k)
    }
}
impl<D, S: Slots, Sym, K> Mul<Tensor<D, S, Sym, K, f64>> for f64 {
    type Output = Tensor<D, S, Sym, K, f64>;
    #[inline(always)]
    fn mul(self, t: Tensor<D, S, Sym, K, f64>) -> Self::Output {
        t.map(|x| self * x)
    }
}
impl<D, S: Slots, Sym, K> Div<f64> for Tensor<D, S, Sym, K, f64> {
    type Output = Self;
    #[inline(always)]
    fn div(self, k: f64) -> Self {
        self.map(|x| x / k)
    }
}

impl<D: DimInfo, S: Slots, Sym, K, N: Scalar> fmt::Debug for Tensor<D, S, Sym, K, N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if S::RANK == 0 {
            write!(f, "{:?}", self.flat(0))?;
        } else {
            write!(f, "[")?;
            for i in 0..S::LEN {
                if i > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{:?}", self.flat(i))?;
            }
            write!(f, "]")?;
        }
        write_dimension::<D>(f)
    }
}

/// Writes ` kg m s^-2` style base units for a dimension.
fn write_dimension<D: DimInfo>(f: &mut fmt::Formatter<'_>) -> fmt::Result {
    let mut exps = [0i32; 16];
    D::exponents(&mut exps[..D::COUNT]);
    for (e, sym) in exps[..D::COUNT].iter().zip(D::SYMBOLS) {
        match e {
            0 => {}
            1 => write!(f, " {sym}")?,
            _ => write!(f, " {sym}^{e}")?,
        }
    }
    Ok(())
}

impl<D, N: Scalar + PartialOrd> PartialOrd for Quantity<D, N> {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        self.c.partial_cmp(&other.c)
    }
}
impl<D, S: Slots, Sym, K, N: Scalar> PartialEq for Tensor<D, S, Sym, K, N> {
    fn eq(&self, other: &Self) -> bool {
        (0..S::LEN).all(|i| self.flat(i) == other.flat(i))
    }
}
