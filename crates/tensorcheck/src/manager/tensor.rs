//! Tensors and quantities, and the operations that assemble expressions (CMP-M01, DD-05).

use core::fmt;
use core::marker::PhantomData;
use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::access::component::Extent;
use crate::access::component::SlotExtent;
use crate::access::component::Slots;
use crate::engine::dimension::{DimDiv, DimInfo, DimMul, Prod, Quot, Unit, from_base, to_base};
use crate::engine::index::{
    AllAbsent, At, Chart, Contract, Contracted, Down, IntoLabels, KindMul, LabelIds, Slot, Space,
    StripLabels, Tensorial, Up, ZipLabels, contract,
};
use crate::engine::index::{
    ChartTransform, MAX_RANK, MapChart, Point, SlotAt, SlotVariances, apply_at_slot,
    inverse_jacobian_at, lower, upper,
};
use crate::engine::metric::{Metric, inverse_at, metric_at};
use crate::engine::symmetry::NoSym;
use crate::utility::diagnostics::{CoordOf, IsTensor, LabelCount, MetricOn, VarianceIs};
use crate::utility::diagnostics::{SameDimension, SameKind, SameSlots, SameSymmetry, UnitMatches};
use crate::utility::linear::LinearError;
use crate::utility::scalar::Scalar;
use crate::utility::typelevel::{Concat, Cons, Nil};

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
/// A vector (one upper index) in chart `C` (API-02).
pub type Vector<D, C, N = f64> = Tensor<D, Cons<Up<C>, Nil>, NoSym, Tensorial, N>;
/// A covector (one lower index) in chart `C` (API-02).
pub type Covector<D, C, N = f64> = Tensor<D, Cons<Down<C>, Nil>, NoSym, Tensorial, N>;
/// A rank-2 tensor with slots `A`, `B` (each `Up<C>` or `Down<C>`).
pub type Rank2<D, A, B, Sym = NoSym, K = Tensorial, N = f64> =
    Tensor<D, Cons<A, Cons<B, Nil>>, Sym, K, N>;
/// A rank-3 tensor with slots `A`, `B`, `C`.
pub type Rank3<D, A, B, C, Sym = NoSym, K = Tensorial, N = f64> =
    Tensor<D, Cons<A, Cons<B, Cons<C, Nil>>>, Sym, K, N>;
/// A rank-4 tensor with slots `A`, `B`, `C`, `E`.
pub type Rank4<D, A, B, C, E, Sym = NoSym, K = Tensorial, N = f64> =
    Tensor<D, Cons<A, Cons<B, Cons<C, Cons<E, Nil>>>>, Sym, K, N>;

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
// Addition and subtraction: same dimension (ERR-01), slots (ERR-03), symmetry
// (ERR-13) and tensoriality (ERR-17).

macro_rules! add_like {
    ($tr:ident, $f:ident, $op:tt) => {
        impl<D1, D2, S1, S2, Sym1, Sym2, K1, K2, N> $tr<Tensor<D2, S2, Sym2, K2, N>> for Tensor<D1, S1, Sym1, K1, N>
        where
            S1: Slots + SameSlots<S2>,
            S2: Slots,
            N: Scalar,
            D1: SameDimension<D2>,
            Sym1: SameSymmetry<Sym2>,
            K1: SameKind<K2>,
        {
            type Output = Tensor<D1, S1, Sym1, K1, N>;
            #[inline(always)]
            fn $f(self, rhs: Tensor<D2, S2, Sym2, K2, N>) -> Self::Output {
                Tensor::from_flat_fn(|i| self.flat(i) $op rhs.flat(i))
            }
        }
    };
}
add_like!(Add, add, +);
add_like!(Sub, sub, -);

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

// ---------------------------------------------------------------------------
// Components with units (API-21, API-23)

impl<D, S: Slots, Sym, K> Tensor<D, S, Sym, K, f64> {
    /// Creates a tensor from its components in `unit`, as nested arrays in index
    /// order (API-21). For symmetric tensors use `from_independent`.
    ///
    /// ```
    /// use tensorcheck::Vector;
    /// use tensorcheck::spaces::euclid3::Cartesian;
    /// use tensorcheck::units::si::{m, s, Velocity};
    /// let v: Vector<Velocity, Cartesian> = Vector::from_components([1.0, 2.0, 3.0], m / s);
    /// assert_eq!(v.components_in(m / s), [1.0, 2.0, 3.0]);
    /// ```
    pub fn from_components<U: Unit + UnitMatches<D>>(c: S::Store<f64>, _unit: U) -> Self
    where
        Sym: NonSymmetric,
    {
        let t = Self::from_store(c);
        t.map(to_base::<U>)
    }

    /// The components in `unit`, as nested arrays (API-23).
    pub fn components_in<U: Unit + UnitMatches<D>>(self, _unit: U) -> S::Store<f64> {
        self.map(from_base::<U>).c
    }
}

/// Symmetry markers that place no constraint on components.
pub trait NonSymmetric {}
impl NonSymmetric for NoSym {}

// ---------------------------------------------------------------------------
// Outer product (API-26)

impl<D, S: Slots, Sym, K, N: Scalar> Tensor<D, S, Sym, K, N> {
    /// The outer product `A ⊗ B`; the slots of `B` follow those of `A` (API-26).
    pub fn outer<D2, S2: Slots, Sym2, K2>(
        self,
        other: Tensor<D2, S2, Sym2, K2, N>,
    ) -> Tensor<Prod<D, D2>, <S as Concat<S2>>::Out, NoSym, <K as KindMul<K2>>::Out, N>
    where
        D: DimMul<D2>,
        S: Concat<S2>,
        <S as Concat<S2>>::Out: Slots,
        K: KindMul<K2>,
    {
        Tensor::from_flat_fn(|i| self.flat(i / S2::LEN) * other.flat(i % S2::LEN))
    }

    /// Attaches labels to the indices for Einstein-convention expressions (API-27).
    ///
    /// ```
    /// use tensorcheck::labels::{a, b};
    /// use tensorcheck::spaces::euclid3::Cartesian;
    /// use tensorcheck::units::si::{m, one};
    /// use tensorcheck::{Tensor, Up, Down, Vector};
    /// use tensorcheck::utility::typelevel::{Cons, Nil};
    /// use tensorcheck::units::si::{Dimensionless, Length};
    /// let t: Tensor<Dimensionless, Cons<Up<Cartesian>, Cons<Down<Cartesian>, Nil>>> =
    ///     Tensor::from_components([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]], one);
    /// let v: Vector<Length, Cartesian> = Vector::from_components([1.0, 1.0, 1.0], m);
    /// let w = (t.idx((a, b)) * v.idx(b)).into_tensor();
    /// assert_eq!(w.components_in(m), [1.0, 2.0, 3.0]);
    /// ```
    pub fn idx<L: IntoLabels>(self, _labels: L) -> Labeled<Self, <S as ZipLabels<L::List>>::Out>
    where
        S: ZipLabels<L::List>,
        <S as ZipLabels<L::List>>::Equal: LabelCount<L::List>,
    {
        Labeled {
            t: self,
            _l: PhantomData,
        }
    }
}

/// A tensor with labels on its indices (API-27). Multiply labeled tensors to
/// contract repeated labels; call [`Labeled::into_tensor`] to get the result.
///
/// `U` lists the labels already contracted inside this expression; using one of
/// them again is an error (ERR-06).
pub struct Labeled<T, L, U = Nil> {
    t: T,
    _l: PhantomData<(L, U)>,
}

impl<D, S: Slots, Sym, K, N: Scalar, L, U> Labeled<Tensor<D, S, Sym, K, N>, L, U> {
    /// Contracts repeated labels within this expression and returns the tensor. The
    /// remaining indices keep their order of appearance (API-29, API-30).
    ///
    /// The result has no symmetry marker (2b §CMP-E04).
    pub fn into_tensor(self) -> Tensor<D, <<L as Contract>::Out as StripLabels>::Slots, NoSym, K, N>
    where
        L: Contract + LabelIds,
        <L as Contract>::Out: StripLabels,
    {
        let one = N::one();
        Tensor::from_store(contract::<
            S,
            Nil,
            <<L as Contract>::Out as StripLabels>::Slots,
            L,
            Nil,
            N,
        >(&self.t.c, &one))
    }
}

impl<D1, S1, Sym1, K1, N, LA, UA, D2, S2, Sym2, K2, LB, UB>
    Mul<Labeled<Tensor<D2, S2, Sym2, K2, N>, LB, UB>>
    for Labeled<Tensor<D1, S1, Sym1, K1, N>, LA, UA>
where
    S1: Slots,
    S2: Slots,
    N: Scalar,
    D1: DimMul<D2>,
    K1: KindMul<K2>,
    LA: Concat<LB> + LabelIds,
    LB: LabelIds,
    <LA as Concat<LB>>::Out: Contract + AllAbsent<<UA as Concat<UB>>::Out>,
    Contracted<LA, LB>: StripLabels,
    UA: Concat<UB>,
    <UA as Concat<UB>>::Out: Concat<<<LA as Concat<LB>>::Out as Contract>::Dummies>,
{
    type Output = Labeled<
        Tensor<
            Prod<D1, D2>,
            <Contracted<LA, LB> as StripLabels>::Slots,
            NoSym,
            <K1 as KindMul<K2>>::Out,
            N,
        >,
        Contracted<LA, LB>,
        <<UA as Concat<UB>>::Out as Concat<<<LA as Concat<LB>>::Out as Contract>::Dummies>>::Out,
    >;
    #[inline(always)]
    fn mul(self, rhs: Labeled<Tensor<D2, S2, Sym2, K2, N>, LB, UB>) -> Self::Output {
        let c = contract::<S1, S2, <Contracted<LA, LB> as StripLabels>::Slots, LA, LB, N>(
            &self.t.c, &rhs.t.c,
        );
        Labeled {
            t: Tensor::from_store(c),
            _l: PhantomData,
        }
    }
}

// ---------------------------------------------------------------------------
// One component with its own dimension (API-24)

/// The dimension of the component selected by coordinates `X` (DD-04).
pub trait ComponentDim<D, X, I> {
    /// `[T] ∏ [x^up] / ∏ [x^down]`.
    type Out;
    /// The flat index of the component.
    fn flat_index() -> usize;
}
impl<D> ComponentDim<D, Nil, Nil> for Nil {
    type Out = D;
    fn flat_index() -> usize {
        0
    }
}
impl<D, Sl, S, X, Xs, I, Is> ComponentDim<D, Cons<X, Xs>, Cons<I, Is>> for Cons<Sl, S>
where
    Sl: Slot + SlotExtent,
    S: Slots + ComponentDim<D, Xs, Is>,
    <Sl::Chart as Chart>::Coords: CoordOf<X, Sl::Chart, I>,
    <Sl::Chart as Chart>::CoordDims: At<I>,
    Sl: VarianceDim<
            <S as ComponentDim<D, Xs, Is>>::Out,
            <<Sl::Chart as Chart>::CoordDims as At<I>>::Out,
        >,
{
    type Out = <Sl as VarianceDim<
        <S as ComponentDim<D, Xs, Is>>::Out,
        <<Sl::Chart as Chart>::CoordDims as At<I>>::Out,
    >>::Out;
    fn flat_index() -> usize {
        <<Sl::Chart as Chart>::Coords as CoordOf<X, Sl::Chart, I>>::POS * S::LEN + S::flat_index()
    }
}

/// Multiplies (upper) or divides (lower) by `[x^a]/[X]`, the coordinate's dimension
/// relative to the space's reference dimension.
pub trait VarianceDim<D, Dx> {
    /// The resulting dimension.
    type Out;
}
impl<C: Chart, D, Dx> VarianceDim<D, Dx> for Up<C>
where
    Dx: DimDiv<<C::Space as Space>::RefDim>,
    D: DimMul<Quot<Dx, <C::Space as Space>::RefDim>>,
{
    type Out = Prod<D, Quot<Dx, <C::Space as Space>::RefDim>>;
}
impl<C: Chart, D, Dx> VarianceDim<D, Dx> for Down<C>
where
    Dx: DimDiv<<C::Space as Space>::RefDim>,
    D: DimDiv<Quot<Dx, <C::Space as Space>::RefDim>>,
{
    type Out = Quot<D, Quot<Dx, <C::Space as Space>::RefDim>>;
}

/// Coordinates passed to `at`: a tuple of coordinate markers.
pub trait IntoCoords {
    /// The coordinates as a type-level list.
    type List;
}
macro_rules! into_coords_tuple {
    ($($x:ident),+) => {
        impl<$($x),+> IntoCoords for ($($x,)+) {
            type List = into_coords_tuple!(@list $($x),+);
        }
    };
    (@list $h:ident $(, $t:ident)*) => { Cons<$h, into_coords_tuple!(@list $($t),*)> };
    (@list) => { Nil };
}
into_coords_tuple!(X1);
into_coords_tuple!(X1, X2);
into_coords_tuple!(X1, X2, X3);
into_coords_tuple!(X1, X2, X3, X4);

impl<D, S: Slots, Sym, K, N: Scalar> Tensor<D, S, Sym, K, N> {
    /// The component selected by coordinate markers, with its own dimension
    /// (API-24). For example `v.at((Theta,))` in polar coordinates has dimension
    /// `[v]·[θ]`.
    pub fn at<X: IntoCoords, I>(
        self,
        _coords: X,
    ) -> Quantity<<S as ComponentDim<D, X::List, I>>::Out, N>
    where
        S: ComponentDim<D, X::List, I>,
    {
        Tensor::from_store(self.flat(<S as ComponentDim<D, X::List, I>>::flat_index()))
    }
}

// ---------------------------------------------------------------------------
// Points from quantities (API-10)

/// Coordinates passed to `Point::new`: a tuple of quantities whose dimensions are
/// the chart's coordinate dimensions, in order.
pub trait CoordQuantities<C: Chart> {
    /// Writes the values in base units.
    fn write(&self, out: &mut [f64]);
}
macro_rules! coord_quantities {
    ($($q:ident $i:tt),+) => {
        impl<C, $($q),+> CoordQuantities<C> for ($(Quantity<$q>,)+)
        where
            C: Chart<CoordDims = coord_quantities!(@list $($q),+)>,
        {
            fn write(&self, out: &mut [f64]) {
                $( out[$i] = self.$i.c; )+
            }
        }
    };
    (@list $h:ident $(, $t:ident)*) => { Cons<$h, coord_quantities!(@list $($t),*)> };
    (@list) => { Nil };
}
coord_quantities!(Q0 0);
coord_quantities!(Q0 0, Q1 1);
coord_quantities!(Q0 0, Q1 1, Q2 2);
coord_quantities!(Q0 0, Q1 1, Q2 2, Q3 3);
coord_quantities!(Q0 0, Q1 1, Q2 2, Q3 3, Q4 4);

impl<C: Chart> Point<C> {
    /// A point from its coordinates as quantities, in the chart's order (API-10).
    ///
    /// ```
    /// use tensorcheck::Point;
    /// use tensorcheck::spaces::plane::{Polar, Theta};
    /// use tensorcheck::units::si::{m, rad};
    /// let p = Point::<Polar>::new((2.0 * m, 0.5 * rad));
    /// assert_eq!(p.coord(Theta), 0.5);
    /// ```
    pub fn new<Q: CoordQuantities<C>>(coords: Q) -> Self {
        let mut buf = [0.0; 16];
        coords.write(&mut buf);
        Point::from_base(C::Space::from_fn(|i| buf[i]))
    }
}

// ---------------------------------------------------------------------------
// Raising, lowering and changing charts (API-31, API-32)

impl<D, S: Slots, Sym, K> Tensor<D, S, Sym, K, f64> {
    fn flat_vec(&self) -> (Vec<f64>, Vec<usize>) {
        (
            (0..S::LEN).map(|i| self.flat(i)).collect(),
            (0..S::RANK).map(S::extent).collect(),
        )
    }

    /// Lowers index number `P` with the metric tensor `g` at `p`:
    /// `T_a = g_ab T^b` (API-31).
    pub fn lower<P, G: Metric>(
        self,
        _slot: P,
        g: &G,
        p: &Point<G::Chart>,
    ) -> Tensor<Prod<D, G::Dim>, <S as SlotAt<P>>::Flipped, NoSym, K>
    where
        S: SlotAt<P>,
        <S as SlotAt<P>>::Flipped: Slots,
        <<S as SlotAt<P>>::Slot as Slot>::Var: VarianceIs<upper, P>,
        G::Chart: MetricOn<<<S as SlotAt<P>>::Slot as Slot>::Chart, G>,
        D: DimMul<G::Dim>,
    {
        let gm = metric_at(g, p);
        let (flat, ext) = self.flat_vec();
        let out = apply_at_slot(&flat, &ext, <S as SlotAt<P>>::POS, &gm);
        Tensor::from_flat_fn(|i| out[i])
    }

    /// Raises index number `P` with the inverse metric tensor at `p`:
    /// `T^a = g^ab T_b` (API-31). Fails if the metric tensor is singular at `p`.
    pub fn raise<P, G: Metric>(
        self,
        _slot: P,
        g: &G,
        p: &Point<G::Chart>,
    ) -> Result<Tensor<Quot<D, G::Dim>, <S as SlotAt<P>>::Flipped, NoSym, K>, LinearError>
    where
        S: SlotAt<P>,
        <S as SlotAt<P>>::Flipped: Slots,
        <<S as SlotAt<P>>::Slot as Slot>::Var: VarianceIs<lower, P>,
        G::Chart: MetricOn<<<S as SlotAt<P>>::Slot as Slot>::Chart, G>,
        D: DimDiv<G::Dim>,
    {
        let ginv = inverse_at(g, p)?;
        let (flat, ext) = self.flat_vec();
        let out = apply_at_slot(&flat, &ext, <S as SlotAt<P>>::POS, &ginv);
        Ok(Tensor::from_flat_fn(|i| out[i]))
    }

    /// Expresses this tensor at `p` in the target chart of `xf` (API-32). Upper
    /// indices transform with `∂x'^a/∂x^b`, lower indices with `∂x^b/∂x'^a`. Every
    /// index must be in the source chart (ERR-10), and only tensors can be
    /// transformed (ERR-16). Fails if the Jacobian is singular at `p`.
    pub fn transform<X: ChartTransform>(
        self,
        xf: &X,
        p: &Point<X::From>,
    ) -> Result<Tensor<D, <S as MapChart<X::From, X::To, X>>::Out, Sym, K>, LinearError>
    where
        S: MapChart<X::From, X::To, X> + SlotVariances,
        <S as MapChart<X::From, X::To, X>>::Out: Slots,
        K: IsTensor,
    {
        let n = <<X::From as Chart>::Space as Extent>::N;
        let jac = xf.jacobian(p);
        let mut ups = [false; MAX_RANK];
        S::write(&mut ups[..S::RANK]);
        let (mut data, ext) = self.flat_vec();
        let lower_mat = if ups[..S::RANK].iter().any(|u| !u) {
            // M[a'][b] = ∂x^b/∂x'^a, the transposed inverse Jacobian.
            let inv = inverse_jacobian_at(xf, p)?;
            let mut m = vec![0.0; n * n];
            for a in 0..n {
                for b in 0..n {
                    m[a * n + b] = inv[b * n + a];
                }
            }
            m
        } else {
            Vec::new()
        };
        for (slot, up) in ups[..S::RANK].iter().enumerate() {
            data = apply_at_slot(&data, &ext, slot, if *up { &jac } else { &lower_mat });
        }
        Ok(Tensor::from_flat_fn(|i| data[i]))
    }
}
