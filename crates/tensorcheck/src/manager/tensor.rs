//! Tensors and quantities, and the operations that assemble expressions (CMP-M01, DD-05).

use core::fmt;
use core::marker::PhantomData;
use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::access::component::Extent;
use crate::access::component::SlotExtent;
use crate::access::component::Slots;
use crate::access::component::{anti_index, sym_index};
use crate::contracts::slot;
use crate::engine::dimension::{DimDiv, DimInfo, DimMul, Prod, Quot, Unit, from_base, to_base};
use crate::engine::dimension::{DimEq, DimPow, Pow};
use crate::engine::index::{
    AllAbsent, At, Chart, Contract, Contracted, Down, IntoLabels, KindMul, LabelIds, Slot, Space,
    StripLabels, Tensorial, Up, ZipLabels, contract,
};
use crate::engine::index::{
    ChartTransform, MAX_RANK, MapChart, Point, SlotAt, SlotVariances, apply_at_slot,
    inverse_jacobian_at, lower, upper,
};
use crate::engine::index::{Here, There};
use crate::engine::metric::{Metric, inverse_at, metric_at};
use crate::engine::symmetry::{Antisymmetric, NoSym, Symmetric};
use crate::utility::diagnostics::{ComponentDimIs, SameSlot};
use crate::utility::diagnostics::{CoordOf, IsTensor, LabelCount, MetricOn, VarianceIs};
use crate::utility::diagnostics::{SameDimension, SameKind, SameSlots, SameSymmetry, UnitMatches};
use crate::utility::linear::LinearError;
use crate::utility::scalar::Scalar;
use crate::utility::typelevel::N1;
use crate::utility::typelevel::{Concat, Cons, Nil};
use crate::utility::vector_space::{PerDim, TimesDim, VectorSpace};

pub use crate::contracts::{Covector, Quantity, Rank2, Rank3, Rank4, Tensor, Vector};

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

impl Quantity<()> {
    /// Creates a quantity whose dimension is the dimension of `unit` (no type
    /// annotation needed). Same as `value * unit`.
    pub fn of<U: Unit>(value: f64, _unit: U) -> Quantity<U::Dim> {
        Tensor::from_store(to_base::<U>(value))
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

// ---------------------------------------------------------------------------
// Symmetric and antisymmetric rank-2 tensors (API-22, API-33, DD-07)

/// A symmetry kind with independent components ([`Symmetric`] or [`Antisymmetric`]).
pub trait SymmetryKind {
    /// The array of independent components for a space with extent `E`.
    type Independent<E: Extent>;
    /// Component `(i, j)` of the full `n × n` matrix from the independent ones.
    fn full<E: Extent>(c: &Self::Independent<E>, i: usize, j: usize) -> f64;
}
impl SymmetryKind for Symmetric {
    type Independent<E: Extent> = E::SymArr<f64>;
    fn full<E: Extent>(c: &E::SymArr<f64>, i: usize, j: usize) -> f64 {
        c.as_ref()[sym_index(i, j, E::N)]
    }
}
impl SymmetryKind for Antisymmetric {
    type Independent<E: Extent> = E::AntiArr<f64>;
    fn full<E: Extent>(c: &E::AntiArr<f64>, i: usize, j: usize) -> f64 {
        match i.cmp(&j) {
            core::cmp::Ordering::Equal => 0.0,
            core::cmp::Ordering::Less => c.as_ref()[anti_index(i, j, E::N)],
            core::cmp::Ordering::Greater => -c.as_ref()[anti_index(j, i, E::N)],
        }
    }
}

impl<D, A, B, Sym: SymmetryKind, K> Tensor<D, Cons<A, Cons<B, Nil>>, Sym, K, f64>
where
    A: SlotExtent + 'static + SameSlot<B>,
    B: SlotExtent + 'static,
{
    /// Creates a symmetric or antisymmetric tensor from its independent components
    /// in `unit` (API-22): `T_ab` with `a <= b` (symmetric) or `a < b`
    /// (antisymmetric), row by row. The rest is filled by the symmetry, so the
    /// result cannot violate it.
    pub fn from_independent<U: Unit + UnitMatches<D>>(
        c: Sym::Independent<A::Ext>,
        _unit: U,
    ) -> Self {
        let n = <A::Ext as Extent>::N;
        Tensor::from_flat_fn(|f| to_base::<U>(Sym::full::<A::Ext>(&c, f / n, f % n)))
    }
}

impl<D, A, B, K, N: Scalar> Tensor<D, Cons<A, Cons<B, Nil>>, NoSym, K, N>
where
    A: SlotExtent + 'static + SameSlot<B>,
    B: SlotExtent + 'static,
{
    /// The symmetric part `(T_ab + T_ba)/2` (API-33).
    pub fn symmetrize(self) -> Tensor<D, Cons<A, Cons<B, Nil>>, Symmetric, K, N> {
        let n = <A::Ext as Extent>::N;
        let half = N::from_f64(0.5);
        Tensor::from_flat_fn(|f| (self.flat(f) + self.flat((f % n) * n + f / n)) * half)
    }
    /// The antisymmetric part `(T_ab − T_ba)/2` (API-33).
    pub fn antisymmetrize(self) -> Tensor<D, Cons<A, Cons<B, Nil>>, Antisymmetric, K, N> {
        let n = <A::Ext as Extent>::N;
        let half = N::from_f64(0.5);
        Tensor::from_flat_fn(|f| (self.flat(f) - self.flat((f % n) * n + f / n)) * half)
    }
}

impl<D, S: Slots, Sym, K, N: Scalar> Tensor<D, S, Sym, K, N> {
    /// Forgets the symmetry marker (the components are unchanged).
    pub fn without_symmetry(self) -> Tensor<D, S, NoSym, K, N> {
        Tensor::from_store(self.c)
    }
}

// ---------------------------------------------------------------------------
// The metric tensor as a tensor at a point

/// The metric tensor `g_ab` of `g` at `p`, as a symmetric tensor.
pub fn metric_tensor<G: Metric>(
    g: &G,
    p: &Point<G::Chart>,
) -> Tensor<G::Dim, Cons<Down<G::Chart>, Cons<Down<G::Chart>, Nil>>, Symmetric> {
    let m = metric_at(g, p);
    Tensor::from_flat_fn(|f| m[f])
}

/// The inverse metric tensor `g^ab` of `g` at `p`, as a symmetric tensor of
/// dimension `1/[g]`. Fails if the metric tensor is singular at `p`.
pub fn inverse_metric_tensor<G: Metric>(
    g: &G,
    p: &Point<G::Chart>,
) -> Result<
    Tensor<Pow<G::Dim, N1>, Cons<Up<G::Chart>, Cons<Up<G::Chart>, Nil>>, Symmetric>,
    LinearError,
>
where
    G::Dim: DimPow<N1>,
{
    let m = inverse_at(g, p)?;
    Ok(Tensor::from_flat_fn(|f| m[f]))
}

// ---------------------------------------------------------------------------
// Vectors from per-component quantities (API-34)

/// A tuple of quantities that can be the components of a vector of dimension `D`
/// in chart `C`.
pub trait VectorParts<D, C: Chart> {
    /// Writes the components in base units.
    fn write(&self, out: &mut [f64]);
}
macro_rules! vector_parts {
    ($($q:ident $i:tt $pos:ty, $p:ty);+) => {
        impl<D, C: Chart, $($q),+> VectorParts<D, C> for ($(Quantity<$q>,)+)
        where
            $(
                C::CoordDims: At<$pos>,
                Up<C>: VarianceDim<D, <C::CoordDims as At<$pos>>::Out>,
                $q: DimEq<<Up<C> as VarianceDim<D, <C::CoordDims as At<$pos>>::Out>>::Out>,
                <$q as DimEq<<Up<C> as VarianceDim<D, <C::CoordDims as At<$pos>>::Out>>::Out>>::Out: ComponentDimIs<$p>,
            )+
        {
            fn write(&self, out: &mut [f64]) {
                $( out[$i] = self.$i.c; )+
            }
        }
    };
}
vector_parts!(Q0 0 Here, slot::_0);
vector_parts!(Q0 0 Here, slot::_0; Q1 1 There<Here>, slot::_1);
vector_parts!(Q0 0 Here, slot::_0; Q1 1 There<Here>, slot::_1; Q2 2 There<There<Here>>, slot::_2);
vector_parts!(Q0 0 Here, slot::_0; Q1 1 There<Here>, slot::_1; Q2 2 There<There<Here>>, slot::_2; Q3 3 There<There<There<Here>>>, slot::_3);

impl<D, C: Chart> Tensor<D, Cons<Up<C>, Nil>, NoSym, Tensorial, f64> {
    /// A vector from its components as quantities, each with its own dimension
    /// (API-34). In polar coordinates a velocity is
    /// `Vector::from_parts((1.0 * (m / s), 0.25 * (rad / s)))`.
    pub fn from_parts<P: VectorParts<D, C>>(parts: P) -> Self {
        let mut buf = [0.0; 16];
        parts.write(&mut buf);
        Tensor::from_flat_fn(|i| buf[i])
    }
}

// ---------------------------------------------------------------------------
// Tensors and points as solver states (U05, DD-11)

impl<D, S: Slots, Sym, K> VectorSpace for Tensor<D, S, Sym, K, f64> {
    const LEN: usize = S::LEN;
    fn component(&self, i: usize) -> f64 {
        self.flat(i)
    }
    fn from_components(f: impl FnMut(usize) -> f64) -> Self {
        Tensor::from_flat_fn(f)
    }
}
impl<D: DimDiv<Dt>, Dt, S: Slots, Sym, K> PerDim<Dt> for Tensor<D, S, Sym, K, f64> {
    type Out = Tensor<Quot<D, Dt>, S, Sym, K, f64>;
}
impl<D: DimMul<Dx>, Dx, S: Slots, Sym, K> TimesDim<Dx> for Tensor<D, S, Sym, K, f64> {
    type Out = Tensor<Prod<D, Dx>, S, Sym, K, f64>;
}

impl<C: Chart> VectorSpace for Point<C> {
    const LEN: usize = <C::Space as Extent>::N;
    fn component(&self, i: usize) -> f64 {
        self.get(i)
    }
    fn from_components(f: impl FnMut(usize) -> f64) -> Self {
        Point::from_base(C::Space::from_fn(f))
    }
}
/// The rate of change of a point is a vector of coordinate velocities, with overall
/// dimension `[X]/[D]` (components `[x^a]/[D]`).
impl<C: Chart, Dt> PerDim<Dt> for Point<C>
where
    <C::Space as Space>::RefDim: DimDiv<Dt>,
{
    type Out = Vector<Quot<<C::Space as Space>::RefDim, Dt>, C>;
}
