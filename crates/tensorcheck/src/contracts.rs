//! Shared contracts: the vocabulary types that every layer may use (DD-01).
//!
//! These are type definitions without behavior beyond trivial accessors: spaces,
//! charts, index slots, labels, points and markers. Keeping them here lets the
//! engines share a vocabulary without calling each other (01 §1.7).

use core::marker::PhantomData;

use crate::access::component::{Extent, SlotExtent, Slots};
use crate::utility::diagnostics::CoordOf;
use crate::utility::scalar::Scalar;
use crate::utility::typelevel::{Cons, Nat, NatEq, Nil};

/// Marks a tensor: it obeys the tensor transformation law.
#[derive(Clone, Copy, Debug)]
pub struct Tensorial;
/// Marks a chart-bound quantity (not a tensor), such as Christoffel symbols.
#[derive(Clone, Copy, Debug)]
pub struct ChartBound;

/// A space that indices range over (EXT-03): spatial, spacetime, field space, ...
pub trait Space: Extent {
    /// The reference dimension `[X]` of coordinates on this space (length for
    /// physical space and spacetime, the field dimension for a field space).
    ///
    /// A component's dimension is the tensor's dimension times `[x^a]/[X]` for each
    /// upper index and divided by it for each lower index, so in a chart whose
    /// coordinates all have dimension `[X]` (Cartesian) the components have the
    /// tensor's own dimension.
    type RefDim;
    /// Name used in documentation and messages.
    const NAME: &'static str;
}

/// The basis of a chart is the coordinate basis `∂/∂x^a`.
#[derive(Clone, Copy, Debug)]
pub struct CoordinateBasis;
/// The basis of a chart is not a coordinate basis (orthonormal frame, tetrad).
#[derive(Clone, Copy, Debug)]
pub struct NonCoordinateBasis;

/// A chart (coordinate system and basis) on a space (EXT-04).
pub trait Chart: 'static {
    /// The space this chart belongs to.
    type Space: Space;
    /// The coordinate markers, as a type-level list.
    type Coords;
    /// The dimension of each coordinate, as a type-level list.
    type CoordDims;
    /// [`CoordinateBasis`] or [`NonCoordinateBasis`].
    type Basis;
    /// Name used in documentation and messages.
    const NAME: &'static str;
}

/// Variance marker: an upper (contravariant) index.
#[allow(non_camel_case_types)]
pub struct upper;
/// Variance marker: a lower (covariant) index.
#[allow(non_camel_case_types)]
pub struct lower;

/// An upper index in chart `C`.
pub struct Up<C>(PhantomData<C>);
/// A lower index in chart `C`.
pub struct Down<C>(PhantomData<C>);

/// An index slot: chart and variance.
pub trait Slot: 'static {
    /// The chart.
    type Chart: Chart;
    /// [`upper`] or [`lower`].
    type Var;
    /// Whether the index is upper.
    const UP: bool;
    /// The same slot with the opposite variance.
    type Flip: Slot<Chart = Self::Chart>;
}
impl<C: Chart> Slot for Up<C> {
    type Chart = C;
    type Var = upper;
    const UP: bool = true;
    type Flip = Down<C>;
}
impl<C: Chart> Slot for Down<C> {
    type Chart = C;
    type Var = lower;
    const UP: bool = false;
    type Flip = Up<C>;
}
impl<C: Chart> SlotExtent for Up<C> {
    type Ext = C::Space;
}
impl<C: Chart> SlotExtent for Down<C> {
    type Ext = C::Space;
}

/// An index label such as `a` or `mu` (API-08). Labels are compared by `Id`.
pub trait Label: Copy + 'static {
    /// The label's number (a type-level natural). Distinct labels need distinct numbers.
    type Id: Nat + NatEq<Self::Id>;
    /// The label's name.
    const NAME: &'static str;
}

/// A slot `Sl` carrying label `L`.
pub struct Ix<L, Sl>(PhantomData<(L, Sl)>);

/// A point given by its coordinates in chart `C`, in base units (API-10).
pub struct Point<C: Chart, N: Scalar = f64> {
    pub(crate) x: <C::Space as Extent>::Arr<N>,
}
impl<C: Chart, N: Scalar> Clone for Point<C, N> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<C: Chart, N: Scalar> Copy for Point<C, N> {}
impl<C: Chart, N: Scalar> core::fmt::Debug for Point<C, N> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}{:?}", C::NAME, self.x.as_ref())
    }
}

impl<C: Chart, N: Scalar> Point<C, N> {
    /// A point from coordinates in base units, in the chart's coordinate order.
    pub fn from_base(x: <C::Space as Extent>::Arr<N>) -> Self {
        Point { x }
    }
    /// Coordinate number `i`, in base units.
    pub fn get(&self, i: usize) -> N {
        self.x.as_ref()[i]
    }
    /// The coordinate marked by `X` (for example `p.coord(Theta)`), in base units.
    pub fn coord<X, I>(&self, _x: X) -> N
    where
        C::Coords: CoordOf<X, C, I>,
    {
        self.x.as_ref()[<C::Coords as CoordOf<X, C, I>>::POS]
    }
    /// All coordinates, in base units.
    pub fn coords(&self) -> &[N] {
        self.x.as_ref()
    }
    /// The same point with coordinates converted to another scalar type.
    pub fn map<M: Scalar>(&self, f: impl Fn(N) -> M) -> Point<C, M> {
        Point {
            x: C::Space::from_fn(|i| f(self.x.as_ref()[i])),
        }
    }
}

/// Slot-position markers for `raise` and `lower` (API-09).
pub mod slot {
    /// The first index.
    #[allow(non_camel_case_types)]
    #[derive(Clone, Copy, Debug)]
    pub struct _0;
    /// The second index.
    #[allow(non_camel_case_types)]
    #[derive(Clone, Copy, Debug)]
    pub struct _1;
    /// The third index.
    #[allow(non_camel_case_types)]
    #[derive(Clone, Copy, Debug)]
    pub struct _2;
    /// The fourth index.
    #[allow(non_camel_case_types)]
    #[derive(Clone, Copy, Debug)]
    pub struct _3;
}

// ---------------------------------------------------------------------------
// Symmetry markers (E04 decides the rules)

/// No symmetry.
#[derive(Clone, Copy, Debug)]
pub struct NoSym;
/// Symmetric in its two indices: `T_ab = T_ba`.
#[derive(Clone, Copy, Debug)]
pub struct Symmetric;
/// Antisymmetric in its two indices: `T_ab = -T_ba`.
#[derive(Clone, Copy, Debug)]
pub struct Antisymmetric;

// ---------------------------------------------------------------------------
// Tensors (M01 and M03 build and combine them)

/// A tensor with dimension `D`, index slots `S`, symmetry `Sym`, tensoriality `K`
/// and component type `N` (API-01).
///
/// Components are stored in the base unit of the unit system. A quantity is a
/// tensor of rank 0 ([`Quantity`]).
pub struct Tensor<D, S: Slots, Sym = NoSym, K = Tensorial, N: Scalar = f64> {
    pub(crate) c: S::Store<N>,
    pub(crate) _m: PhantomData<(D, Sym, K)>,
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

impl<D, S: Slots, N: Scalar> Tensor<D, S, NoSym, Tensorial, N> {
    /// A tensor from its components in base units, for any scalar type `N`. Use it in
    /// tensor fields written generically over `N` (EXT-08); prefer
    /// `from_components` with an explicit unit elsewhere.
    pub fn from_base_components(c: S::Store<N>) -> Self {
        Tensor { c, _m: PhantomData }
    }
}

impl<D, S: Slots, Sym, K, N: Scalar> Tensor<D, S, Sym, K, N> {
    /// The components in base units, as nested arrays.
    pub fn base_components(&self) -> S::Store<N> {
        self.c
    }
}
