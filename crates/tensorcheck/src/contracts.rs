//! Shared contracts: the vocabulary types that every layer may use (DD-01).
//!
//! These are type definitions without behavior beyond trivial accessors: spaces,
//! charts, index slots, labels, points and markers. Keeping them here lets the
//! engines share a vocabulary without calling each other (01 §1.7).

use core::marker::PhantomData;

use crate::access::component::{Extent, SlotExtent};
use crate::utility::diagnostics::CoordOf;
use crate::utility::scalar::Scalar;
use crate::utility::typelevel::{Nat, NatEq};

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
