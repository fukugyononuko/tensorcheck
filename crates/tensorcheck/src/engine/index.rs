//! Index slots, spaces, charts, labels and contraction (CMP-E02, DD-04, DD-06).
//!
//! Type level: which indices can be contracted and what slots remain. Value level:
//! the contraction sums, computed directly without forming outer products.

use core::marker::PhantomData;

use crate::access::component::{Extent, SlotExtent, Slots};
use crate::utility::diagnostics::{
    AppearsOnce, CoordOf, LabelCount, Opposite, SameChart, SameSpace,
};
use crate::utility::scalar::Scalar;
use crate::utility::typelevel::{Bool, Concat, Cons, False, Nat, NatEq, Nil, True};

// ---------------------------------------------------------------------------
// Tensoriality

/// Marks a tensor: it obeys the tensor transformation law.
#[derive(Clone, Copy, Debug)]
pub struct Tensorial;
/// Marks a chart-bound quantity (not a tensor), such as Christoffel symbols.
#[derive(Clone, Copy, Debug)]
pub struct ChartBound;

/// Tensoriality of a product: a tensor only if both factors are tensors.
pub trait KindMul<K2> {
    /// The tensoriality of the product.
    type Out;
}
impl KindMul<Tensorial> for Tensorial {
    type Out = Tensorial;
}
impl KindMul<ChartBound> for Tensorial {
    type Out = ChartBound;
}
impl<K> KindMul<K> for ChartBound {
    type Out = ChartBound;
}

// ---------------------------------------------------------------------------
// Spaces, charts and slots

/// A space that indices range over (EXT-03): spatial, spacetime, field space, ...
pub trait Space: Extent {
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
impl<L> Opposite<lower, L> for upper {}
impl<L> Opposite<upper, L> for lower {}

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

// ---------------------------------------------------------------------------
// Labels

/// An index label such as `a` or `mu` (API-08). Labels are compared by `Id`.
pub trait Label: Copy + 'static {
    /// The label's number (a type-level natural). Distinct labels need distinct numbers.
    type Id: Nat + NatEq<Self::Id>;
    /// The label's name.
    const NAME: &'static str;
}

/// A slot `Sl` carrying label `L`.
pub struct Ix<L, Sl>(PhantomData<(L, Sl)>);

/// Labels passed to `idx`: a single label or a tuple of labels.
pub trait IntoLabels {
    /// The labels as a type-level list.
    type List;
}
macro_rules! into_labels_tuple {
    ($($l:ident),+) => {
        impl<$($l: Label),+> IntoLabels for ($($l,)+) {
            type List = into_labels_tuple!(@list $($l),+);
        }
    };
    (@list $h:ident $(, $t:ident)*) => { Cons<$h, into_labels_tuple!(@list $($t),*)> };
    (@list) => { Nil };
}
impl<A1: Label> IntoLabels for A1 {
    type List = Cons<A1, Nil>;
}
into_labels_tuple!(A1);
into_labels_tuple!(A1, A2);
into_labels_tuple!(A1, A2, A3);
into_labels_tuple!(A1, A2, A3, A4);
into_labels_tuple!(A1, A2, A3, A4, A5);
into_labels_tuple!(A1, A2, A3, A4, A5, A6);

/// Zips slots with labels. Defined for lists of any lengths (the extra elements are
/// dropped) so that a length mismatch is reported by [`LabelCount`] alone.
pub trait ZipLabels<L> {
    /// The slots zipped with their labels.
    type Out;
    /// Whether the lengths are equal.
    type Equal: Bool;
}
impl ZipLabels<Nil> for Nil {
    type Out = Nil;
    type Equal = True;
}
impl<L, Ls> ZipLabels<Cons<L, Ls>> for Nil {
    type Out = Nil;
    type Equal = False;
}
impl<Sl, S> ZipLabels<Nil> for Cons<Sl, S> {
    type Out = Nil;
    type Equal = False;
}
impl<Sl, S: ZipLabels<Ls>, L, Ls> ZipLabels<Cons<L, Ls>> for Cons<Sl, S> {
    type Out = Cons<Ix<L, Sl>, <S as ZipLabels<Ls>>::Out>;
    type Equal = <S as ZipLabels<Ls>>::Equal;
}
impl<L> LabelCount<L> for True {}

// ---------------------------------------------------------------------------
// Contraction at the type level (Einstein convention)
//
// The labeled slots of all factors are concatenated; every label that appears twice
// is contracted, after checking space, chart and variance (ERR-07, ERR-08, ERR-05).

/// Searches a labeled list for label `L` (paired with slot `Sl`) and removes the
/// first match.
pub trait FindRemove<L, Sl> {
    /// [`True`] if found.
    type Found: Bool;
    /// The list without the match.
    type Rest;
}
impl<L, Sl> FindRemove<L, Sl> for Nil {
    type Found = False;
    type Rest = Nil;
}
impl<L: Label, Sl, L2: Label, Sl2, T> FindRemove<L, Sl> for Cons<Ix<L2, Sl2>, T>
where
    L2::Id: NatEq<L::Id>,
    Self: FindBranch<<L2::Id as NatEq<L::Id>>::Out, L, Sl>,
{
    type Found = <Self as FindBranch<<L2::Id as NatEq<L::Id>>::Out, L, Sl>>::Found;
    type Rest = <Self as FindBranch<<L2::Id as NatEq<L::Id>>::Out, L, Sl>>::Rest;
}

/// Helper of [`FindRemove`]: `B` tells whether the head matches.
pub trait FindBranch<B, L, Sl> {
    /// [`True`] if found.
    type Found: Bool;
    /// The list without the match.
    type Rest;
}
impl<L, Sl: Slot, L2, Sl2: Slot, T> FindBranch<True, L, Sl> for Cons<Ix<L2, Sl2>, T>
where
    <Sl::Chart as Chart>::Space: SameSpace<<Sl2::Chart as Chart>::Space, L>,
    Sl::Chart: SameChart<Sl2::Chart, L>,
    Sl::Var: Opposite<Sl2::Var, L>,
{
    type Found = True;
    type Rest = T;
}
impl<L, Sl, H, T: FindRemove<L, Sl>> FindBranch<False, L, Sl> for Cons<H, T> {
    type Found = <T as FindRemove<L, Sl>>::Found;
    type Rest = Cons<H, <T as FindRemove<L, Sl>>::Rest>;
}

/// Whether a labeled list contains label `L`.
pub trait ContainsLabel<L> {
    /// [`True`] if it does.
    type Out: Bool;
}
impl<L> ContainsLabel<L> for Nil {
    type Out = False;
}
impl<L: Label, L2: Label, Sl, T> ContainsLabel<L> for Cons<Ix<L2, Sl>, T>
where
    L2::Id: NatEq<L::Id>,
    T: ContainsLabel<L>,
    <L2::Id as NatEq<L::Id>>::Out: crate::utility::typelevel::Or<<T as ContainsLabel<L>>::Out>,
{
    type Out = <<L2::Id as NatEq<L::Id>>::Out as crate::utility::typelevel::Or<
        <T as ContainsLabel<L>>::Out,
    >>::Out;
}

/// Contracts every repeated label in a labeled list; `Out` is the free part and
/// `Dummies` the contracted labels.
pub trait Contract {
    /// The labeled slots that remain.
    type Out;
    /// The labels that were contracted (as labeled slots).
    type Dummies;
}
impl Contract for Nil {
    type Out = Nil;
    type Dummies = Nil;
}
impl<L, Sl, T> Contract for Cons<Ix<L, Sl>, T>
where
    T: FindRemove<L, Sl>,
    Self: ContractBranch<<T as FindRemove<L, Sl>>::Found>,
{
    type Out = <Self as ContractBranch<<T as FindRemove<L, Sl>>::Found>>::Out;
    type Dummies = <Self as ContractBranch<<T as FindRemove<L, Sl>>::Found>>::Dummies;
}

/// Helper of [`Contract`]: `B` tells whether the head's label appears again.
pub trait ContractBranch<B> {
    /// The labeled slots that remain.
    type Out;
    /// The labels that were contracted.
    type Dummies;
}
impl<L, Sl, T> ContractBranch<True> for Cons<Ix<L, Sl>, T>
where
    T: FindRemove<L, Sl>,
    <T as FindRemove<L, Sl>>::Rest: ContainsLabel<L> + Contract,
    <<T as FindRemove<L, Sl>>::Rest as ContainsLabel<L>>::Out: AppearsOnce<L>,
{
    type Out = <<T as FindRemove<L, Sl>>::Rest as Contract>::Out;
    type Dummies = Cons<Ix<L, Sl>, <<T as FindRemove<L, Sl>>::Rest as Contract>::Dummies>;
}
impl<L, Sl, T: Contract> ContractBranch<False> for Cons<Ix<L, Sl>, T> {
    type Out = Cons<Ix<L, Sl>, <T as Contract>::Out>;
    type Dummies = <T as Contract>::Dummies;
}

/// Every label of `Self` is absent from `U` (the labels already contracted in a
/// sub-expression). Otherwise the label would appear more than twice (ERR-06).
pub trait AllAbsent<U> {}
impl<U> AllAbsent<U> for Nil {}
impl<U, L, Sl, T> AllAbsent<U> for Cons<Ix<L, Sl>, T>
where
    U: ContainsLabel<L>,
    <U as ContainsLabel<L>>::Out: AppearsOnce<L>,
    T: AllAbsent<U>,
{
}

/// The slots of a labeled list, without the labels.
pub trait StripLabels {
    /// The slots.
    type Slots: Slots;
}
impl StripLabels for Nil {
    type Slots = Nil;
}
impl<L, Sl: SlotExtent + 'static, T: StripLabels> StripLabels for Cons<Ix<L, Sl>, T> {
    type Slots = Cons<Sl, T::Slots>;
}

/// Run-time information of a labeled list: the label number of each slot.
pub trait LabelIds {
    /// Number of slots.
    const RANK: usize;
    /// Writes the label numbers into `out`.
    fn ids(out: &mut [usize]);
}
impl LabelIds for Nil {
    const RANK: usize = 0;
    fn ids(_: &mut [usize]) {}
}
impl<L: Label, Sl, T: LabelIds> LabelIds for Cons<Ix<L, Sl>, T> {
    const RANK: usize = T::RANK + 1;
    fn ids(out: &mut [usize]) {
        out[0] = <L::Id as Nat>::VALUE;
        T::ids(&mut out[1..]);
    }
}

/// The free part of the concatenation of two labeled lists.
pub type Contracted<LA, LB> = <<LA as Concat<LB>>::Out as Contract>::Out;

// ---------------------------------------------------------------------------
// Contraction at the value level

/// Largest total rank (both factors together) that one contraction may involve.
pub const MAX_RANK: usize = 12;

/// Computes `Σ A_(…) B_(…)` over every repeated label.
///
/// `SA`/`SB` are the factors' slots, `SO` the result's slots; `LA`/`LB` give the
/// label numbers. The result has the free indices of A, then those of B, in order.
#[inline(always)]
pub fn contract<SA: Slots, SB: Slots, SO: Slots, LA: LabelIds, LB: LabelIds, N: Scalar>(
    a: &SA::Store<N>,
    b: &SB::Store<N>,
) -> SO::Store<N> {
    let ra = SA::RANK;
    let rb = SB::RANK;
    let r = ra + rb;
    assert!(r <= MAX_RANK, "contraction of total rank above {MAX_RANK}");
    let mut ids = [0usize; MAX_RANK];
    LA::ids(&mut ids[..ra]);
    LB::ids(&mut ids[ra..r]);
    let mut ext = [0usize; MAX_RANK];
    for (i, e) in ext.iter_mut().enumerate().take(r) {
        *e = if i < ra {
            SA::extent(i)
        } else {
            SB::extent(i - ra)
        };
    }
    // Role of each position: free output index k, or contracted pair j.
    let mut role_free = [usize::MAX; MAX_RANK];
    let mut role_pair = [usize::MAX; MAX_RANK];
    let mut pair_ext = [0usize; MAX_RANK];
    let (mut nfree, mut npair) = (0, 0);
    for i in 0..r {
        if role_pair[i] != usize::MAX {
            continue;
        }
        match (i + 1..r).find(|&j| ids[j] == ids[i]) {
            Some(j) => {
                role_pair[i] = npair;
                role_pair[j] = npair;
                pair_ext[npair] = ext[i];
                npair += 1;
            }
            None => {
                role_free[i] = nfree;
                nfree += 1;
            }
        }
    }
    // Strides of A and B (row-major).
    let mut stride = [0usize; MAX_RANK];
    let mut s = 1;
    for i in (0..ra).rev() {
        stride[i] = s;
        s *= ext[i];
    }
    s = 1;
    for i in (ra..r).rev() {
        stride[i] = s;
        s *= ext[i];
    }
    let mut out_ext = [0usize; MAX_RANK];
    for i in 0..r {
        if role_free[i] != usize::MAX {
            out_ext[role_free[i]] = ext[i];
        }
    }
    let pair_len: usize = pair_ext[..npair].iter().product();
    SO::from_fn(|o| {
        let mut out_idx = [0usize; MAX_RANK];
        let mut rem = o;
        for k in (0..nfree).rev() {
            out_idx[k] = rem % out_ext[k];
            rem /= out_ext[k];
        }
        let mut acc = N::zero();
        for c in 0..pair_len {
            let mut pair_idx = [0usize; MAX_RANK];
            let mut rem = c;
            for j in (0..npair).rev() {
                pair_idx[j] = rem % pair_ext[j];
                rem /= pair_ext[j];
            }
            let (mut fa, mut fb) = (0, 0);
            for i in 0..r {
                let v = if role_free[i] != usize::MAX {
                    out_idx[role_free[i]]
                } else {
                    pair_idx[role_pair[i]]
                };
                if i < ra {
                    fa += v * stride[i];
                } else {
                    fb += v * stride[i];
                }
            }
            acc = acc + SA::get(a, fa) * SB::get(b, fb);
        }
        acc
    })
}

// ---------------------------------------------------------------------------
// Coordinates and component dimensions

/// Position marker: the coordinate is the head of the list.
pub struct Here;
/// Position marker: the coordinate is in the tail, at position `I`.
pub struct There<I>(PhantomData<I>);

impl<X, C, T> CoordOf<X, C, Here> for Cons<X, T> {
    const POS: usize = 0;
}
impl<X, C, Y, T: CoordOf<X, C, I>, I> CoordOf<X, C, There<I>> for Cons<Y, T> {
    const POS: usize = T::POS + 1;
}

/// The element at position `I` of a type-level list.
pub trait At<I> {
    /// The element.
    type Out;
}
impl<H, T> At<Here> for Cons<H, T> {
    type Out = H;
}
impl<H, T: At<I>, I> At<There<I>> for Cons<H, T> {
    type Out = <T as At<I>>::Out;
}
