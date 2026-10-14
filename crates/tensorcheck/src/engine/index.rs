//! Index slots, spaces, charts, labels and contraction (CMP-E02, DD-04, DD-06).
//!
//! Type level: which indices can be contracted and what slots remain. Value level:
//! the contraction sums, computed directly without forming outer products.

use core::marker::PhantomData;

use crate::access::component::{Extent, SlotExtent, Slots};
pub use crate::contracts::Ix;
pub use crate::contracts::{
    Chart, ChartBound, CoordinateBasis, Down, Label, NonCoordinateBasis, Point, Slot, Space,
    Tensorial, Up, lower, slot, upper,
};
use crate::utility::diagnostics::{
    AppearsOnce, CoordOf, LabelCount, Opposite, SameChart, SameSpace,
};
use crate::utility::scalar::Scalar;
use crate::utility::typelevel::{Bool, Concat, Cons, False, Nat, NatEq, Nil, True};

// ---------------------------------------------------------------------------
// Tensoriality

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

impl<L> Opposite<lower, L> for upper {}
impl<L> Opposite<upper, L> for lower {}

// ---------------------------------------------------------------------------
// Labels

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

/// The label number of each slot of a labeled list.
pub trait LabelIds {
    /// Number of slots.
    const RANK: usize;
    /// Label numbers, padded with zeros, as a constant.
    const IDS: [usize; MAX_RANK];
    /// Writes the label numbers into `out`.
    fn ids(out: &mut [usize]) {
        out.copy_from_slice(&Self::IDS[..out.len()]);
    }
}
impl LabelIds for Nil {
    const RANK: usize = 0;
    const IDS: [usize; MAX_RANK] = [0; MAX_RANK];
}
impl<L: Label, Sl, T: LabelIds> LabelIds for Cons<Ix<L, Sl>, T> {
    const RANK: usize = T::RANK + 1;
    const IDS: [usize; MAX_RANK] = crate::utility::extent::prepend(<L::Id as Nat>::VALUE, T::IDS);
}

/// The free part of the concatenation of two labeled lists.
pub type Contracted<LA, LB> = <<LA as Concat<LB>>::Out as Contract>::Out;

// ---------------------------------------------------------------------------
// Contraction at the value level

/// Largest total rank (both factors together) that one contraction may involve.
pub use crate::utility::extent::MAX_RANK;

/// Largest table in a constant contraction plan (outputs or contracted index
/// combinations). Larger contractions use the run-time algorithm.
pub const PLAN_TABLE: usize = 1024;

/// A contraction computed once per combination of types, at compile time:
/// `out[o] = Σ_c A[a_off[o] + pa[c]] · B[b_off[o] + pb[c]]` on flat slices.
pub struct Plan {
    /// Whether the tables were large enough (otherwise use the run-time algorithm).
    pub fits: bool,
    /// Number of output components.
    pub out_len: usize,
    /// Number of combinations of contracted indices.
    pub pair_len: usize,
    /// Offset in A of output component `o` (free indices only).
    pub a_off: [u32; PLAN_TABLE],
    /// Offset in B of output component `o`.
    pub b_off: [u32; PLAN_TABLE],
    /// Offset in A of contracted combination `c`.
    pub pa: [u32; PLAN_TABLE],
    /// Offset in B of contracted combination `c`.
    pub pb: [u32; PLAN_TABLE],
    /// If `pa[c] = c · step_a` and `pb[c] = c · step_b` for every `c`, the steps
    /// (the inner sum is then a strided dot product); otherwise 0.
    pub step_a: usize,
    /// See `step_a`.
    pub step_b: usize,
    /// Whether the contracted offsets are an arithmetic progression.
    pub affine: bool,
}

impl Plan {
    /// Builds the plan from the factors' ranks, extents and label numbers.
    pub const fn make(
        ra: usize,
        rb: usize,
        ext_a: [usize; MAX_RANK],
        ext_b: [usize; MAX_RANK],
        ids_a: [usize; MAX_RANK],
        ids_b: [usize; MAX_RANK],
    ) -> Plan {
        let r = ra + rb;
        let mut ids = [0usize; 2 * MAX_RANK];
        let mut ext = [0usize; 2 * MAX_RANK];
        let mut i = 0;
        while i < ra {
            ids[i] = ids_a[i];
            ext[i] = ext_a[i];
            i += 1;
        }
        i = 0;
        while i < rb {
            ids[ra + i] = ids_b[i];
            ext[ra + i] = ext_b[i];
            i += 1;
        }
        // Strides within A and within B (row-major).
        let mut stride = [0usize; 2 * MAX_RANK];
        let mut s = 1;
        i = ra;
        while i > 0 {
            i -= 1;
            stride[i] = s;
            s *= ext[i];
        }
        s = 1;
        i = r;
        while i > ra {
            i -= 1;
            stride[i] = s;
            s *= ext[i];
        }
        // Roles: free output index number, or contracted pair number.
        const NONE: usize = usize::MAX;
        let mut free = [NONE; 2 * MAX_RANK];
        let mut pair = [NONE; 2 * MAX_RANK];
        let mut pair_ext = [0usize; MAX_RANK];
        let mut out_ext = [0usize; 2 * MAX_RANK];
        let (mut nfree, mut npair) = (0, 0);
        i = 0;
        while i < r {
            if pair[i] == NONE {
                let mut j = i + 1;
                let mut found = NONE;
                while j < r {
                    if ids[j] == ids[i] && found == NONE {
                        found = j;
                    }
                    j += 1;
                }
                if found == NONE {
                    free[i] = nfree;
                    out_ext[nfree] = ext[i];
                    nfree += 1;
                } else {
                    pair[i] = npair;
                    pair[found] = npair;
                    pair_ext[npair] = ext[i];
                    npair += 1;
                }
            }
            i += 1;
        }
        let mut out_len = 1;
        i = 0;
        while i < nfree {
            out_len *= out_ext[i];
            i += 1;
        }
        let mut pair_len = 1;
        i = 0;
        while i < npair {
            pair_len *= pair_ext[i];
            i += 1;
        }
        let mut plan = Plan {
            fits: out_len <= PLAN_TABLE && pair_len <= PLAN_TABLE,
            out_len,
            pair_len,
            a_off: [0; PLAN_TABLE],
            b_off: [0; PLAN_TABLE],
            pa: [0; PLAN_TABLE],
            pb: [0; PLAN_TABLE],
            step_a: 0,
            step_b: 0,
            affine: false,
        };
        if !plan.fits {
            return plan;
        }
        // Offsets of every output component.
        let mut o = 0;
        while o < out_len {
            let mut idx = [0usize; 2 * MAX_RANK];
            let mut rem = o;
            let mut k = nfree;
            while k > 0 {
                k -= 1;
                idx[k] = rem % out_ext[k];
                rem /= out_ext[k];
            }
            let (mut fa, mut fb) = (0, 0);
            i = 0;
            while i < r {
                if free[i] != NONE {
                    if i < ra {
                        fa += idx[free[i]] * stride[i];
                    } else {
                        fb += idx[free[i]] * stride[i];
                    }
                }
                i += 1;
            }
            plan.a_off[o] = fa as u32;
            plan.b_off[o] = fb as u32;
            o += 1;
        }
        // Offsets of every combination of contracted indices.
        let mut c = 0;
        while c < pair_len {
            let mut idx = [0usize; MAX_RANK];
            let mut rem = c;
            let mut k = npair;
            while k > 0 {
                k -= 1;
                idx[k] = rem % pair_ext[k];
                rem /= pair_ext[k];
            }
            let (mut da, mut db) = (0, 0);
            i = 0;
            while i < r {
                if pair[i] != NONE {
                    if i < ra {
                        da += idx[pair[i]] * stride[i];
                    } else {
                        db += idx[pair[i]] * stride[i];
                    }
                }
                i += 1;
            }
            plan.pa[c] = da as u32;
            plan.pb[c] = db as u32;
            c += 1;
        }
        // Arithmetic progression? (true for contracting trailing/leading indices).
        let (sa, sb) = if pair_len > 1 {
            (plan.pa[1] as usize, plan.pb[1] as usize)
        } else {
            (0, 0)
        };
        let mut affine = true;
        c = 0;
        while c < pair_len {
            if plan.pa[c] as usize != c * sa || plan.pb[c] as usize != c * sb {
                affine = false;
            }
            c += 1;
        }
        plan.affine = affine;
        plan.step_a = sa;
        plan.step_b = sb;
        plan
    }
}

/// The constant plan for contracting `SA` (labels `LA`) with `SB` (labels `LB`).
///
/// The scalar facts are separate constants so that loop bounds and strides are
/// compile-time constants the optimizer can unroll and vectorize.
pub trait ContractPlan {
    /// The plan, evaluated at compile time for each combination of types.
    const PLAN: Plan;
    /// `PLAN.fits`.
    const FITS: bool;
    /// `PLAN.out_len`.
    const OUT_LEN: usize;
    /// `PLAN.pair_len`.
    const PAIR_LEN: usize;
    /// Whether contracted offsets are `c · step`.
    const AFFINE: bool;
    /// Contracted step in A.
    const STEP_A: usize;
    /// Contracted step in B.
    const STEP_B: usize;
    /// Whether output offsets are `o · step` in both factors.
    const OUT_AFFINE: bool;
    /// Output step in A.
    const OUT_STEP_A: usize;
    /// Output step in B.
    const OUT_STEP_B: usize;
}
impl<SA: Slots, SB: Slots, LA: LabelIds, LB: LabelIds> ContractPlan for (SA, SB, LA, LB) {
    const PLAN: Plan = Plan::make(
        SA::RANK,
        SB::RANK,
        SA::EXTENTS,
        SB::EXTENTS,
        LA::IDS,
        LB::IDS,
    );
    const FITS: bool = Self::PLAN.fits;
    const OUT_LEN: usize = Self::PLAN.out_len;
    const PAIR_LEN: usize = Self::PLAN.pair_len;
    const AFFINE: bool = Self::PLAN.affine;
    const STEP_A: usize = Self::PLAN.step_a;
    const STEP_B: usize = Self::PLAN.step_b;
    const OUT_AFFINE: bool = out_affine(&Self::PLAN).0;
    const OUT_STEP_A: usize = out_affine(&Self::PLAN).1;
    const OUT_STEP_B: usize = out_affine(&Self::PLAN).2;
}

/// Whether `a_off[o] = o · sa` and `b_off[o] = o · sb` for every output `o`.
const fn out_affine(p: &Plan) -> (bool, usize, usize) {
    if !p.fits {
        return (false, 0, 0);
    }
    let (sa, sb) = if p.out_len > 1 {
        (p.a_off[1] as usize, p.b_off[1] as usize)
    } else {
        (0, 0)
    };
    let mut o = 0;
    while o < p.out_len {
        if p.a_off[o] as usize != o * sa || p.b_off[o] as usize != o * sb {
            return (false, 0, 0);
        }
        o += 1;
    }
    (true, sa, sb)
}

/// Computes `Σ A_(…) B_(…)` over every repeated label.
///
/// `SA`/`SB` are the factors' slots, `SO` the result's slots; `LA`/`LB` give the
/// label numbers. The result has the free indices of A, then those of B, in order.
#[inline(always)]
pub fn contract<SA: Slots, SB: Slots, SO: Slots, LA: LabelIds, LB: LabelIds, N: Scalar>(
    a: &SA::Store<N>,
    b: &SB::Store<N>,
) -> SO::Store<N> {
    type P<SA, SB, LA, LB> = (SA, SB, LA, LB);
    if <P<SA, SB, LA, LB> as ContractPlan>::FITS {
        let (fa, fb) = (SA::flat(a), SB::flat(b));
        let len = <P<SA, SB, LA, LB> as ContractPlan>::PAIR_LEN;
        let (sa, sb) = (
            <P<SA, SB, LA, LB> as ContractPlan>::STEP_A,
            <P<SA, SB, LA, LB> as ContractPlan>::STEP_B,
        );
        if <P<SA, SB, LA, LB> as ContractPlan>::AFFINE
            && <P<SA, SB, LA, LB> as ContractPlan>::OUT_AFFINE
        {
            // Both loops are affine with constant bounds and steps: the fast path.
            let (oa, ob) = (
                <P<SA, SB, LA, LB> as ContractPlan>::OUT_STEP_A,
                <P<SA, SB, LA, LB> as ContractPlan>::OUT_STEP_B,
            );
            return SO::from_fn(|o| {
                let (a0, b0) = (o * oa, o * ob);
                let mut acc = N::zero();
                for c in 0..len {
                    acc = acc + fa[a0 + c * sa] * fb[b0 + c * sb];
                }
                acc
            });
        }
        let plan = &<P<SA, SB, LA, LB> as ContractPlan>::PLAN;
        if <P<SA, SB, LA, LB> as ContractPlan>::AFFINE {
            return SO::from_fn(|o| {
                let (a0, b0) = (plan.a_off[o] as usize, plan.b_off[o] as usize);
                let mut acc = N::zero();
                for c in 0..len {
                    acc = acc + fa[a0 + c * sa] * fb[b0 + c * sb];
                }
                acc
            });
        }
        return SO::from_fn(|o| {
            let (a0, b0) = (plan.a_off[o] as usize, plan.b_off[o] as usize);
            let mut acc = N::zero();
            for c in 0..len {
                acc = acc + fa[a0 + plan.pa[c] as usize] * fb[b0 + plan.pb[c] as usize];
            }
            acc
        });
    }
    contract_runtime::<SA, SB, SO, LA, LB, N>(a, b)
}

/// The run-time contraction algorithm, for contractions too large for a plan table.
fn contract_runtime<SA: Slots, SB: Slots, SO: Slots, LA: LabelIds, LB: LabelIds, N: Scalar>(
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

// ---------------------------------------------------------------------------
// Points and chart transforms

/// A change of chart on one space (EXT-05).
pub trait ChartTransform: 'static {
    /// Source chart.
    type From: Chart;
    /// Target chart (on the same space).
    type To: Chart<Space = <Self::From as Chart>::Space>;
    /// The coordinates of `p` in the target chart.
    fn map<N: Scalar>(&self, p: &Point<Self::From, N>) -> Point<Self::To, N>;
    /// The Jacobian `∂x'^a/∂x^b` at `p`, row-major (`a` is the row), in base units.
    fn jacobian<N: Scalar>(&self, p: &Point<Self::From, N>) -> Vec<N>;
    /// The inverse Jacobian `∂x^b/∂x'^a` at `p`, row-major (`b` is the row), if known
    /// in closed form.
    fn inverse_jacobian<N: Scalar>(&self, _p: &Point<Self::From, N>) -> Option<Vec<N>> {
        None
    }
}

/// The inverse Jacobian `∂x^b/∂x'^a` at `p`: the closed form if given, otherwise
/// the numerical inverse (U04).
pub fn inverse_jacobian_at<X: ChartTransform, N: Scalar>(
    xf: &X,
    p: &Point<X::From, N>,
) -> Result<Vec<N>, crate::utility::linear::LinearError> {
    if let Some(j) = xf.inverse_jacobian(p) {
        return Ok(j);
    }
    let n = <<X::From as Chart>::Space as Extent>::N;
    let j = xf.jacobian(p);
    let mut inv = vec![N::zero(); n * n];
    crate::utility::linear::inverse(&j, n, &mut inv)?;
    Ok(inv)
}

/// Applies the matrix `mat` (row-major, `n × n`) to slot number `slot` of a flat
/// component array: `out[.., i, ..] = Σ_k mat[i][k] · in[.., k, ..]`.
pub fn apply_at_slot<N: Scalar>(data: &[N], extents: &[usize], slot: usize, mat: &[N]) -> Vec<N> {
    let n = extents[slot];
    let inner: usize = extents[slot + 1..].iter().product();
    let outer: usize = extents[..slot].iter().product();
    let mut out = vec![N::zero(); data.len()];
    for o in 0..outer {
        for i in 0..n {
            for r in 0..inner {
                let mut acc = N::zero();
                for k in 0..n {
                    acc = acc + mat[i * n + k] * data[(o * n + k) * inner + r];
                }
                out[(o * n + i) * inner + r] = acc;
            }
        }
    }
    out
}

/// Type-level access to slot number `P` of a slot list, and replacement of it.
pub trait SlotAt<P> {
    /// The slot at `P`.
    type Slot: Slot;
    /// The list with the slot at `P` flipped (upper ↔ lower).
    type Flipped;
    /// The position as a number.
    const POS: usize;
}
impl<H: Slot, T> SlotAt<slot::_0> for Cons<H, T> {
    type Slot = H;
    type Flipped = Cons<H::Flip, T>;
    const POS: usize = 0;
}
impl<H, T: SlotAt<slot::_0>> SlotAt<slot::_1> for Cons<H, T> {
    type Slot = <T as SlotAt<slot::_0>>::Slot;
    type Flipped = Cons<H, <T as SlotAt<slot::_0>>::Flipped>;
    const POS: usize = 1;
}
impl<H, T: SlotAt<slot::_1>> SlotAt<slot::_2> for Cons<H, T> {
    type Slot = <T as SlotAt<slot::_1>>::Slot;
    type Flipped = Cons<H, <T as SlotAt<slot::_1>>::Flipped>;
    const POS: usize = 2;
}
impl<H, T: SlotAt<slot::_2>> SlotAt<slot::_3> for Cons<H, T> {
    type Slot = <T as SlotAt<slot::_2>>::Slot;
    type Flipped = Cons<H, <T as SlotAt<slot::_2>>::Flipped>;
    const POS: usize = 3;
}

/// Replaces the chart of every slot by `To` (all slots must be in `From`, ERR-10).
pub trait MapChart<From, To, X> {
    /// The slots in the target chart.
    type Out;
}
impl<From, To, X> MapChart<From, To, X> for Nil {
    type Out = Nil;
}
impl<From, To, X, C, T: MapChart<From, To, X>> MapChart<From, To, X> for Cons<Up<C>, T>
where
    C: crate::utility::diagnostics::TransformFrom<From, X>,
{
    type Out = Cons<Up<To>, <T as MapChart<From, To, X>>::Out>;
}
impl<From, To, X, C, T: MapChart<From, To, X>> MapChart<From, To, X> for Cons<Down<C>, T>
where
    C: crate::utility::diagnostics::TransformFrom<From, X>,
{
    type Out = Cons<Down<To>, <T as MapChart<From, To, X>>::Out>;
}

/// Run-time variance of each slot (`true` = upper).
pub trait SlotVariances {
    /// Writes the variances into `out`.
    fn write(out: &mut [bool]);
}
impl SlotVariances for Nil {
    fn write(_: &mut [bool]) {}
}
impl<H: Slot, T: SlotVariances> SlotVariances for Cons<H, T> {
    fn write(out: &mut [bool]) {
        out[0] = H::UP;
        T::write(&mut out[1..]);
    }
}

impl<P> crate::utility::diagnostics::VarianceIs<lower, P> for lower {}
impl<P> crate::utility::diagnostics::VarianceIs<upper, P> for upper {}
