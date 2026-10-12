//! Compile-error messages (CMP-U03, DD-16, 2a §2a.4).
//!
//! Each kind of mistake has its own marker trait carrying the message. Engines decide
//! whether a combination is allowed by requiring these traits; when the requirement
//! fails, the compiler prints the message below as the headline.
//!
//! Every check trait also has a "decoy" impl for [`Decoy`]. With a single impl such as
//! `impl<D> SameDimension<D> for D`, the compiler would use it to infer the other
//! operand's type and report a plain type mismatch instead of our message. A second
//! candidate keeps inference from committing early, so the unmet bound is reported
//! with the headline below (DD-16).

use core::marker::PhantomData;

/// Never constructed; exists only to give each check trait a second impl.
pub struct Decoy<T>(PhantomData<T>, core::convert::Infallible);

/// ERR-01: both operands of `+` / `-` must have the same dimension.
#[diagnostic::on_unimplemented(
    message = "cannot add or subtract quantities with different dimensions",
    label = "this operand has a different physical dimension",
    note = "both operands must have the same physical dimension"
)]
pub trait SameDimension<D2> {}
impl<D> SameDimension<D> for D {}
impl<D> SameDimension<Decoy<D>> for D {}

/// The unit `Self` can be used for a quantity of dimension `D` (checked through
/// [`UnitDimensionIs`]).
pub trait UnitMatches<D> {}

/// ERR-02: implemented only for `True`, so a dimension mismatch reports this message.
#[diagnostic::on_unimplemented(
    message = "unit `{U}` does not have the dimension of this quantity",
    label = "wrong unit for this quantity",
    note = "use a unit of the same dimension, or convert the quantity first"
)]
pub trait UnitDimensionIs<U> {}
impl<U> UnitDimensionIs<U> for crate::utility::typelevel::True {}

/// ERR-03: tensors added or subtracted must have the same index structure.
#[diagnostic::on_unimplemented(
    message = "cannot add tensors with different index structure",
    label = "this tensor has different indices",
    note = "rank, index positions (upper/lower) and coordinate systems must all match"
)]
pub trait SameSlots<S2> {}
impl<S> SameSlots<S> for S {}
impl<S> SameSlots<Decoy<S>> for S {}

/// ERR-13: tensors added or subtracted must have the same symmetry.
#[diagnostic::on_unimplemented(
    message = "cannot add tensors with different symmetry",
    label = "this tensor has a different symmetry",
    note = "symmetrize or antisymmetrize one of them first"
)]
pub trait SameSymmetry<Sym2> {}
impl<S> SameSymmetry<S> for S {}
impl<S> SameSymmetry<Decoy<S>> for S {}

/// ERR-17: a tensor cannot be added to a chart-bound (non-tensor) quantity.
#[diagnostic::on_unimplemented(
    message = "cannot add a tensor and a chart-bound (non-tensor) quantity",
    label = "one operand is not a tensor",
    note = "use geometry::covariant_derivative to combine partial derivatives and Christoffel symbols"
)]
pub trait SameKind<K2> {}
impl<K> SameKind<K> for K {}
impl<K> SameKind<Decoy<K>> for K {}

/// ERR-04: one label per index. Implemented only for `True` (lengths equal).
#[diagnostic::on_unimplemented(
    message = "the number of index labels does not match the rank of the tensor",
    label = "give exactly one label per index"
)]
pub trait LabelCount<L> {}

/// ERR-05: a repeated label must be upper in one place and lower in the other.
#[diagnostic::on_unimplemented(
    message = "cannot contract index `{L}`: it is {Self} in both places",
    label = "repeated index with the same position",
    note = "a repeated index must be upper in one place and lower in the other; raise or lower one with a metric first"
)]
pub trait Opposite<V2, L> {}

/// ERR-06: a label may appear at most twice.
#[diagnostic::on_unimplemented(
    message = "index `{L}` appears more than twice",
    label = "repeated index",
    note = "use a different label for each pair of contracted indices"
)]
pub trait AppearsOnce<L> {}
impl<L> AppearsOnce<L> for crate::utility::typelevel::False {}

/// ERR-07: contracted indices must belong to the same space.
#[diagnostic::on_unimplemented(
    message = "index `{L}` belongs to different spaces",
    label = "indices from different spaces",
    note = "contracted indices must belong to the same space (e.g. both spacetime, or both field space)"
)]
pub trait SameSpace<Sp2, L> {}
impl<X, L> SameSpace<X, L> for X {}
impl<X, L> SameSpace<Decoy<X>, L> for X {}

/// ERR-08: contracted indices must be in the same chart.
#[diagnostic::on_unimplemented(
    message = "index `{L}` is in chart `{Self}` in one place and `{C2}` in the other",
    label = "indices in different charts",
    note = "transform one tensor to the other chart first"
)]
pub trait SameChart<C2, L> {}
impl<X, L> SameChart<X, L> for X {}
impl<X, L> SameChart<Decoy<X>, L> for X {}

/// ERR-11: the coordinate must belong to the chart of the index.
#[diagnostic::on_unimplemented(
    message = "coordinate `{X}` does not belong to chart `{C}`",
    label = "unknown coordinate for this index"
)]
pub trait CoordOf<X, C, I> {
    /// Position of the coordinate in the chart.
    const POS: usize;
}
