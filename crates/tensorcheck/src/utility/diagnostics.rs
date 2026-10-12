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
