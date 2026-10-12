//! The numeric type of tensor components (CMP-U02, DD-09).

use core::fmt::Debug;
use core::ops::{Add, Div, Mul, Neg, Sub};

/// A number that tensor components can be made of.
///
/// `f64` is the ordinary choice. The automatic-differentiation types also implement it,
/// so that a metric or tensor field written generically over `N: Scalar` can be
/// evaluated for values and for derivatives from the same definition.
pub trait Scalar:
    Copy
    + Debug
    + PartialEq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + 'static
{
    /// Converts a constant.
    fn from_f64(x: f64) -> Self;
    /// The value part (for `f64`, the number itself).
    fn re(self) -> f64;
    /// Zero.
    fn zero() -> Self {
        Self::from_f64(0.0)
    }
    /// One.
    fn one() -> Self {
        Self::from_f64(1.0)
    }
    /// Square root.
    fn sqrt(self) -> Self;
    /// Exponential.
    fn exp(self) -> Self;
    /// Natural logarithm.
    fn ln(self) -> Self;
    /// Sine.
    fn sin(self) -> Self;
    /// Cosine.
    fn cos(self) -> Self;
    /// Tangent.
    fn tan(self) -> Self;
    /// Hyperbolic sine.
    fn sinh(self) -> Self;
    /// Hyperbolic cosine.
    fn cosh(self) -> Self;
    /// Hyperbolic tangent.
    fn tanh(self) -> Self;
    /// Arcsine.
    fn asin(self) -> Self;
    /// Arccosine.
    fn acos(self) -> Self;
    /// Arctangent.
    fn atan(self) -> Self;
    /// Integer power.
    fn powi(self, k: i32) -> Self;
    /// Real power.
    fn powf(self, p: f64) -> Self;
    /// Absolute value.
    fn abs(self) -> Self;
    /// Scales by a constant.
    fn scale(self, k: f64) -> Self {
        self * Self::from_f64(k)
    }
}

impl Scalar for f64 {
    fn from_f64(x: f64) -> Self {
        x
    }
    fn re(self) -> f64 {
        self
    }
    fn sqrt(self) -> Self {
        f64::sqrt(self)
    }
    fn exp(self) -> Self {
        f64::exp(self)
    }
    fn ln(self) -> Self {
        f64::ln(self)
    }
    fn sin(self) -> Self {
        f64::sin(self)
    }
    fn cos(self) -> Self {
        f64::cos(self)
    }
    fn tan(self) -> Self {
        f64::tan(self)
    }
    fn sinh(self) -> Self {
        f64::sinh(self)
    }
    fn cosh(self) -> Self {
        f64::cosh(self)
    }
    fn tanh(self) -> Self {
        f64::tanh(self)
    }
    fn asin(self) -> Self {
        f64::asin(self)
    }
    fn acos(self) -> Self {
        f64::acos(self)
    }
    fn atan(self) -> Self {
        f64::atan(self)
    }
    fn powi(self, k: i32) -> Self {
        f64::powi(self, k)
    }
    fn powf(self, p: f64) -> Self {
        f64::powf(self, p)
    }
    fn abs(self) -> Self {
        f64::abs(self)
    }
    fn scale(self, k: f64) -> Self {
        self * k
    }
}
