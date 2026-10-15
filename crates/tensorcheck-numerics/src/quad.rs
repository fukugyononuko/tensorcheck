//! Definite integrals that keep dimensions (API-14, API-46 to API-48).
//!
//! ```
//! use tensorcheck::Quantity;
//! use tensorcheck::units::si::{Energy, Force, J, Length, N, m};
//! use tensorcheck_numerics::quad::{GaussKronrod, Tolerance, integrate};
//!
//! // Work of a spring force F(x) = −k x from 0 to 0.1 m with k = 200 N/m.
//! let force = |x: Quantity<Length>| -> Quantity<Force> { Quantity::new(-200.0 * x.value_in(m), N) };
//! let w = integrate(force, (0.0 * m)..(0.1 * m), GaussKronrod, Tolerance { rtol: 1e-12, atol: 1e-15 * J }).unwrap();
//! let work: Quantity<Energy> = w.value;
//! assert!((work.value_in(J) + 1.0).abs() < 1e-12);
//! ```

use core::fmt;
use core::ops::Range;

use tensorcheck::utility::diagnostics::ToleranceFor;
use tensorcheck::{Quantity, TimesDim, VectorSpace};

use crate::manager::solver;

/// Composite trapezoid rule with `n` intervals (error estimate from `n` vs `2n`).
#[derive(Clone, Copy, Debug)]
pub struct Trapezoid {
    /// Number of intervals, at least 1 (else [`QuadError::InvalidCount`]).
    pub n: usize,
}
/// Composite Simpson rule with `n` intervals (rounded up to even).
#[derive(Clone, Copy, Debug)]
pub struct Simpson {
    /// Number of intervals, at least 1 (else [`QuadError::InvalidCount`]).
    pub n: usize,
}
/// Romberg integration (trapezoid with Richardson extrapolation, up to 20 levels).
#[derive(Clone, Copy, Debug)]
pub struct Romberg;
/// `n`-point Gauss–Legendre rule (error estimate from `n` vs `n + 1` points).
#[derive(Clone, Copy, Debug)]
pub struct GaussLegendre {
    /// Number of points, at least 1 (else [`QuadError::InvalidCount`]).
    pub n: usize,
}
/// Adaptive 7-point Gauss / 15-point Kronrod rule (QUADPACK `qag` style).
#[derive(Clone, Copy, Debug)]
pub struct GaussKronrod;

/// A quadrature method.
pub trait Method: Copy {
    #[doc(hidden)]
    fn kind(self) -> solver::QuadKind;
}
impl Method for Trapezoid {
    fn kind(self) -> solver::QuadKind {
        solver::QuadKind::Trapezoid(self.n)
    }
}
impl Method for Simpson {
    fn kind(self) -> solver::QuadKind {
        solver::QuadKind::Simpson(self.n)
    }
}
impl Method for Romberg {
    fn kind(self) -> solver::QuadKind {
        solver::QuadKind::Romberg
    }
}
impl Method for GaussLegendre {
    fn kind(self) -> solver::QuadKind {
        solver::QuadKind::GaussLegendre(self.n)
    }
}
impl Method for GaussKronrod {
    fn kind(self) -> solver::QuadKind {
        solver::QuadKind::GaussKronrod
    }
}

/// Error tolerance: relative `rtol` (a plain number) and absolute `atol` with the
/// type of the result (API-47), or of the state for ODEs (API-43).
///
/// `rtol` and every component of `atol` must be finite and nonnegative, and each
/// component needs `atol > 0` or `rtol > 0`. Anything else is rejected before the
/// integrand or right-hand side is evaluated ([`QuadError::InvalidTolerance`],
/// [`SolveError::InvalidTolerance`](crate::ode::SolveError::InvalidTolerance)).
#[derive(Clone, Copy, Debug)]
pub struct Tolerance<A> {
    /// Relative tolerance.
    pub rtol: f64,
    /// Absolute tolerance, one per component, with the result's dimension.
    pub atol: A,
}

/// A computed integral (API-48).
#[derive(Clone, Copy, Debug)]
pub struct Integral<V> {
    /// The integral.
    pub value: V,
    /// The estimated absolute error, per component, in base units (largest).
    pub error_estimate: f64,
    /// Number of integrand evaluations.
    pub evaluations: usize,
}

/// Why an integral failed (API-48).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum QuadError {
    /// The estimated error stayed above the tolerance.
    ToleranceNotMet {
        /// The scaled error estimate (above 1 means above tolerance).
        scaled_error: f64,
    },
    /// The integrand returned NaN or an infinity.
    NonFinite,
    /// The tolerance is outside its domain (see [`Tolerance`]); nothing was evaluated.
    InvalidTolerance,
    /// The method's count of intervals or points is unsupported (zero, or so large
    /// that the refined count overflows); nothing was evaluated.
    InvalidCount {
        /// The count given.
        n: usize,
    },
    /// An end of the range is NaN or infinite; nothing was evaluated.
    InvalidRange,
}
impl fmt::Display for QuadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            QuadError::ToleranceNotMet { scaled_error } => {
                write!(
                    f,
                    "tolerance not met (scaled error estimate {scaled_error:.3e})"
                )
            }
            QuadError::NonFinite => f.write_str("the integrand returned a non-finite value"),
            QuadError::InvalidTolerance => f.write_str(
                "invalid tolerance (rtol and atol must be finite and nonnegative, not both zero)",
            ),
            QuadError::InvalidCount { n } => {
                write!(f, "unsupported number of intervals or points: {n}")
            }
            QuadError::InvalidRange => f.write_str("the integration range is not finite"),
        }
    }
}
impl std::error::Error for QuadError {}

/// `∫_a^b f(x) dx` (API-46). The result has the dimension of `f` times that of `x`;
/// `atol` must have that type (ERR-15). The tolerance, the method's count and the
/// range are checked before `f` is called.
pub fn integrate<X, Y, F, A, M>(
    f: F,
    range: Range<Quantity<X>>,
    method: M,
    tol: Tolerance<A>,
) -> Result<Integral<<Y as TimesDim<X>>::Out>, QuadError>
where
    F: FnMut(Quantity<X>) -> Y,
    Y: VectorSpace + TimesDim<X>,
    A: VectorSpace + ToleranceFor<<Y as TimesDim<X>>::Out>,
    M: Method,
{
    solver::integrate(f, range, method.kind(), tol.rtol, &tol.atol)
}
