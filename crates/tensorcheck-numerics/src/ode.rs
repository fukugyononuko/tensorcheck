//! Ordinary differential equations that keep dimensions (API-13, API-40 to API-45).
//!
//! One entry point, [`solve_ivp`], for systems `dy/dt = f(t, y)` of any size. Write
//! higher-order equations as first-order systems: for `x'' = a(x)`, use the state
//! `(x, v)` and return `(v, a(x))`.
//!
//! ```
//! use tensorcheck::Quantity;
//! use tensorcheck::units::si::{Frequency, Time, one, s};
//! use tensorcheck_numerics::ode::{DormandPrince54, Tolerance, solve_ivp};
//!
//! // dy/dt = −y / τ with τ = 2 s.
//! let tau = 2.0 * s;
//! let rhs = |_t: Quantity<Time>, y: &Quantity<tensorcheck::units::si::Dimensionless>| -> Quantity<Frequency> { -*y / tau };
//! let sol = solve_ivp(rhs, (0.0 * s)..(1.0 * s), 1.0 * one, DormandPrince54, Tolerance { rtol: 1e-10, atol: 1e-12 * one }).unwrap();
//! assert!((sol.y.value_in(one) - (-0.5f64).exp()).abs() < 1e-9);
//! ```

use core::fmt;
use core::ops::Range;

use tensorcheck::utility::diagnostics::{RhsReturns, ToleranceFor};
use tensorcheck::{PerDim, Quantity, VectorSpace};

use crate::manager::solver;

/// The classical fourth-order Runge–Kutta method with a fixed step `h` (in base
/// units of the independent variable). No error control.
#[derive(Clone, Copy, Debug)]
pub struct Rk4 {
    /// Step size in base units: a finite positive magnitude (else
    /// [`SolveError::InvalidStep`]). The direction comes from the interval.
    pub h: f64,
}
/// Dormand–Prince 5(4): explicit, adaptive step size, for non-stiff problems.
#[derive(Clone, Copy, Debug)]
pub struct DormandPrince54;

/// RODAS (Hairer & Wanner): linearly implicit Rosenbrock method of order 4(3),
/// L-stable and stiffly accurate, adaptive step size, for stiff problems.
#[derive(Clone, Copy, Debug)]
pub struct Rosenbrock {
    /// How the Jacobian `∂f/∂y` is obtained.
    pub jacobian: Jacobian,
}

/// How the Jacobian is obtained for [`Rosenbrock`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Jacobian {
    /// Central differences at every step (DD-15).
    Numerical,
}

/// An ODE method.
pub trait Method: Copy {
    #[doc(hidden)]
    fn kind(self) -> solver::OdeKind;
    #[doc(hidden)]
    fn step_limit(self) -> usize {
        MAX_STEPS
    }

    /// The same method with at most `n` steps instead of [`MAX_STEPS`] (API-57).
    /// Adaptive methods count rejected steps too. Reaching the limit returns
    /// [`SolveError::MaxSteps`] with the time reached; with `n = 0` that is the start
    /// time (unless the interval is empty).
    ///
    /// ```
    /// use tensorcheck::Quantity;
    /// use tensorcheck::units::si::{Dimensionless, Frequency, Time, one, s};
    /// use tensorcheck_numerics::ode::{Method, Rk4, SolveError, Tolerance, solve_ivp};
    ///
    /// let rhs = |_t: Quantity<Time>, _y: &Quantity<Dimensionless>| -> Quantity<Frequency> { 1.0 * (one / s) };
    /// let tol = Tolerance { rtol: 1e-6, atol: 1e-9 * one };
    /// // Four steps of 0.25 s cover 1 s; a limit of 2 stops at 0.5 s.
    /// let r = solve_ivp(rhs, (0.0 * s)..(1.0 * s), 0.0 * one, Rk4 { h: 0.25 }.max_steps(2), tol);
    /// assert!(matches!(r, Err(SolveError::MaxSteps { t }) if t.value_in(s) == 0.5));
    /// ```
    fn max_steps(self, n: usize) -> StepLimit<Self> {
        StepLimit {
            method: self,
            max_steps: n,
        }
    }
}

/// A method with its own step limit, made by [`Method::max_steps`].
#[derive(Clone, Copy, Debug)]
pub struct StepLimit<M> {
    /// The underlying method.
    pub method: M,
    /// The largest number of steps.
    pub max_steps: usize,
}
impl<M: Method> Method for StepLimit<M> {
    fn kind(self) -> solver::OdeKind {
        self.method.kind()
    }
    fn step_limit(self) -> usize {
        self.max_steps
    }
}
impl Method for Rk4 {
    fn kind(self) -> solver::OdeKind {
        solver::OdeKind::Rk4(self.h)
    }
}
impl Method for DormandPrince54 {
    fn kind(self) -> solver::OdeKind {
        solver::OdeKind::Dopri5
    }
}
impl Method for Rosenbrock {
    fn kind(self) -> solver::OdeKind {
        solver::OdeKind::Rodas
    }
}

/// Error tolerance: relative `rtol` and absolute `atol` with the state's type
/// (API-43).
pub use crate::quad::Tolerance;

/// Counters of a solve.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Stats {
    /// Accepted steps.
    pub steps: usize,
    /// Rejected steps.
    pub rejected: usize,
    /// Right-hand-side evaluations.
    pub evaluations: usize,
}

/// The state at the end of the interval (API-44).
pub struct Solution<T, Y> {
    /// The final time.
    pub t: Quantity<T>,
    /// The final state.
    pub y: Y,
    /// Counters.
    pub stats: Stats,
}

/// Why a solve failed, with the time it happened (API-44).
pub enum SolveError<T> {
    /// The step limit was reached.
    MaxSteps {
        /// Time reached.
        t: Quantity<T>,
    },
    /// The step size became too small for the arithmetic (`|h| < 16 ε |t|`).
    StepTooSmall {
        /// Time reached.
        t: Quantity<T>,
        /// The step size in base units.
        h: f64,
    },
    /// The state or its derivative became NaN or infinite.
    NonFinite {
        /// Time reached.
        t: Quantity<T>,
    },
    /// A linear system in an implicit method was singular.
    Singular {
        /// Time reached.
        t: Quantity<T>,
    },
    /// The tolerance is outside its domain (see [`Tolerance`]); nothing was evaluated.
    InvalidTolerance,
    /// The fixed step of [`Rk4`] is not finite and positive; nothing was evaluated.
    InvalidStep {
        /// The step given, in base units.
        h: f64,
    },
    /// An end of the interval is NaN or infinite; nothing was evaluated.
    InvalidSpan,
}
impl<T> fmt::Display for SolveError<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (what, t) = match self {
            SolveError::MaxSteps { t } => ("step limit reached", t.base_value()),
            SolveError::StepTooSmall { t, .. } => ("step size too small", t.base_value()),
            SolveError::NonFinite { t } => ("non-finite value", t.base_value()),
            SolveError::Singular { t } => ("singular linear system", t.base_value()),
            SolveError::InvalidTolerance => {
                return f.write_str(
                    "invalid tolerance (rtol and atol must be finite and nonnegative, not both zero)",
                );
            }
            SolveError::InvalidStep { h } => {
                return write!(
                    f,
                    "invalid fixed step h = {h:e} (must be finite and positive)"
                );
            }
            SolveError::InvalidSpan => return f.write_str("the time interval is not finite"),
        };
        write!(f, "{what} at t = {t:e} (base units)")
    }
}
impl<T> Clone for SolveError<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for SolveError<T> {}
impl<T> PartialEq for SolveError<T> {
    fn eq(&self, o: &Self) -> bool {
        match (self, o) {
            (SolveError::MaxSteps { t: a }, SolveError::MaxSteps { t: b })
            | (SolveError::NonFinite { t: a }, SolveError::NonFinite { t: b })
            | (SolveError::Singular { t: a }, SolveError::Singular { t: b }) => {
                a.base_value() == b.base_value()
            }
            (SolveError::StepTooSmall { t: a, h: x }, SolveError::StepTooSmall { t: b, h: y }) => {
                a.base_value() == b.base_value() && x == y
            }
            (SolveError::InvalidTolerance, SolveError::InvalidTolerance)
            | (SolveError::InvalidSpan, SolveError::InvalidSpan) => true,
            (SolveError::InvalidStep { h: x }, SolveError::InvalidStep { h: y }) => x == y,
            _ => false,
        }
    }
}
impl<T> fmt::Debug for SolveError<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}
impl<T> std::error::Error for SolveError<T> {}

impl<T> SolveError<T> {
    /// The time at which the solve failed; `None` when the input was rejected before
    /// the solve started.
    pub fn time(&self) -> Option<Quantity<T>> {
        match *self {
            SolveError::MaxSteps { t }
            | SolveError::NonFinite { t }
            | SolveError::Singular { t } => Some(t),
            SolveError::StepTooSmall { t, .. } => Some(t),
            SolveError::InvalidTolerance
            | SolveError::InvalidStep { .. }
            | SolveError::InvalidSpan => None,
        }
    }
}

impl<T, Y: Clone> Clone for Solution<T, Y> {
    fn clone(&self) -> Self {
        Solution {
            t: self.t,
            y: self.y.clone(),
            stats: self.stats,
        }
    }
}
impl<T, Y: Copy> Copy for Solution<T, Y> {}
impl<T, Y: fmt::Debug> fmt::Debug for Solution<T, Y> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Solution")
            .field("t (base units)", &self.t.base_value())
            .field("y", &self.y)
            .field("stats", &self.stats)
            .finish()
    }
}

/// Default largest number of steps (DD-13); change it with [`Method::max_steps`].
pub const MAX_STEPS: usize = 100_000;

/// Solves `dy/dt = f(t, y)` from `span.start` to `span.end` (which may be earlier)
/// (API-40). `f` must return the state type divided by the time dimension (ERR-14)
/// and `atol` must have the state type (ERR-15).
///
/// Before `f` is called, the interval must be finite ([`SolveError::InvalidSpan`]),
/// the tolerance in its domain ([`SolveError::InvalidTolerance`], see
/// [`Tolerance`]), the step of [`Rk4`] finite and positive
/// ([`SolveError::InvalidStep`]) and `y0` finite ([`SolveError::NonFinite`] at the
/// start time).
pub fn solve_ivp<T, Y, R, F, A, M>(
    f: F,
    span: Range<Quantity<T>>,
    y0: Y,
    method: M,
    tol: Tolerance<A>,
) -> Result<Solution<T, Y>, SolveError<T>>
where
    F: FnMut(Quantity<T>, &Y) -> R,
    Y: VectorSpace + PerDim<T>,
    R: VectorSpace + RhsReturns<<Y as PerDim<T>>::Out>,
    A: VectorSpace + ToleranceFor<Y>,
    M: Method,
{
    solver::solve_ivp(
        f,
        span,
        y0,
        method.kind(),
        method.step_limit(),
        tol.rtol,
        &tol.atol,
        &mut |_, _| {},
    )
}

/// Like [`solve_ivp`], calling `on_step(t, &y)` after every accepted step (API-45).
/// Rejected input (see [`solve_ivp`]) calls neither `f` nor `on_step`.
pub fn solve_ivp_with<T, Y, R, F, A, M, O>(
    f: F,
    span: Range<Quantity<T>>,
    y0: Y,
    method: M,
    tol: Tolerance<A>,
    mut on_step: O,
) -> Result<Solution<T, Y>, SolveError<T>>
where
    F: FnMut(Quantity<T>, &Y) -> R,
    Y: VectorSpace + PerDim<T>,
    R: VectorSpace + RhsReturns<<Y as PerDim<T>>::Out>,
    A: VectorSpace + ToleranceFor<Y>,
    M: Method,
    O: FnMut(Quantity<T>, &Y),
{
    solver::solve_ivp(
        f,
        span,
        y0,
        method.kind(),
        method.step_limit(),
        tol.rtol,
        &tol.atol,
        &mut on_step,
    )
}
