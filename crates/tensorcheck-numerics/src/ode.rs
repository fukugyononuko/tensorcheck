//! Ordinary differential equations that keep dimensions (API-13, API-40 to API-45,
//! API-57 to API-59).
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
use core::marker::PhantomData;
use core::ops::Range;

use tensorcheck::utility::diagnostics::{JacobianReturns, RhsReturns, ToleranceFor};
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
pub struct Rosenbrock<J = Jacobian> {
    /// How the Jacobian `∂f/∂y` is obtained: [`Jacobian::Numerical`] or
    /// [`Jacobian::given`].
    pub jacobian: J,
}

/// How the Jacobian is obtained for [`Rosenbrock`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Jacobian {
    /// Central differences at every step (DD-15).
    Numerical,
}

impl Jacobian {
    /// The Jacobian given by the user (API-60, FR-30), as the product
    /// `jv(t, &y, &dy)` = `J(t, y)·dy`, the sum over `j` of `∂f/∂y_j · dy_j`. `dy`
    /// has the type of the state and the result the type of the right-hand side
    /// (ERR-23), so every entry `(i, j)` must have the dimension `[f_i]/[y_j]`: a
    /// term of another dimension cannot be added to the sum (ERR-01). This holds
    /// when the state and the result are types of this library (quantities,
    /// tensors, points and tuples of them). A state type of your own implements
    /// `VectorSpace`, and implementing `JacobianReturns` for it replaces the check,
    /// as `RhsReturns` does for the right-hand side (issue #55).
    ///
    /// `jv` must be linear in `dy`; the solver builds the matrix from the products
    /// with the unit vectors of the state in base units, at the start of every
    /// step. `∂f/∂t` is still found by central differences (DD-15). A product with
    /// a NaN or infinite component fails the solve with [`SolveError::NonFinite`]
    /// at the step's start time. [`Stats::evaluations`] counts `f`, not `jv`.
    ///
    /// ```
    /// use tensorcheck::Quantity;
    /// use tensorcheck::units::si::{Dimensionless, Frequency, Time, one, s};
    /// use tensorcheck_numerics::ode::{Jacobian, Rosenbrock, Tolerance, solve_ivp};
    ///
    /// // dy/dt = −k y with k = 1000 /s: J·dy = −k dy.
    /// let k = 1000.0 * (one / s);
    /// let rhs = |_t: Quantity<Time>, y: &Quantity<Dimensionless>| -> Quantity<Frequency> { -k * *y };
    /// let jv = |_t: Quantity<Time>, _y: &Quantity<Dimensionless>, dy: &Quantity<Dimensionless>| -> Quantity<Frequency> { -k * *dy };
    /// let method = Rosenbrock { jacobian: Jacobian::given(jv) };
    /// let sol = solve_ivp(rhs, (0.0 * s)..(0.01 * s), 1.0 * one, method, Tolerance { rtol: 1e-8, atol: 1e-12 * one }).unwrap();
    /// assert!((sol.y.value_in(one) - (-10.0f64).exp()).abs() < 1e-8);
    /// ```
    pub fn given<T, Y, R, J>(jv: J) -> GivenJacobian<J, T, Y, R>
    where
        J: FnMut(Quantity<T>, &Y, &Y) -> R,
        Y: VectorSpace + PerDim<T>,
        R: VectorSpace + JacobianReturns<<Y as PerDim<T>>::Out>,
    {
        GivenJacobian {
            jv,
            types: PhantomData,
        }
    }
}

/// A Jacobian given by the user, made by [`Jacobian::given`].
pub struct GivenJacobian<J, T, Y, R> {
    jv: J,
    types: PhantomData<fn(Quantity<T>, Y) -> R>,
}
impl<J: Clone, T, Y, R> Clone for GivenJacobian<J, T, Y, R> {
    fn clone(&self) -> Self {
        GivenJacobian {
            jv: self.jv.clone(),
            types: PhantomData,
        }
    }
}
impl<J: Copy, T, Y, R> Copy for GivenJacobian<J, T, Y, R> {}
impl<J, T, Y, R> fmt::Debug for GivenJacobian<J, T, Y, R> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("GivenJacobian")
    }
}

mod sealed {
    /// ERR-24: only the Jacobians of this crate, whose dimensions are checked, can
    /// be used; a matrix filled outside would skip the checks (PR #52 review).
    #[diagnostic::on_unimplemented(
        message = "this is not a Jacobian for this solve",
        label = "not `Jacobian::Numerical` or `Jacobian::given`",
        note = "use `Jacobian::Numerical`, or `Jacobian::given(|t, y, dy| …)` with the time and state types of the right-hand side"
    )]
    pub trait Sealed {}
}
impl sealed::Sealed for Jacobian {}
impl<J, T, Y, R> sealed::Sealed for GivenJacobian<J, T, Y, R> {}

/// ERR-24: the `jacobian` of [`Rosenbrock`] must be [`Jacobian::Numerical`] or a
/// [`Jacobian::given`] taking the time and state types of the solve. Sealed: it is
/// implemented only for those two.
#[doc(hidden)]
#[diagnostic::on_unimplemented(
    message = "this is not a Jacobian for this solve",
    label = "not `Jacobian::Numerical` or a `Jacobian::given` for this time and state",
    note = "use `Jacobian::Numerical`, or `Jacobian::given(|t, y, dy| …)` with the time and state types of the right-hand side"
)]
pub trait JacobianFor<T, Y>: sealed::Sealed {
    /// Fills `out` (row-major, base units) with `∂f/∂y` at `(t, y)`; `false` means
    /// that the solver uses central differences instead.
    fn fill(&mut self, t: f64, y: &[f64], out: &mut [f64]) -> bool;
}
impl<T, Y> JacobianFor<T, Y> for Jacobian {
    fn fill(&mut self, _: f64, _: &[f64], _: &mut [f64]) -> bool {
        false
    }
}
impl<J, T, Y, R> JacobianFor<T, Y> for GivenJacobian<J, T, Y, R>
where
    J: FnMut(Quantity<T>, &Y, &Y) -> R,
    Y: VectorSpace,
    R: VectorSpace,
{
    fn fill(&mut self, t: f64, y: &[f64], out: &mut [f64]) -> bool {
        let n = Y::LEN;
        let y = Y::from_components(|i| y[i]);
        for j in 0..n {
            // Column j is the product with the unit vector e_j.
            let e = Y::from_components(|i| if i == j { 1.0 } else { 0.0 });
            let col = (self.jv)(Quantity::from_base_components(t), &y, &e);
            for i in 0..n {
                out[i * n + j] = col.component(i);
            }
        }
        true
    }
}

/// An ODE method. Not required to be `Copy`, so that a [`Jacobian::given`] may own
/// and mutate its state (PR #52 review).
pub trait Method: Sized {
    #[doc(hidden)]
    type Jac;
    #[doc(hidden)]
    fn jacobian(self) -> Self::Jac;
    #[doc(hidden)]
    fn kind(&self) -> solver::OdeKind;
    #[doc(hidden)]
    fn step_limit(&self) -> usize {
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
    type Jac = M::Jac;
    fn jacobian(self) -> M::Jac {
        self.method.jacobian()
    }
    fn kind(&self) -> solver::OdeKind {
        self.method.kind()
    }
    fn step_limit(&self) -> usize {
        self.max_steps
    }
}
impl Method for Rk4 {
    type Jac = Jacobian;
    fn jacobian(self) -> Jacobian {
        Jacobian::Numerical
    }
    fn kind(&self) -> solver::OdeKind {
        solver::OdeKind::Rk4(self.h)
    }
}
impl Method for DormandPrince54 {
    type Jac = Jacobian;
    fn jacobian(self) -> Jacobian {
        Jacobian::Numerical
    }
    fn kind(&self) -> solver::OdeKind {
        solver::OdeKind::Dopri5
    }
}
impl<J> Method for Rosenbrock<J> {
    type Jac = J;
    fn jacobian(self) -> J {
        self.jacobian
    }
    fn kind(&self) -> solver::OdeKind {
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

/// The Jacobian of `method` on plain arrays, as the solver calls it.
fn given_jacobian<T, Y, M>(method: M) -> impl FnMut(f64, &[f64], &mut [f64]) -> bool
where
    M: Method,
    M::Jac: JacobianFor<T, Y>,
{
    let mut jac = method.jacobian();
    move |t, y, out| jac.fill(t, y, out)
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
    M::Jac: JacobianFor<T, Y>,
{
    solver::solve_ivp(
        f,
        span,
        y0,
        method.kind(),
        method.step_limit(),
        &mut given_jacobian::<T, Y, M>(method),
        tol.rtol,
        &tol.atol,
        &mut |_, _| {},
        &mut (),
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
    M::Jac: JacobianFor<T, Y>,
    O: FnMut(Quantity<T>, &Y),
{
    solver::solve_ivp(
        f,
        span,
        y0,
        method.kind(),
        method.step_limit(),
        &mut given_jacobian::<T, Y, M>(method),
        tol.rtol,
        &tol.atol,
        &mut on_step,
        &mut (),
    )
}

/// A condition `g(t, y) = 0` to find while solving (API-58, FR-31). `g` returns a
/// quantity of any dimension. By default, every sign change is recorded and the
/// solve goes on.
///
/// `g` must be continuous along the interpolated solution. Only the sign of `g` is
/// examined, so for a discontinuous `g` an event may be reported where `g` is not
/// zero. In every case the event point is either a point where `g` was evaluated
/// as exactly zero, or the end `θ_b` on the side of the step end (after the sign
/// change in the direction of integration, so the earlier time when solving
/// backward) of an interval of `θ` at most 4ε wide across which the sign of `g`
/// changes (for a jump of sign, the interval holds the jump). Here `θ` runs from 0
/// to 1 over a step from `t_n` to `t_n + h`, and `g` is evaluated at, and the event
/// is reported at, the time `t_n + θ h` rounded to `f64` (the step end time when
/// `θ = 1`). Only the width in `θ` is bounded: the times of the two ends can be
/// further apart than 4ε times the step, for example one `f64` spacing
/// (1.2e-4 s) apart near 1e12 s (issue #48). If `g` is evaluated as NaN or
/// infinity, for example exactly at a pole, the solve fails with
/// [`SolveError::NonFinite`].
///
/// ```
/// use tensorcheck::Quantity;
/// use tensorcheck::units::si::{Dimensionless, Frequency, Time, one, s};
/// use tensorcheck_numerics::ode::{DormandPrince54, Event, Tolerance, solve_ivp_events};
///
/// // y = t: stop where y reaches 0.5.
/// let rhs = |_t: Quantity<Time>, _y: &Quantity<Dimensionless>| -> Quantity<Frequency> { 1.0 * (one / s) };
/// let half = Event::new(|_t: Quantity<Time>, y: &Quantity<Dimensionless>| *y - 0.5 * one).terminal();
/// let tol = Tolerance { rtol: 1e-10, atol: 1e-12 * one };
/// let sol = solve_ivp_events(rhs, (0.0 * s)..(1.0 * s), 0.0 * one, DormandPrince54, tol, half).unwrap();
/// assert!((sol.t.value_in(s) - 0.5).abs() < 1e-12);
/// assert_eq!(sol.events.len(), 1);
/// ```
#[derive(Clone, Copy, Debug)]
pub struct Event<G> {
    g: G,
    rising: bool,
    falling: bool,
    terminal: bool,
}
impl<G> Event<G> {
    /// Records every zero crossing of `g` and goes on.
    pub fn new(g: G) -> Self {
        Event {
            g,
            rising: true,
            falling: true,
            terminal: false,
        }
    }
    /// Only crossings from negative to zero or positive.
    pub fn rising(self) -> Self {
        Event {
            falling: false,
            ..self
        }
    }
    /// Only crossings from positive to zero or negative.
    pub fn falling(self) -> Self {
        Event {
            rising: false,
            ..self
        }
    }
    /// Stops the solve at the first event.
    pub fn terminal(self) -> Self {
        Event {
            terminal: true,
            ..self
        }
    }
}

/// Where an event happened (API-59).
pub struct Crossing<T, Y> {
    /// The time of the event.
    pub t: Quantity<T>,
    /// The state at that time, interpolated within the step (DD-12): the
    /// continuous extension of order 4 for [`DormandPrince54`], cubic Hermite
    /// interpolation for [`Rk4`] and [`Rosenbrock`].
    pub y: Y,
}

/// The result of [`solve_ivp_events`] (API-59): the state at the end, or at the
/// event when the event is terminal, and the events in time order.
pub struct EventSolution<T, Y> {
    /// The final time: the end of the interval, or the time of a terminal event.
    pub t: Quantity<T>,
    /// The final state.
    pub y: Y,
    /// Counters.
    pub stats: Stats,
    /// The events found.
    pub events: Vec<Crossing<T, Y>>,
}

impl<T, Y: Clone> Clone for Crossing<T, Y> {
    fn clone(&self) -> Self {
        Crossing {
            t: self.t,
            y: self.y.clone(),
        }
    }
}
impl<T, Y: fmt::Debug> fmt::Debug for Crossing<T, Y> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Crossing")
            .field("t (base units)", &self.t.base_value())
            .field("y", &self.y)
            .finish()
    }
}
impl<T, Y: Clone> Clone for EventSolution<T, Y> {
    fn clone(&self) -> Self {
        EventSolution {
            t: self.t,
            y: self.y.clone(),
            stats: self.stats,
            events: self.events.clone(),
        }
    }
}
impl<T, Y: fmt::Debug> fmt::Debug for EventSolution<T, Y> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("EventSolution")
            .field("t (base units)", &self.t.base_value())
            .field("y", &self.y)
            .field("stats", &self.stats)
            .field("events", &self.events)
            .finish()
    }
}

/// An [`Event`] as the solver sees it, collecting the crossings.
struct Watcher<G, T, Y> {
    event: Event<G>,
    found: Vec<Crossing<T, Y>>,
}
impl<G, T, Y, D> solver::Watch for Watcher<G, T, Y>
where
    G: FnMut(Quantity<T>, &Y) -> Quantity<D>,
    Y: VectorSpace,
{
    const ACTIVE: bool = true;
    fn value(&mut self, t: f64, y: &[f64]) -> f64 {
        (self.event.g)(
            Quantity::from_base_components(t),
            &Y::from_components(|i| y[i]),
        )
        .base_value()
    }
    fn wants(&self, g_before: f64) -> bool {
        (self.event.rising && g_before < 0.0) || (self.event.falling && g_before > 0.0)
    }
    fn terminal(&self) -> bool {
        self.event.terminal
    }
    fn record(&mut self, t: f64, y: &[f64]) {
        self.found.push(Crossing {
            t: Quantity::from_base_components(t),
            y: Y::from_components(|i| y[i]),
        });
    }
}

/// Like [`solve_ivp`], finding where `g(t, y)` of `event` crosses zero (API-59,
/// FR-31). The sign of `g` is checked at the end of every accepted step; in a step
/// where it changes, the time is found on the interpolant of the step (no
/// re-solving; see [`Crossing::y`]). Two sign changes within one step are missed, and a zero exactly
/// at the start is not an event. `g` is assumed continuous; see [`Event`] for what
/// the event point guarantees otherwise. `g` returning NaN or infinity gives
/// [`SolveError::NonFinite`]. Rejected input calls neither `f` nor `g`.
pub fn solve_ivp_events<T, Y, R, F, A, M, G, D>(
    f: F,
    span: Range<Quantity<T>>,
    y0: Y,
    method: M,
    tol: Tolerance<A>,
    event: Event<G>,
) -> Result<EventSolution<T, Y>, SolveError<T>>
where
    F: FnMut(Quantity<T>, &Y) -> R,
    Y: VectorSpace + PerDim<T>,
    R: VectorSpace + RhsReturns<<Y as PerDim<T>>::Out>,
    A: VectorSpace + ToleranceFor<Y>,
    M: Method,
    M::Jac: JacobianFor<T, Y>,
    G: FnMut(Quantity<T>, &Y) -> Quantity<D>,
{
    solve_ivp_events_with(f, span, y0, method, tol, event, |_, _| {})
}

/// Like [`solve_ivp_events`], calling `on_step(t, &y)` after every accepted step
/// (API-59). A terminal event ends the calls with the event point.
pub fn solve_ivp_events_with<T, Y, R, F, A, M, G, D, O>(
    f: F,
    span: Range<Quantity<T>>,
    y0: Y,
    method: M,
    tol: Tolerance<A>,
    event: Event<G>,
    mut on_step: O,
) -> Result<EventSolution<T, Y>, SolveError<T>>
where
    F: FnMut(Quantity<T>, &Y) -> R,
    Y: VectorSpace + PerDim<T>,
    R: VectorSpace + RhsReturns<<Y as PerDim<T>>::Out>,
    A: VectorSpace + ToleranceFor<Y>,
    M: Method,
    M::Jac: JacobianFor<T, Y>,
    G: FnMut(Quantity<T>, &Y) -> Quantity<D>,
    O: FnMut(Quantity<T>, &Y),
{
    let mut watcher = Watcher {
        event,
        found: Vec::new(),
    };
    let sol = solver::solve_ivp(
        f,
        span,
        y0,
        method.kind(),
        method.step_limit(),
        &mut given_jacobian::<T, Y, M>(method),
        tol.rtol,
        &tol.atol,
        &mut on_step,
        &mut watcher,
    )?;
    Ok(EventSolution {
        t: sol.t,
        y: sol.y,
        stats: sol.stats,
        events: watcher.found,
    })
}
