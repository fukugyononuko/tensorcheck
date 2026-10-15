//! TC-U28: invalid input is rejected before anything is evaluated (Codex review
//! 2026-10-14, findings 3 and 4, and the hardening of quadrature counts).

use core::cell::Cell;

use tensorcheck::Quantity;
use tensorcheck::units::si::{Dimensionless, Length, Time, Velocity, m, one, s};
use tensorcheck_numerics::ode::{
    DormandPrince54, Jacobian, Rk4, Rosenbrock, SolveError, solve_ivp, solve_ivp_with,
};
use tensorcheck_numerics::quad::{
    GaussKronrod, GaussLegendre, QuadError, Romberg, Simpson, Tolerance, Trapezoid, integrate,
};

type Q = Quantity<Dimensionless>;
type State = (Quantity<Length>, Quantity<Velocity>);

const BAD: [f64; 4] = [-1e-6, f64::NAN, f64::INFINITY, f64::NEG_INFINITY];

/// Integrates x² on [0, 1] and returns the result with the number of integrand calls.
fn quad<M: tensorcheck_numerics::quad::Method>(
    method: M,
    range: core::ops::Range<Q>,
    rtol: f64,
    atol: f64,
) -> (Result<f64, QuadError>, usize) {
    let calls = Cell::new(0);
    let r = integrate(
        |x: Q| {
            calls.set(calls.get() + 1);
            x * x
        },
        range,
        method,
        Tolerance {
            rtol,
            atol: atol * one,
        },
    );
    (r.map(|i| i.value.value_in(one)), calls.get())
}

/// Integrates a vector-valued integrand (x, x²) with per-component `atol`.
fn quad_vec(rtol: f64, atol: (f64, f64)) -> (Result<(), QuadError>, usize) {
    let calls = Cell::new(0);
    let r = integrate(
        |x: Quantity<Time>| {
            calls.set(calls.get() + 1);
            (
                Quantity::of(x.value_in(s), m),
                Quantity::of(x.value_in(s), one),
            )
        },
        (0.0 * s)..(1.0 * s),
        GaussKronrod,
        Tolerance {
            rtol,
            atol: (Quantity::of(atol.0, m * s), Quantity::of(atol.1, s)),
        },
    );
    (r.map(|_| ()), calls.get())
}

#[test]
fn quadrature_rejects_invalid_tolerances() {
    let r = (0.0 * one)..(1.0 * one);
    for bad in BAD {
        for method_is_gk in [true, false] {
            let (res, calls) = if method_is_gk {
                quad(GaussKronrod, r.clone(), bad, 1e-12)
            } else {
                quad(Romberg, r.clone(), 1e-8, bad)
            };
            assert_eq!(res, Err(QuadError::InvalidTolerance), "{bad}");
            assert_eq!(calls, 0, "{bad}");
        }
        // One bad component of a vector tolerance is enough.
        assert_eq!(
            quad_vec(1e-8, (1e-12, bad)),
            (Err(QuadError::InvalidTolerance), 0)
        );
        assert_eq!(
            quad_vec(1e-8, (bad, 1e-12)),
            (Err(QuadError::InvalidTolerance), 0)
        );
    }
    // Before the fix a negative atol made the scaled error negative and a coarse
    // trapezoid rule was accepted as Ok with an error far above any tolerance.
    let (res, calls) = quad(Trapezoid { n: 2 }, r.clone(), 1e-12, -1e-9);
    assert_eq!((res, calls), (Err(QuadError::InvalidTolerance), 0));
    // Policy for zeros: rtol = 0 is fine with atol > 0, atol = 0 with rtol > 0, but
    // not both zero for a component.
    assert!(quad(GaussKronrod, r.clone(), 0.0, 1e-12).0.is_ok());
    assert!(quad(GaussKronrod, r.clone(), 1e-12, 0.0).0.is_ok());
    assert_eq!(
        quad(GaussKronrod, r.clone(), 0.0, 0.0),
        (Err(QuadError::InvalidTolerance), 0)
    );
    assert_eq!(
        quad_vec(0.0, (1e-12, 0.0)),
        (Err(QuadError::InvalidTolerance), 0)
    );
}

#[test]
fn quadrature_rejects_unsupported_counts_and_ranges() {
    let r = (0.0 * one)..(1.0 * one);
    let tol = (1e-3, 1e-6);
    for (res, n) in [
        (quad(Trapezoid { n: 0 }, r.clone(), tol.0, tol.1), 0),
        (quad(Simpson { n: 0 }, r.clone(), tol.0, tol.1), 0),
        (quad(GaussLegendre { n: 0 }, r.clone(), tol.0, tol.1), 0),
        // 2n, an even n and n + 1 overflow.
        (
            quad(
                Trapezoid {
                    n: usize::MAX / 2 + 1,
                },
                r.clone(),
                tol.0,
                tol.1,
            ),
            usize::MAX / 2 + 1,
        ),
        (
            quad(Simpson { n: usize::MAX }, r.clone(), tol.0, tol.1),
            usize::MAX,
        ),
        (
            quad(
                Simpson {
                    n: usize::MAX / 2 + 1,
                },
                r.clone(),
                tol.0,
                tol.1,
            ),
            usize::MAX / 2 + 1,
        ),
        (
            quad(GaussLegendre { n: usize::MAX }, r.clone(), tol.0, tol.1),
            usize::MAX,
        ),
    ] {
        assert_eq!(res, (Err(QuadError::InvalidCount { n }), 0));
    }
    // The smallest counts work: one interval, Simpson rounded up to two, one point.
    assert!(quad(Trapezoid { n: 1 }, r.clone(), 1.0, 1.0).0.is_ok());
    assert!(quad(Simpson { n: 1 }, r.clone(), 1e-12, 1e-12).0.is_ok());
    assert!(quad(GaussLegendre { n: 1 }, r.clone(), 1.0, 1.0).0.is_ok());
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            quad(GaussKronrod, (0.0 * one)..(bad * one), tol.0, tol.1),
            (Err(QuadError::InvalidRange), 0)
        );
        assert_eq!(
            quad(Romberg, (bad * one)..(1.0 * one), tol.0, tol.1),
            (Err(QuadError::InvalidRange), 0)
        );
    }
}

/// The harmonic oscillator, counting right-hand-side and `on_step` calls.
fn ode<M: tensorcheck_numerics::ode::Method>(
    method: M,
    span: core::ops::Range<Quantity<Time>>,
    y0: State,
    rtol: f64,
    atol: (f64, f64),
) -> (Result<(), SolveError<Time>>, usize, usize) {
    let (calls, steps) = (Cell::new(0), Cell::new(0));
    let r = solve_ivp_with(
        |_t: Quantity<Time>, y: &State| {
            calls.set(calls.get() + 1);
            (y.1, -y.0 / (1.0 * (s * s)))
        },
        span,
        y0,
        method,
        Tolerance {
            rtol,
            atol: (Quantity::of(atol.0, m), Quantity::of(atol.1, m / s)),
        },
        |_, _| steps.set(steps.get() + 1),
    );
    (r.map(|_| ()), calls.get(), steps.get())
}

fn y0() -> State {
    (1.0 * m, 0.0 * (m / s))
}

#[test]
fn ode_rejects_invalid_tolerances() {
    let span = (0.0 * s)..(1.0 * s);
    for bad in BAD {
        for (rtol, atol) in [
            (bad, (1e-9, 1e-9)),
            (1e-6, (bad, 1e-9)),
            (1e-6, (1e-9, bad)),
        ] {
            let want = (Err(SolveError::InvalidTolerance), 0, 0);
            assert_eq!(ode(DormandPrince54, span.clone(), y0(), rtol, atol), want);
            assert_eq!(
                ode(
                    Rosenbrock {
                        jacobian: Jacobian::Numerical
                    },
                    span.clone(),
                    y0(),
                    rtol,
                    atol
                ),
                want
            );
            // The same contract for the fixed-step method, which does not use it.
            assert_eq!(ode(Rk4 { h: 0.1 }, span.clone(), y0(), rtol, atol), want);
        }
    }
    assert_eq!(
        ode(DormandPrince54, span.clone(), y0(), 0.0, (1e-9, 0.0)),
        (Err(SolveError::InvalidTolerance), 0, 0)
    );
    assert!(
        ode(DormandPrince54, span, y0(), 0.0, (1e-9, 1e-9))
            .0
            .is_ok()
    );
}

#[test]
fn rk4_rejects_steps_that_are_not_positive() {
    // Before the fix h = 0 or h < 0 ran 100,000 meaningless steps (each reported to
    // on_step) and then returned MaxSteps.
    for h in [0.0, -0.1, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for span in [(0.0 * s)..(1.0 * s), (1.0 * s)..(0.0 * s)] {
            let (res, calls, steps) = ode(Rk4 { h }, span, y0(), 1e-6, (1e-9, 1e-9));
            assert!(matches!(res, Err(SolveError::InvalidStep { .. })), "{h}");
            assert_eq!((calls, steps), (0, 0), "{h}");
        }
    }
    // Backward integration with a positive magnitude: x(t) = cos t from t = 1 to 0.
    let start = (1f64.cos() * m, -1f64.sin() * (m / s));
    let r = solve_ivp(
        |_t: Quantity<Time>, y: &State| (y.1, -y.0 / (1.0 * (s * s))),
        (1.0 * s)..(0.0 * s),
        start,
        Rk4 { h: 0.01 },
        Tolerance {
            rtol: 1e-6,
            atol: (1e-9 * m, 1e-9 * (m / s)),
        },
    )
    .unwrap();
    assert!((r.y.0.value_in(m) - 1.0).abs() < 1e-9);
    assert_eq!(r.t.value_in(s), 0.0);
}

#[test]
fn ode_rejects_non_finite_spans_and_states() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for span in [(0.0 * s)..(bad * s), (bad * s)..(1.0 * s)] {
            let r = ode(DormandPrince54, span, y0(), 1e-6, (1e-9, 1e-9));
            assert_eq!(r, (Err(SolveError::InvalidSpan), 0, 0), "{bad}");
        }
        // A non-finite initial state fails at the start time, even on an empty interval
        // (which used to return it as the solution).
        for span in [(0.0 * s)..(1.0 * s), (0.5 * s)..(0.5 * s)] {
            let start = span.start;
            let r = ode(
                DormandPrince54,
                span,
                (bad * m, 0.0 * (m / s)),
                1e-6,
                (1e-9, 1e-9),
            );
            assert_eq!(r, (Err(SolveError::NonFinite { t: start }), 0, 0), "{bad}");
        }
    }
    // Rejected input has no failure time.
    assert!(SolveError::<Time>::InvalidSpan.time().is_none());
    assert!(SolveError::<Time>::InvalidStep { h: 0.0 }.time().is_none());
    assert!(SolveError::<Time>::InvalidTolerance.time().is_none());
}
