//! A Jacobian given by the user to the Rosenbrock method (issue #23, FR-30). TC-U30.

use std::cell::RefCell;
use std::rc::Rc;

use tensorcheck::units::si::{
    Acceleration, Dimensionless, Frequency, Length, Time, Velocity, m, one, s,
};
use tensorcheck::{Quantity, Vector, VectorSpace};
use tensorcheck_numerics::ode::{
    Jacobian, Method, Rosenbrock, SolveError, Tolerance, solve_ivp, solve_ivp_with,
};

type Osc = (Quantity<Length>, Quantity<Velocity>);
type OscRate = (Quantity<Velocity>, Quantity<Acceleration>);
/// A call of the Jacobian: (t, y, dy) in base units.
type Call = (f64, [f64; 2], [f64; 2]);

const K: f64 = 1.0; // 1/s²
const C: f64 = 1000.0; // 1/s

/// A stiff damped oscillator x'' = −k x − c v, with components of different
/// dimensions: J·dy = (dv, −k dx − c dv).
fn damped(_t: Quantity<Time>, y: &Osc) -> OscRate {
    let k = K * (one / (s * s));
    let c = C * (one / s);
    (y.1, -(k * y.0) - c * y.1)
}
fn damped_jv(_t: Quantity<Time>, _y: &Osc, dy: &Osc) -> OscRate {
    let k = K * (one / (s * s));
    let c = C * (one / s);
    (dy.1, -(k * dy.0) - c * dy.1)
}

/// x(t) for x(0) = x0, v(0) = 0, from the eigenvalues of [[0, 1], [−k, −c]].
fn damped_exact(t: f64, x0: f64) -> f64 {
    let d = (C * C - 4.0 * K).sqrt();
    let (l1, l2) = ((-C + d) / 2.0, (-C - d) / 2.0);
    // x = a e^{l1 t} + b e^{l2 t}, a + b = x0, a l1 + b l2 = 0.
    let a = x0 * l2 / (l2 - l1);
    let b = x0 - a;
    a * (l1 * t).exp() + b * (l2 * t).exp()
}

fn tol(scale: f64) -> Tolerance<Osc> {
    Tolerance {
        rtol: 1e-9,
        atol: (scale * 1e-12 * m, scale * 1e-12 * (m / s)),
    }
}

#[test]
fn the_jacobian_is_evaluated_at_each_step_start_on_the_unit_vectors() {
    // 1e6 s: the times are rounded, but the step stays above the limit 16 ε |t|.
    for t0 in [0.0, 1e6] {
        let calls: RefCell<Vec<Call>> = RefCell::new(Vec::new());
        let jv = |t: Quantity<Time>, y: &Osc, dy: &Osc| -> OscRate {
            calls.borrow_mut().push((
                t.base_value(),
                [y.component(0), y.component(1)],
                [dy.component(0), dy.component(1)],
            ));
            damped_jv(t, y, dy)
        };
        let y0: Osc = (1.0 * m, 0.0 * (m / s));
        let mut points = vec![(t0, [1.0, 0.0])];
        let sol = solve_ivp_with(
            damped,
            (t0 * s)..((t0 + 1.0) * s),
            y0,
            Rosenbrock {
                jacobian: Jacobian::given(jv),
            },
            tol(1.0),
            |t, y: &Osc| points.push((t.base_value(), [y.component(0), y.component(1)])),
        )
        .unwrap();
        let calls = calls.into_inner();
        // One product per unit vector, per attempted step (accepted or rejected).
        assert_eq!(
            calls.len(),
            2 * (sol.stats.steps + sol.stats.rejected),
            "{t0}"
        );
        for (k, (t, y, dy)) in calls.iter().enumerate() {
            let unit = if k % 2 == 0 { [1.0, 0.0] } else { [0.0, 1.0] };
            assert_eq!(*dy, unit, "{t0}: call {k}");
            // (t, y) is the start of a step: the initial point or an accepted point,
            // exactly.
            assert!(
                points.iter().any(|(tp, yp)| tp == t && yp == y),
                "{t0}: call {k} at t = {t}"
            );
        }
    }
}

#[test]
fn any_fnmut_accepted_by_given_can_be_used() {
    // A closure that owns non-Copy state and mutates it (PR #52, second review).
    let count = Rc::new(RefCell::new(0usize));
    let seen = Rc::clone(&count);
    let mut log: Vec<f64> = Vec::new();
    let jv = move |t: Quantity<Time>, y: &Osc, dy: &Osc| -> OscRate {
        log.push(t.base_value());
        *seen.borrow_mut() = log.len();
        damped_jv(t, y, dy)
    };
    let y0: Osc = (1.0 * m, 0.0 * (m / s));
    let method = Rosenbrock {
        jacobian: Jacobian::given(jv),
    };
    // Also through the step limit and the variant with on_step.
    let sol = solve_ivp_with(
        damped,
        (0.0 * s)..(1.0 * s),
        y0,
        method.max_steps(10_000),
        tol(1.0),
        |_, _| {},
    )
    .unwrap();
    assert_eq!(*count.borrow(), 2 * (sol.stats.steps + sol.stats.rejected));
}

#[test]
fn the_given_matrix_is_the_one_used() {
    // A wrong Jacobian (zero) must change the solve; the right one matches the exact
    // solution.
    let y0: Osc = (1.0 * m, 0.0 * (m / s));
    let zero = |_t: Quantity<Time>, _y: &Osc, dy: &Osc| -> OscRate {
        (0.0 * dy.1, 0.0 * (dy.1 / (1.0 * s)))
    };
    let right = solve_ivp(
        damped,
        (0.0 * s)..(1.0 * s),
        y0,
        Rosenbrock {
            jacobian: Jacobian::given(damped_jv),
        },
        tol(1.0),
    )
    .unwrap();
    let wrong = solve_ivp(
        damped,
        (0.0 * s)..(1.0 * s),
        y0,
        Rosenbrock {
            jacobian: Jacobian::given(zero),
        },
        tol(1.0),
    );
    let exact = damped_exact(1.0, 1.0);
    assert!((right.y.0.value_in(m) - exact).abs() < 1e-7 * exact.abs());
    match wrong {
        Ok(w) => assert!(
            w.stats.steps + w.stats.rejected > 2 * (right.stats.steps + right.stats.rejected),
            "zero Jacobian: {:?}, right: {:?}",
            w.stats,
            right.stats
        ),
        Err(e) => panic!("zero Jacobian failed: {e}"),
    }
}

#[test]
fn the_solution_is_right_at_any_scale_of_the_state() {
    // 1e-310 is subnormal.
    for scale in [1.0, 1e-300, 1e300, 1e-310] {
        let y0: Osc = (scale * m, 0.0 * (m / s));
        let sol = solve_ivp(
            damped,
            (0.0 * s)..(1.0 * s),
            y0,
            Rosenbrock {
                jacobian: Jacobian::given(damped_jv),
            },
            tol(scale),
        )
        .unwrap();
        let exact = damped_exact(1.0, scale);
        let got = sol.y.0.value_in(m);
        assert!(got.is_finite(), "{scale}");
        assert!(
            ((got - exact) / exact).abs() < 1e-6,
            "{scale}: {got:e} vs {exact:e}"
        );
    }
}

#[test]
fn backward_integration_matches_the_numerical_jacobian() {
    // From the exact state at t = 1 s back to t = 0.
    let x1 = damped_exact(1.0, 1.0);
    let d = (C * C - 4.0 * K).sqrt();
    let (l1, l2) = ((-C + d) / 2.0, (-C - d) / 2.0);
    let a = l2 / (l2 - l1);
    let b = 1.0 - a;
    let v1 = a * l1 * l1.exp() + b * l2 * l2.exp();
    let y1: Osc = (x1 * m, v1 * (m / s));
    let given = solve_ivp(
        damped,
        (1.0 * s)..(0.0 * s),
        y1,
        Rosenbrock {
            jacobian: Jacobian::given(damped_jv),
        },
        tol(1.0),
    )
    .unwrap();
    let numerical = solve_ivp(
        damped,
        (1.0 * s)..(0.0 * s),
        y1,
        Rosenbrock {
            jacobian: Jacobian::Numerical,
        },
        tol(1.0),
    )
    .unwrap();
    assert_eq!(given.t.value_in(s), 0.0);
    let (xg, xn) = (given.y.0.value_in(m), numerical.y.0.value_in(m));
    assert!(((xg - xn) / xn).abs() < 1e-9, "{xg} vs {xn}");
    // The fast mode b e^{l2 t} (b ≈ −1e-6) is e^{−999} at t = 1 and cannot be
    // recovered backward, so x(0) is 1 m only up to |b|.
    assert!((xg - 1.0).abs() < 1e-5, "{:?}", given.y);
}

#[test]
fn a_non_finite_jacobian_fails_at_the_step_start() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        // Only the entry (1, 0) is bad.
        let jv = move |_t: Quantity<Time>, _y: &Osc, dy: &Osc| -> OscRate {
            (dy.1, bad * (dy.0 / (1.0 * (s * s))))
        };
        let y0: Osc = (1.0 * m, 0.0 * (m / s));
        let r = solve_ivp(
            damped,
            (0.5 * s)..(1.0 * s),
            y0,
            Rosenbrock {
                jacobian: Jacobian::given(jv),
            },
            tol(1.0),
        );
        assert_eq!(
            r.unwrap_err(),
            SolveError::NonFinite { t: 0.5 * s },
            "{bad}"
        );
    }
}

#[test]
fn the_step_limit_applies_with_a_given_jacobian() {
    let y0: Osc = (1.0 * m, 0.0 * (m / s));
    let method = Rosenbrock {
        jacobian: Jacobian::given(damped_jv),
    };
    let r = solve_ivp(
        damped,
        (0.0 * s)..(1.0 * s),
        y0,
        method.max_steps(3),
        tol(1.0),
    );
    assert!(matches!(r, Err(SolveError::MaxSteps { .. })), "{r:?}");
    let r = solve_ivp(
        damped,
        (0.0 * s)..(1.0 * s),
        y0,
        method.max_steps(0),
        tol(1.0),
    );
    assert_eq!(r.unwrap_err(), SolveError::MaxSteps { t: 0.0 * s });
}

tensorcheck::space!(Species, 3, Dimensionless);
tensorcheck::chart!(Conc on Species { Y1: Dimensionless, Y2: Dimensionless, Y3: Dimensionless });

#[test]
fn robertson_matches_the_numerical_jacobian_with_fewer_evaluations() {
    let rhs = |_t: Quantity<Time>, y: &Vector<Dimensionless, Conc>| -> Vector<Frequency, Conc> {
        let [y1, y2, y3] = y.components_in(one);
        Vector::from_components(
            [
                -0.04 * y1 + 1e4 * y2 * y3,
                0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2 * y2,
                3e7 * y2 * y2,
            ],
            one / s,
        )
    };
    let jv = |_t: Quantity<Time>,
              y: &Vector<Dimensionless, Conc>,
              dy: &Vector<Dimensionless, Conc>|
     -> Vector<Frequency, Conc> {
        let [_, y2, y3] = y.components_in(one);
        let [d1, d2, d3] = dy.components_in(one);
        Vector::from_components(
            [
                -0.04 * d1 + 1e4 * (y3 * d2 + y2 * d3),
                0.04 * d1 - 1e4 * (y3 * d2 + y2 * d3) - 6e7 * y2 * d2,
                6e7 * y2 * d2,
            ],
            one / s,
        )
    };
    let y0 = Vector::from_components([1.0, 0.0, 0.0], one);
    let tol = Tolerance {
        rtol: 1e-8,
        atol: Vector::from_components([1e-14, 1e-20, 1e-14], one),
    };
    let span = (0.0 * s)..(1e11 * s);
    let numerical = solve_ivp(
        rhs,
        span.clone(),
        y0,
        Rosenbrock {
            jacobian: Jacobian::Numerical,
        },
        tol,
    )
    .unwrap();
    let given = solve_ivp(
        rhs,
        span,
        y0,
        Rosenbrock {
            jacobian: Jacobian::given(jv),
        },
        tol,
    )
    .unwrap();
    let (yn, yg) = (numerical.y.components_in(one), given.y.components_in(one));
    for i in 0..3 {
        assert!(
            ((yg[i] - yn[i]) / yn[i]).abs() < 1e-7,
            "{i}: {yg:?} vs {yn:?}"
        );
    }
    // Central differences cost 2n = 6 evaluations of f per attempted step.
    assert!(
        given.stats.evaluations < numerical.stats.evaluations,
        "{:?} vs {:?}",
        given.stats,
        numerical.stats
    );
}
