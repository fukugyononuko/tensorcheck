//! Phase 7: explicit ODE solvers. TC-U23, TC-I06, TC-I08 (geodesic), TC-U27.

use tensorcheck::geometry::christoffel;
use tensorcheck::labels::{a, b, c};
use tensorcheck::spaces::euclid3::Cartesian;
use tensorcheck::spaces::sphere2::{Angles, SphereMetric};
use tensorcheck::units::si::{Dimensionless, Frequency, Length, Time, Velocity, m, one, rad, s};
use tensorcheck::{Point, Quantity, Vector, metric_tensor};
use tensorcheck_numerics::ode::{
    DormandPrince54, Rk4, SolveError, Tolerance, solve_ivp, solve_ivp_with,
};

type D0 = Dimensionless;

/// Harmonic oscillator x'' = −x (ω = 1/s): state (x, v); exact x = cos t.
fn oscillator(
    _t: Quantity<Time>,
    y: &(Quantity<Length>, Quantity<Velocity>),
) -> (
    Quantity<Velocity>,
    Quantity<tensorcheck::units::si::Acceleration>,
) {
    let omega2 = 1.0 * (one / (s * s));
    (
        y.1,
        -(y.0 * Quantity::of(omega2.value_in(one / (s * s)), one / (s * s))),
    )
}

fn osc_error(method: Rk4) -> f64 {
    let y0 = (1.0 * m, 0.0 * (m / s));
    let atol = (1e-12 * m, 1e-12 * (m / s));
    let sol = solve_ivp(
        oscillator,
        (0.0 * s)..(1.0 * s),
        y0,
        method,
        Tolerance { rtol: 1e-12, atol },
    )
    .unwrap();
    (sol.y.0.value_in(m) - 1f64.cos()).abs()
}

#[test]
fn rk4_is_fourth_order() {
    // TC-U23: halving h divides the error by 16.
    let order = (osc_error(Rk4 { h: 0.1 }) / osc_error(Rk4 { h: 0.05 })).log2();
    assert!((order - 4.0).abs() <= 0.2, "order {order}");
}

#[test]
fn dormand_prince_meets_its_tolerance() {
    // TC-U23: local error control; the global error may be a few times the
    // tolerance (DD-13), so the test allows 10 × tol.
    for rtol in [1e-6, 1e-9, 1e-12] {
        let y0 = (1.0 * m, 0.0 * (m / s));
        let atol = (rtol * m, rtol * (m / s));
        let sol = solve_ivp(
            oscillator,
            (0.0 * s)..(10.0 * s),
            y0,
            DormandPrince54,
            Tolerance { rtol, atol },
        )
        .unwrap();
        let err = (sol.y.0.value_in(m) - 10f64.cos()).abs();
        assert!(err <= 10.0 * rtol, "rtol {rtol}: error {err}");
        assert!(sol.stats.steps > 0 && sol.stats.evaluations > 6 * sol.stats.steps);
    }
    // Backwards in time.
    let y0 = (10f64.cos() * m, -10f64.sin() * (m / s));
    let sol = solve_ivp(
        oscillator,
        (10.0 * s)..(0.0 * s),
        y0,
        DormandPrince54,
        Tolerance {
            rtol: 1e-10,
            atol: (1e-10 * m, 1e-10 * (m / s)),
        },
    )
    .unwrap();
    assert!((sol.y.0.value_in(m) - 1.0).abs() <= 1e-9);
}

#[test]
fn small_amplitude_pendulum() {
    // TC-I06 / CC-07: θ'' = −(g/L) sin θ with g/L = 4 /s²; for θ0 = 1e-3 rad,
    // θ(t) ≈ θ0 cos(2t) (relative error O(θ0²) ≈ 1e-7 from the nonlinearity).
    let rhs = |_t: Quantity<Time>, y: &(Quantity<D0>, Quantity<Frequency>)| {
        let th = y.0.value_in(rad);
        (y.1, Quantity::of(-4.0 * th.sin(), rad / (s * s)))
    };
    let y0 = (1e-3 * rad, 0.0 * (rad / s));
    let mut samples = 0;
    let sol = solve_ivp_with(
        rhs,
        (0.0 * s)..(3.0 * s),
        y0,
        DormandPrince54,
        Tolerance {
            rtol: 1e-10,
            atol: (1e-14 * rad, 1e-14 * (rad / s)),
        },
        |_, _| samples += 1,
    )
    .unwrap();
    let want = 1e-3 * (2.0f64 * 3.0).cos();
    assert!((sol.y.0.value_in(rad) - want).abs() <= 1e-3 * 1e-6);
    assert_eq!(samples, sol.stats.steps);
}

#[test]
fn geodesics_on_the_sphere_are_great_circles() {
    // TC-I08 / CC-09: dx/dλ = u, du/dλ = −Γ^a_bc u^b u^c with λ in seconds.
    let a_r = 2.0;
    let g = SphereMetric::new(a_r * m);
    type U = Vector<Velocity, Angles>;
    let rhs = |_l: Quantity<Time>, y: &(Point<Angles>, U)| {
        let gamma = christoffel(&g, &y.0).unwrap();
        let acc = -(gamma.idx((a, b, c)) * y.1.idx(b) * y.1.idx(c)).into_tensor();
        (y.1, acc)
    };
    let p0 = Point::<Angles>::new((0.7 * rad, 0.2 * rad));
    let u0: U = Vector::from_parts((0.3 * (rad / s), 0.5 * (rad / s)));
    let tol = Tolerance {
        rtol: 1e-11,
        atol: (
            Point::<Angles>::new((1e-12 * rad, 1e-12 * rad)),
            Vector::from_parts((1e-12 * (rad / s), 1e-12 * (rad / s))),
        ),
    };
    let speed2 = |p: &Point<Angles>, u: &U| {
        let gm = metric_tensor(&g, p);
        (gm.idx((a, b)) * u.idx(a) * u.idx(b))
            .into_tensor()
            .value_in(m * m / (s * s))
    };
    let embed = |p: &Point<Angles>| {
        let (th, ph) = (p.get(0), p.get(1));
        [th.sin() * ph.cos(), th.sin() * ph.sin(), th.cos()]
    };
    let x0 = embed(&p0);
    // Normal of the great circle: x0 × (dx/dλ at λ = 0).
    let (th, ph) = (0.7f64, 0.2f64);
    let dx = [
        th.cos() * ph.cos() * 0.3 - th.sin() * ph.sin() * 0.5,
        th.cos() * ph.sin() * 0.3 + th.sin() * ph.cos() * 0.5,
        -th.sin() * 0.3,
    ];
    let nrm = [
        x0[1] * dx[2] - x0[2] * dx[1],
        x0[2] * dx[0] - x0[0] * dx[2],
        x0[0] * dx[1] - x0[1] * dx[0],
    ];
    let v0 = speed2(&p0, &u0);
    let mut worst_plane: f64 = 0.0;
    let mut worst_speed: f64 = 0.0;
    let sol = solve_ivp_with(
        rhs,
        (0.0 * s)..(4.0 * s),
        (p0, u0),
        DormandPrince54,
        tol,
        |_, y| {
            let x = embed(&y.0);
            let off = (x[0] * nrm[0] + x[1] * nrm[1] + x[2] * nrm[2]).abs();
            worst_plane = worst_plane.max(off);
            worst_speed = worst_speed.max((speed2(&y.0, &y.1) - v0).abs() / v0);
        },
    )
    .unwrap();
    // Allow 10 × rtol for both (local error control, DD-13).
    assert!(
        worst_plane <= 1e-10,
        "off the great circle by {worst_plane}"
    );
    assert!(worst_speed <= 1e-10, "speed changed by {worst_speed}");
    assert!(sol.stats.steps > 5);
    let _ = Cartesian;
}

#[test]
fn failures_carry_the_time() {
    // TC-U27. y' = y² from y(0) = 1 blows up at t = 1.
    let rhs = |_t: Quantity<Time>, y: &Quantity<D0>| Quantity::of(y.value_in(one).powi(2), one / s);
    let r = solve_ivp(
        rhs,
        (0.0 * s)..(2.0 * s),
        1.0 * one,
        DormandPrince54,
        Tolerance {
            rtol: 1e-8,
            atol: 1e-8 * one,
        },
    );
    let e = r.expect_err("must fail near the blow-up");
    assert!(
        matches!(
            e,
            SolveError::StepTooSmall { .. } | SolveError::NonFinite { .. }
        ),
        "{e}"
    );
    let tf = e.time().expect("a solve failure has a time").value_in(s);
    assert!((tf - 1.0).abs() < 1e-3, "failed at t = {tf}");

    // MaxSteps: a fixed step far too small for the interval.
    let decay = |_t: Quantity<Time>, y: &Quantity<D0>| -*y / (1.0 * s);
    let r = solve_ivp(
        decay,
        (0.0 * s)..(1.0 * s),
        1.0 * one,
        Rk4 { h: 1e-6 },
        Tolerance {
            rtol: 1e-6,
            atol: 1e-6 * one,
        },
    );
    assert!(matches!(r, Err(SolveError::MaxSteps { .. })));

    // NonFinite: the right-hand side returns NaN.
    let bad = |_t: Quantity<Time>, _y: &Quantity<D0>| Quantity::of(f64::NAN, one / s);
    let r = solve_ivp(
        bad,
        (0.0 * s)..(1.0 * s),
        1.0 * one,
        Rk4 { h: 0.1 },
        Tolerance {
            rtol: 1e-6,
            atol: 1e-6 * one,
        },
    );
    assert!(matches!(r, Err(SolveError::NonFinite { .. })));
}
