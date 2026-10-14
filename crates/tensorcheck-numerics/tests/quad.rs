//! TC-U22 (adaptive rules), TC-I07 (CC-08), failures of TC-U27.

use tensorcheck::spaces::euclid3::Cartesian;
use tensorcheck::units::si::{Dimensionless, Energy, Force, J, Length, N, m, one, s};
use tensorcheck::{Quantity, Vector};
use tensorcheck_numerics::quad::{
    GaussKronrod, GaussLegendre, QuadError, Romberg, Simpson, Tolerance, Trapezoid, integrate,
};

type Q = Quantity<Dimensionless>;

fn check<M: tensorcheck_numerics::quad::Method>(
    method: M,
    f: fn(f64) -> f64,
    a: f64,
    b: f64,
    exact: f64,
    rtol: f64,
) {
    let tol = Tolerance {
        rtol,
        atol: 1e-300 * one,
    };
    let r = integrate(
        |x: Q| Quantity::of(f(x.value_in(one)), one),
        (a * one)..(b * one),
        method,
        tol,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let got = r.value.value_in(one);
    let actual = (got - exact).abs();
    assert!(
        actual <= rtol * exact.abs(),
        "error {actual} above tolerance {}",
        rtol * exact.abs()
    );
    // The estimate must not understate the error (an understated estimate is unsafe).
    assert!(
        r.error_estimate >= actual,
        "estimate {} below actual {actual}",
        r.error_estimate
    );
}

#[test]
fn adaptive_rules_meet_their_tolerance() {
    let e1 = 1f64.exp() - 1.0;
    check(GaussKronrod, f64::exp, 0.0, 1.0, e1, 1e-12);
    check(Romberg, f64::exp, 0.0, 1.0, e1, 1e-12);
    check(GaussKronrod, f64::sqrt, 0.0, 1.0, 2.0 / 3.0, 1e-10);
    check(Romberg, f64::sqrt, 0.0, 1.0, 2.0 / 3.0, 1e-6);
    let pi = std::f64::consts::PI;
    check(
        GaussKronrod,
        |x| (50.0 * x).sin().powi(2),
        0.0,
        pi,
        pi / 2.0,
        1e-10,
    );
    check(
        Romberg,
        |x| (50.0 * x).sin().powi(2),
        0.0,
        pi,
        pi / 2.0,
        1e-10,
    );
    check(GaussLegendre { n: 12 }, f64::exp, 0.0, 1.0, e1, 1e-12);
    check(Simpson { n: 64 }, f64::exp, 0.0, 1.0, e1, 1e-8);
}

#[test]
fn work_of_a_spring_has_dimension_of_energy() {
    // TC-I07 / CC-08: W = ∫ F dx with F = −k x, k = 200 N/m, from 0 to 0.1 m: −1 J.
    let force =
        |x: Quantity<Length>| -> Quantity<Force> { Quantity::new(-200.0 * x.value_in(m), N) };
    let w = integrate(
        force,
        (0.0 * m)..(0.1 * m),
        GaussKronrod,
        Tolerance {
            rtol: 1e-12,
            atol: 1e-15 * J,
        },
    )
    .unwrap();
    let work: Quantity<Energy> = w.value;
    assert!((work.value_in(J) + 1.0).abs() <= 1e-12);
}

#[test]
fn vector_valued_integrand() {
    // ∫_0^1 (t, t², 1) m dt = (1/2, 1/3, 1) m·s.
    let f = |t: Quantity<tensorcheck::units::si::Time>| -> Vector<Length, Cartesian> {
        let t = t.value_in(s);
        Vector::from_components([t, t * t, 1.0], m)
    };
    let atol: Vector<tensorcheck::Prod<Length, tensorcheck::units::si::Time>, Cartesian> =
        Vector::from_components([1e-14; 3], m * s);
    let r = integrate(
        f,
        (0.0 * s)..(1.0 * s),
        GaussKronrod,
        Tolerance { rtol: 1e-12, atol },
    )
    .unwrap();
    let v = r.value.components_in(m * s);
    for (got, want) in v.iter().zip([0.5, 1.0 / 3.0, 1.0]) {
        assert!((got - want).abs() <= 1e-13);
    }
}

#[test]
fn failures_are_reported_not_hidden() {
    // A coarse rule with a tight tolerance: the estimate exceeds the tolerance.
    let r = integrate(
        |x: Q| x * x,
        (0.0 * one)..(1.0 * one),
        Trapezoid { n: 2 },
        Tolerance {
            rtol: 1e-12,
            atol: 1e-300 * one,
        },
    );
    assert!(matches!(r, Err(QuadError::ToleranceNotMet { .. })), "{r:?}");
    // 1/x at x = 0 is infinite.
    let r = integrate(
        |x: Q| Quantity::of(1.0 / x.value_in(one), one),
        (0.0 * one)..(1.0 * one),
        Trapezoid { n: 4 },
        Tolerance {
            rtol: 1e-6,
            atol: 1e-9 * one,
        },
    );
    assert_eq!(r.err(), Some(QuadError::NonFinite));
}
