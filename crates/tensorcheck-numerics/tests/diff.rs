//! TC-U26: finite differences.
use tensorcheck::Quantity;
use tensorcheck::Vector;
use tensorcheck::spaces::plane::Cartesian2;
use tensorcheck::units::si::{Dimensionless, Length, Time, Velocity, m, one, s};
use tensorcheck_numerics::diff::{Central, Richardson, derivative, jacobian};

#[test]
fn central_and_richardson_derivatives() {
    let x = |t: Quantity<Time>| -> Quantity<Length> {
        Quantity::new(t.value_in(s).exp() * t.value_in(s).sin(), m)
    };
    let exact = |t: f64| t.exp() * (t.sin() + t.cos());
    for t0 in [0.3, 1.7, 4.0] {
        let c = derivative(x, t0 * s, Central);
        let r = derivative(x, t0 * s, Richardson);
        let (vc, vr): (Quantity<Velocity>, Quantity<Velocity>) = (c.value, r.value);
        let want = exact(t0);
        let ec = (vc.value_in(m / s) - want).abs();
        let er = (vr.value_in(m / s) - want).abs();
        // Central: error ≈ ε^(2/3) relative (rounding/truncation balance), allow 100×.
        assert!(
            ec <= 100.0 * f64::EPSILON.powf(2.0 / 3.0) * want.abs().max(1.0),
            "central {ec}"
        );
        // Ridders: much more accurate, and its estimate is not below the actual error.
        assert!(er <= 1e-9 * want.abs().max(1.0), "richardson {er}");
        assert!(
            r.error_estimate >= er,
            "estimate {} below {er}",
            r.error_estimate
        );
    }
}

#[test]
fn central_difference_is_second_order_in_the_step() {
    // Order 2 checked on the engine formula with explicit steps.
    let f = |x: f64| x.sin();
    let d = |h: f64| (f(1.0 + h) - f(1.0 - h)) / (2.0 * h);
    let e1 = (d(1e-2) - 1f64.cos()).abs();
    let e2 = (d(5e-3) - 1f64.cos()).abs();
    assert!(((e1 / e2).log2() - 2.0).abs() <= 0.2);
}

#[test]
fn jacobian_of_a_vector_function() {
    // f(x, y) = (x² y, sin x): J = [[2xy, x²], [cos x, 0]].
    let f = |v: &Vector<Dimensionless, Cartesian2>| -> Vector<Dimensionless, Cartesian2> {
        let c = v.components_in(one);
        Vector::from_components([c[0] * c[0] * c[1], c[0].sin()], one)
    };
    let v: Vector<Dimensionless, Cartesian2> = Vector::from_components([1.3, -0.4], one);
    let j = jacobian(f, &v);
    let want = [[2.0 * 1.3 * -0.4, 1.3 * 1.3], [1.3f64.cos(), 0.0]];
    for (i, row) in want.iter().enumerate() {
        for (k, w) in row.iter().enumerate() {
            assert!((j.get(i, k) - w).abs() <= 1e-9, "{i}{k}");
        }
    }
}
