//! TC-U26: Ridders' error estimate is not below the actual error on 4 functions ×
//! 200 points (run `ratio_of_estimate_to_actual_error` with --ignored --nocapture
//! to see the measured ratios).
use tensorcheck::Quantity;
use tensorcheck::units::si::{Length, Time, m, s};
use tensorcheck_numerics::diff::{Richardson, derivative};

/// A real function.
type Fun = fn(f64) -> f64;

#[test]
#[ignore]
fn ratio_of_estimate_to_actual_error() {
    let cases: [(&str, Fun, Fun); 4] = [
        (
            "e^t sin t",
            |t| t.exp() * t.sin(),
            |t| t.exp() * (t.sin() + t.cos()),
        ),
        ("sin 3t", |t| (3.0 * t).sin(), |t| 3.0 * (3.0 * t).cos()),
        ("t^5", |t| t.powi(5), |t| 5.0 * t.powi(4)),
        (
            "ln(1+t^2)",
            |t| (1.0 + t * t).ln(),
            |t| 2.0 * t / (1.0 + t * t),
        ),
    ];
    for (name, f, df) in cases {
        let (mut min_ratio, mut max_ratio) = (f64::INFINITY, 0.0f64);
        for k in 0..200 {
            let t0 = -5.0 + 10.0 * (k as f64 + 0.5) / 200.0;
            let d = derivative(
                |t: Quantity<Time>| -> Quantity<Length> { Quantity::new(f(t.value_in(s)), m) },
                t0 * s,
                Richardson,
            );
            let actual = (d.value.value_in(m / s) - df(t0)).abs();
            if actual > 0.0 {
                let r = d.error_estimate / actual;
                min_ratio = min_ratio.min(r);
                max_ratio = max_ratio.max(r);
            }
        }
        println!("{name:12} estimate/actual: min {min_ratio:.2}  max {max_ratio:.1e}");
    }
}

#[test]
fn ridders_estimate_is_not_below_the_actual_error() {
    let cases: [(Fun, Fun); 4] = [
        (|t| t.exp() * t.sin(), |t| t.exp() * (t.sin() + t.cos())),
        (|t| (3.0 * t).sin(), |t| 3.0 * (3.0 * t).cos()),
        (|t| t.powi(5), |t| 5.0 * t.powi(4)),
        (|t| (1.0 + t * t).ln(), |t| 2.0 * t / (1.0 + t * t)),
    ];
    for (f, df) in cases {
        for k in 0..200 {
            let t0 = -5.0 + 10.0 * (k as f64 + 0.5) / 200.0;
            let d = derivative(
                |t: Quantity<Time>| -> Quantity<Length> { Quantity::new(f(t.value_in(s)), m) },
                t0 * s,
                Richardson,
            );
            let actual = (d.value.value_in(m / s) - df(t0)).abs();
            assert!(
                d.error_estimate >= actual,
                "t0 = {t0}: estimate {} below {actual}",
                d.error_estimate
            );
        }
    }
}
