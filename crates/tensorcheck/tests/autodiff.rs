//! TC-U09: dual numbers give exact derivatives up to rounding.
use tensorcheck::Scalar;
use tensorcheck::spaces::plane::Plane;
use tensorcheck::utility::autodiff::{Dual, Dual2};

type D2 = Dual2<Plane>;

/// (name, f, f', f'', domain)
#[allow(clippy::type_complexity)]
fn cases() -> Vec<(
    &'static str,
    fn(D2) -> D2,
    fn(f64) -> f64,
    fn(f64) -> f64,
    fn(f64) -> f64,
    (f64, f64),
)> {
    vec![
        (
            "sqrt",
            |x| x.sqrt(),
            f64::sqrt,
            |x| 0.5 / x.sqrt(),
            |x| -0.25 * x.powf(-1.5),
            (0.1, 5.0),
        ),
        (
            "exp",
            |x| x.exp(),
            f64::exp,
            f64::exp,
            f64::exp,
            (-3.0, 3.0),
        ),
        (
            "ln",
            |x| x.ln(),
            f64::ln,
            |x| 1.0 / x,
            |x| -1.0 / (x * x),
            (0.1, 5.0),
        ),
        (
            "sin",
            |x| x.sin(),
            f64::sin,
            f64::cos,
            |x| -x.sin(),
            (-3.0, 3.0),
        ),
        (
            "cos",
            |x| x.cos(),
            f64::cos,
            |x| -x.sin(),
            |x| -x.cos(),
            (-3.0, 3.0),
        ),
        (
            "tan",
            |x| x.tan(),
            f64::tan,
            |x| 1.0 / x.cos().powi(2),
            |x| 2.0 * x.tan() / x.cos().powi(2),
            (-1.2, 1.2),
        ),
        (
            "sinh",
            |x| x.sinh(),
            f64::sinh,
            f64::cosh,
            f64::sinh,
            (-3.0, 3.0),
        ),
        (
            "cosh",
            |x| x.cosh(),
            f64::cosh,
            f64::sinh,
            f64::cosh,
            (-3.0, 3.0),
        ),
        (
            "tanh",
            |x| x.tanh(),
            f64::tanh,
            |x| 1.0 / x.cosh().powi(2),
            |x| -2.0 * x.tanh() / x.cosh().powi(2),
            (-3.0, 3.0),
        ),
        (
            "asin",
            |x| x.asin(),
            f64::asin,
            |x| 1.0 / (1.0 - x * x).sqrt(),
            |x| x / (1.0 - x * x).powf(1.5),
            (-0.9, 0.9),
        ),
        (
            "acos",
            |x| x.acos(),
            f64::acos,
            |x| -1.0 / (1.0 - x * x).sqrt(),
            |x| -x / (1.0 - x * x).powf(1.5),
            (-0.9, 0.9),
        ),
        (
            "atan",
            |x| x.atan(),
            f64::atan,
            |x| 1.0 / (1.0 + x * x),
            |x| -2.0 * x / (1.0 + x * x).powi(2),
            (-3.0, 3.0),
        ),
        (
            "powi3",
            |x| x.powi(3),
            |x| x.powi(3),
            |x| 3.0 * x * x,
            |x| 6.0 * x,
            (-3.0, 3.0),
        ),
        (
            "powf2.5",
            |x| x.powf(2.5),
            |x| x.powf(2.5),
            |x| 2.5 * x.powf(1.5),
            |x| 3.75 * x.powf(0.5),
            (0.1, 3.0),
        ),
        (
            "recip",
            |x| D2::from_f64(1.0) / x,
            |x| 1.0 / x,
            |x| -1.0 / (x * x),
            |x| 2.0 / (x * x * x),
            (0.2, 4.0),
        ),
    ]
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

/// Worst relative error over 20 points, in units of ε, for value, f', f''.
fn worst() -> Vec<(&'static str, f64, f64, f64)> {
    cases()
        .into_iter()
        .map(|(name, f, f0, f1, f2, (lo, hi))| {
            let (mut w0, mut w1, mut w2) = (0.0f64, 0.0f64, 0.0f64);
            for k in 0..20 {
                let x = lo + (hi - lo) * (k as f64 + 0.5) / 20.0;
                let y = f(D2::variable(x, 0));
                w0 = w0.max(rel(y.re, f0(x)));
                w1 = w1.max(rel(y.d[0], f1(x)));
                w2 = w2.max(rel(y.second(0, 0), f2(x)));
            }
            (
                name,
                w0 / f64::EPSILON,
                w1 / f64::EPSILON,
                w2 / f64::EPSILON,
            )
        })
        .collect()
}

#[test]
#[ignore]
fn measure() {
    for (name, w0, w1, w2) in worst() {
        println!("{name:8} value {w0:6.2} ε   f' {w1:6.2} ε   f'' {w2:6.2} ε");
    }
}

#[test]
fn derivatives_match_analytic_formulas() {
    // Threshold set from the measurement (see `measure`); rationale in the comment
    // of THRESHOLD_EPS.
    for (name, w0, w1, w2) in worst() {
        assert!(
            w0 <= THRESHOLD_EPS && w1 <= THRESHOLD_EPS && w2 <= THRESHOLD_EPS,
            "{name}: {w0} {w1} {w2}"
        );
    }
}

/// Dual numbers have no truncation error, only rounding: each derivative is a
/// handful of floating-point operations on the value. Measured worst case over 20
/// points per function (2026-10-14): tan 2.74 ε (from 1 + tan²), tanh 0.86 ε after
/// computing sech² as 1/cosh² (1 − tanh² had cancelled to 48 ε), all others ≤ 1.5 ε.
/// The design value in 04 (TC-U09) is 4 ε.
const THRESHOLD_EPS: f64 = 4.0;

#[test]
fn mixed_partials_of_a_two_variable_function() {
    // f(x, y) = x² sin y: ∂x = 2x sin y, ∂y = x² cos y, ∂x∂y = 2x cos y, ∂y∂y = −x² sin y.
    let (x0, y0) = (1.3, 0.7);
    let x = D2::variable(x0, 0);
    let y = D2::variable(y0, 1);
    let f = x * x * y.sin();
    let tol = 8.0 * f64::EPSILON;
    assert!(rel(f.d[0], 2.0 * x0 * y0.sin()) <= tol);
    assert!(rel(f.d[1], x0 * x0 * y0.cos()) <= tol);
    assert!(rel(f.second(0, 1), 2.0 * x0 * y0.cos()) <= tol);
    assert!(rel(f.second(1, 1), -x0 * x0 * y0.sin()) <= tol);
    let g = Dual::<Plane>::variable(x0, 0) * Dual::<Plane>::variable(y0, 1);
    assert_eq!(g.d, [y0, x0]);
}

#[test]
fn atan2_derivatives_in_every_quadrant() {
    for &(x0, y0) in &[
        (1.0, 0.5),
        (-1.0, 0.5),
        (-1.0, -0.5),
        (1.0, -0.5),
        (0.1, 2.0),
        (-0.1, -2.0),
    ] {
        let x = D2::variable(x0, 0);
        let y = D2::variable(y0, 1);
        let t = y.atan2(x);
        let r2: f64 = x0 * x0 + y0 * y0;
        assert!(rel(t.re, f64::atan2(y0, x0)) <= 2.0 * f64::EPSILON);
        assert!(rel(t.d[0], -y0 / r2) <= 8.0 * f64::EPSILON);
        assert!(rel(t.d[1], x0 / r2) <= 8.0 * f64::EPSILON);
        assert!(rel(t.second(0, 0), 2.0 * x0 * y0 / (r2 * r2)) <= 16.0 * f64::EPSILON);
    }
}

/// Every component of a `Dual2` of the plane (value, 2 first and 3 second derivatives).
fn parts2(v: D2) -> [f64; 6] {
    [
        v.re(),
        v.d[0],
        v.d[1],
        v.second(0, 0),
        v.second(0, 1),
        v.second(1, 1),
    ]
}

#[test]
fn powers_zero_and_one_are_regular_at_zero() {
    // TC-U09 (Codex review 2026-10-14, finding 2). x^0 = 1 and x^1 = x are smooth,
    // but the general chain-rule formulas evaluate 0 · 0^-1 at x = 0 and gave NaN.
    // The input is x(u, v) = u at u = 0 with ∂_v x = 0 too, so every component counts.
    let x2 = D2::variable(0.0, 0);
    let one = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let same = parts2(x2);
    assert_eq!(parts2(x2.powi(0)), one);
    assert_eq!(parts2(x2.powi(1)), same);
    assert_eq!(parts2(x2.powf(0.0)), one);
    assert_eq!(parts2(x2.powf(1.0)), same);
    let x1 = Dual::<Plane>::variable(0.0, 0);
    for (got, want) in [
        (x1.powi(0), [1.0, 0.0, 0.0]),
        (x1.powi(1), [0.0, 1.0, 0.0]),
        (x1.powf(0.0), [1.0, 0.0, 0.0]),
        (x1.powf(1.0), [0.0, 1.0, 0.0]),
    ] {
        assert_eq!([got.re(), got.d[0], got.d[1]], want);
    }
    // x² at 0 was already regular and stays so; a negative power at 0 is singular.
    assert_eq!(parts2(x2.powi(2)), [0.0, 0.0, 0.0, 2.0, 0.0, 0.0]);
    assert!(!x2.powi(-1).re().is_finite());
}

#[test]
fn powi_accepts_every_i32_exponent() {
    // TC-U09: k − 1 and k − 2 overflowed i32 for k = i32::MIN (a panic in debug builds).
    // At x = 1 every power is 1, the derivatives are k and k (k − 1).
    let k = i32::MIN;
    let v = D2::variable(1.0, 0).powi(k);
    let kf = k as f64;
    assert_eq!(parts2(v), [1.0, kf, 0.0, kf * (kf - 1.0), 0.0, 0.0]);
    let v = D2::variable(1.0, 0).powi(i32::MIN + 1);
    assert_eq!(v.re(), 1.0);
}
