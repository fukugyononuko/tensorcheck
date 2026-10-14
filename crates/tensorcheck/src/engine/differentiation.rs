//! Derivatives (CMP-E08, DD-09, DD-15).
//!
//! Automatic differentiation: a point is turned into dual numbers (one variable per
//! coordinate), a function written generically over `N: Scalar` is evaluated once,
//! and the values and derivatives are read back. Finite differences for functions
//! that are not generic come in a later phase.

use crate::contracts::{Chart, Point};
use crate::utility::autodiff::{Dual, Dual2};
use crate::utility::extent::Extent;

/// The space of chart `C`, as an extent.
pub type SpaceOf<C> = <C as Chart>::Space;

/// The point `p` with each coordinate an independent first-order dual variable.
pub fn seed1<C: Chart>(p: &Point<C>) -> Point<C, Dual<SpaceOf<C>>> {
    Point::from_base(<SpaceOf<C> as Extent>::from_fn(|i| {
        Dual::variable(p.get(i), i)
    }))
}

/// The point `p` with each coordinate an independent second-order dual variable.
pub fn seed2<C: Chart>(p: &Point<C>) -> Point<C, Dual2<SpaceOf<C>>> {
    Point::from_base(<SpaceOf<C> as Extent>::from_fn(|i| {
        Dual2::variable(p.get(i), i)
    }))
}

/// Values and first derivatives of a list of first-order duals:
/// `values[k]`, `first[c * len + k] = ∂_c v_k`.
pub fn unpack1<E: Extent>(v: &[Dual<E>]) -> (Vec<f64>, Vec<f64>) {
    let len = v.len();
    let n = E::N;
    let values = v.iter().map(|x| x.re).collect();
    let mut first = vec![0.0; n * len];
    for (k, x) in v.iter().enumerate() {
        for c in 0..n {
            first[c * len + k] = x.d.as_ref()[c];
        }
    }
    (values, first)
}

/// Values, first and second derivatives of a list of second-order duals:
/// `first[c * len + k] = ∂_c v_k`, `second[(c * n + d) * len + k] = ∂_c∂_d v_k`.
pub fn unpack2<E: Extent>(v: &[Dual2<E>]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let len = v.len();
    let n = E::N;
    let values = v.iter().map(|x| x.re).collect();
    let mut first = vec![0.0; n * len];
    let mut second = vec![0.0; n * n * len];
    for (k, x) in v.iter().enumerate() {
        for c in 0..n {
            first[c * len + k] = x.d.as_ref()[c];
            for d in 0..n {
                second[(c * n + d) * len + k] = x.second(c, d);
            }
        }
    }
    (values, first, second)
}

// ---------------------------------------------------------------------------
// Finite differences (DD-15), for functions that are not generic over `N: Scalar`.

/// Central difference `(f(x+h) − f(x−h)) / 2h` with `h = ε^(1/3) max(|x|, scale)`,
/// which balances rounding (`ε/h`) against truncation (`h²`).
pub fn central(f: &mut dyn FnMut(f64) -> f64, x: f64, scale: f64) -> f64 {
    let h = f64::EPSILON.cbrt() * x.abs().max(scale);
    let (xp, xm) = (x + h, x - h);
    (f(xp) - f(xm)) / (xp - xm)
}

/// Ridders' method: central differences with a shrinking step, extrapolated in a
/// Neville tableau; stops when the error starts to grow. Returns `(derivative,
/// error estimate)`. `h0` is the initial step.
pub fn ridders(f: &mut dyn FnMut(f64) -> f64, x: f64, h0: f64) -> (f64, f64) {
    const CON: f64 = 1.4;
    const CON2: f64 = CON * CON;
    const NTAB: usize = 10;
    const SAFE: f64 = 2.0;
    let mut a = [[0.0f64; NTAB]; NTAB];
    let mut h = h0;
    // Largest |f| among the points used: the rounding of a difference scales with it.
    let mut fmax: f64 = 0.0;
    let diff = |f: &mut dyn FnMut(f64) -> f64, h: f64, fmax: &mut f64| {
        let (fp, fm) = (f(x + h), f(x - h));
        *fmax = fmax.max(fp.abs()).max(fm.abs());
        (fp - fm) / (2.0 * h)
    };
    a[0][0] = diff(f, h, &mut fmax);
    let mut err = f64::MAX;
    let mut ans = a[0][0];
    let mut h_ans = h;
    for i in 1..NTAB {
        h /= CON;
        a[0][i] = diff(f, h, &mut fmax);
        let mut fac = CON2;
        for j in 1..=i {
            a[j][i] = (a[j - 1][i] * fac - a[j - 1][i - 1]) / (fac - 1.0);
            fac *= CON2;
            let errt = (a[j][i] - a[j - 1][i])
                .abs()
                .max((a[j][i] - a[j - 1][i - 1]).abs());
            if errt <= err {
                err = errt;
                ans = a[j][i];
                h_ans = h;
            }
        }
        if (a[i][i] - a[i - 1][i - 1]).abs() >= SAFE * err {
            break;
        }
    }
    // The tableau sees truncation error only; add the rounding floor so the
    // estimate does not fall below the actual error (TC-U26). A central difference
    // with step h carries about 2 ε |f| / h of rounding, and each extrapolation level
    // multiplies the coefficient sum by (fac + 1)/(fac − 1) for fac = 1.96, 3.84,
    // 7.53, …, a product of about 9: hence 20 ε |f| / h.
    let floor = 20.0 * f64::EPSILON * fmax.max(f64::MIN_POSITIVE) / h_ans;
    // The tableau difference is a heuristic for the truncation error, not a bound:
    // over 4 functions × 200 points it fell to 0.54 of the actual error. A safety
    // factor of 4 keeps the estimate above the actual error on that set (TC-U26).
    (ans, (4.0 * err).max(floor))
}

/// The Jacobian `∂f_i/∂y_j` of `f: R^n → R^m` by central differences, row-major
/// (`m × n`). `scale[j]` is a typical magnitude of `y_j` (used when `y_j ≈ 0`).
pub fn jacobian_fd(
    f: &mut dyn FnMut(&[f64], &mut [f64]),
    y: &[f64],
    m: usize,
    scale: &[f64],
) -> Vec<f64> {
    let n = y.len();
    let mut jac = vec![0.0; m * n];
    let mut yp = y.to_vec();
    let (mut fp, mut fm) = (vec![0.0; m], vec![0.0; m]);
    for j in 0..n {
        let h = f64::EPSILON.cbrt() * y[j].abs().max(scale[j]);
        yp[j] = y[j] + h;
        f(&yp, &mut fp);
        let hp = yp[j] - y[j];
        yp[j] = y[j] - h;
        f(&yp, &mut fm);
        let hm = y[j] - yp[j];
        yp[j] = y[j];
        for i in 0..m {
            jac[i * n + j] = (fp[i] - fm[i]) / (hp + hm);
        }
    }
    jac
}
