//! The solver manager (CMP-M02): runs the numerical procedures and converts between
//! dimensioned states and the engines' plain arrays (2b §CMP-M02, CC-07, CC-08).

use core::cell::Cell;
use core::ops::Range;

use tensorcheck::{Quantity, TimesDim, VectorSpace};

use crate::engine::error_control::{judge, max_norm, next_step, pick_interval, rms_norm};
use crate::engine::ode_step::{Rhs, Work, dopri5_step, rk4};
use crate::engine::quadrature::{gauss_kronrod15, gauss_legendre, romberg, simpson, trapezoid};
use crate::ode::{Solution, SolveError, Stats};
use crate::quad::{Integral, QuadError};

/// The quadrature method, as seen by the manager.
#[doc(hidden)]
#[derive(Clone, Copy, Debug)]
pub enum QuadKind {
    Trapezoid(usize),
    Simpson(usize),
    Romberg,
    GaussLegendre(usize),
    GaussKronrod,
}

/// Largest number of subintervals for adaptive Gauss–Kronrod (DD-13).
pub const MAX_INTERVALS: usize = 1000;

pub(crate) fn integrate<X, Y, F, A>(
    mut f: F,
    range: Range<Quantity<X>>,
    kind: QuadKind,
    rtol: f64,
    atol: &A,
) -> Result<Integral<<Y as TimesDim<X>>::Out>, QuadError>
where
    F: FnMut(Quantity<X>) -> Y,
    Y: VectorSpace + TimesDim<X>,
    A: VectorSpace,
{
    let m = Y::LEN;
    let (a, b) = (range.start.base_value(), range.end.base_value());
    let atol: Vec<f64> = (0..m).map(|i| atol.component(i)).collect();
    let evals = Cell::new(0usize);
    let finite = Cell::new(true);
    let mut g = |x: f64, out: &mut [f64]| {
        evals.set(evals.get() + 1);
        let y = f(Quantity::from_base_components(x));
        for (k, o) in out.iter_mut().enumerate() {
            *o = y.component(k);
            if !o.is_finite() {
                finite.set(false);
            }
        }
    };
    let (value, err) = match kind {
        QuadKind::Trapezoid(n) => {
            let coarse = trapezoid(&mut g, a, b, n, m);
            let fine = trapezoid(&mut g, a, b, 2 * n, m);
            let e: Vec<f64> = fine
                .iter()
                .zip(&coarse)
                .map(|(x, y)| (x - y).abs())
                .collect();
            (fine, e)
        }
        QuadKind::Simpson(n) => {
            let coarse = simpson(&mut g, a, b, n, m);
            let fine = simpson(&mut g, a, b, 2 * n, m);
            let e: Vec<f64> = fine
                .iter()
                .zip(&coarse)
                .map(|(x, y)| (x - y).abs())
                .collect();
            (fine, e)
        }
        QuadKind::GaussLegendre(n) => {
            let lo = gauss_legendre(&mut g, a, b, n, m);
            let hi = gauss_legendre(&mut g, a, b, n + 1, m);
            let e: Vec<f64> = hi.iter().zip(&lo).map(|(x, y)| (x - y).abs()).collect();
            (lo, e)
        }
        QuadKind::Romberg => {
            let mut done = |best: &[f64], e: &[f64]| judge(max_norm(e, best, &atol, rtol));
            let (v, e, _levels) = romberg(&mut g, a, b, m, 20, &mut done);
            (v, e)
        }
        QuadKind::GaussKronrod => {
            // Adaptive bisection of the interval with the largest scaled error.
            let mut intervals: Vec<(f64, f64, Vec<f64>, Vec<f64>)> = Vec::new();
            let (r, e) = gauss_kronrod15(&mut g, a, b, m);
            intervals.push((a, b, r, e));
            loop {
                let total: Vec<f64> = (0..m)
                    .map(|k| intervals.iter().map(|iv| iv.2[k]).sum())
                    .collect();
                let total_err: Vec<f64> = (0..m)
                    .map(|k| intervals.iter().map(|iv| iv.3[k]).sum())
                    .collect();
                if judge(max_norm(&total_err, &total, &atol, rtol))
                    || intervals.len() >= MAX_INTERVALS
                    || !finite.get()
                {
                    break (total, total_err);
                }
                let scaled: Vec<f64> = intervals
                    .iter()
                    .map(|iv| max_norm(&iv.3, &total, &atol, rtol))
                    .collect();
                let (lo, hi, _, _) = intervals.swap_remove(pick_interval(&scaled));
                let mid = 0.5 * (lo + hi);
                let (r1, e1) = gauss_kronrod15(&mut g, lo, mid, m);
                let (r2, e2) = gauss_kronrod15(&mut g, mid, hi, m);
                intervals.push((lo, mid, r1, e1));
                intervals.push((mid, hi, r2, e2));
            }
        }
    };
    if !finite.get() || value.iter().any(|v| !v.is_finite()) {
        return Err(QuadError::NonFinite);
    }
    let scaled = max_norm(&err, &value, &atol, rtol);
    if !judge(scaled) {
        return Err(QuadError::ToleranceNotMet {
            scaled_error: scaled,
        });
    }
    Ok(Integral {
        value: <Y as TimesDim<X>>::Out::from_components(|k| value[k]),
        error_estimate: err.iter().fold(0.0, |acc, e| acc.max(*e)),
        evaluations: evals.get(),
    })
}

// ---------------------------------------------------------------------------
// ODEs (CC-07)

/// The ODE method, as seen by the manager.
#[doc(hidden)]
#[derive(Clone, Copy, Debug)]
pub enum OdeKind {
    Rk4(f64),
    Dopri5,
}

pub(crate) fn solve_ivp<T, Y, R, F, A>(
    mut f: F,
    span: Range<Quantity<T>>,
    y0: Y,
    kind: OdeKind,
    rtol: f64,
    atol: &A,
    on_step: &mut dyn FnMut(Quantity<T>, &Y),
) -> Result<Solution<T, Y>, SolveError<T>>
where
    F: FnMut(Quantity<T>, &Y) -> R,
    Y: VectorSpace,
    R: VectorSpace,
    A: VectorSpace,
{
    let n = Y::LEN;
    let (t0, t1) = (span.start.base_value(), span.end.base_value());
    let dir = if t1 >= t0 { 1.0 } else { -1.0 };
    let atol: Vec<f64> = (0..n).map(|i| atol.component(i)).collect();
    let mut stats = Stats::default();
    let evals = Cell::new(0usize);
    let mut rhs = |t: f64, y: &[f64], out: &mut [f64]| {
        evals.set(evals.get() + 1);
        let r = f(
            Quantity::from_base_components(t),
            &Y::from_components(|i| y[i]),
        );
        for (k, o) in out.iter_mut().enumerate() {
            *o = r.component(k);
        }
    };
    let at = |t: f64| Quantity::from_base_components(t);
    let mut y: Vec<f64> = (0..n).map(|i| y0.component(i)).collect();
    let mut y_new = vec![0.0; n];
    let mut t = t0;
    let finite = |v: &[f64]| v.iter().all(|x| x.is_finite());
    match kind {
        OdeKind::Rk4(h_abs) => {
            let mut w = Work::new(4, n);
            while dir * (t1 - t) > 0.0 {
                if stats.steps >= crate::ode::MAX_STEPS {
                    return Err(SolveError::MaxSteps { t: at(t) });
                }
                let h = dir * h_abs.min(dir * (t1 - t));
                rk4(&mut rhs, t, &y, h, &mut w, &mut y_new);
                if !finite(&y_new) {
                    return Err(SolveError::NonFinite { t: at(t) });
                }
                t += h;
                core::mem::swap(&mut y, &mut y_new);
                stats.steps += 1;
                on_step(at(t), &Y::from_components(|i| y[i]));
            }
        }
        OdeKind::Dopri5 => {
            let mut w = Work::new(7, n);
            let mut err = vec![0.0; n];
            rhs(t, &y, &mut w.k[0]);
            let mut h = initial_step(&mut rhs, t, &y, &w.k[0], dir, 4, &atol, rtol);
            let mut rejected_last = false;
            while dir * (t1 - t) > 0.0 {
                if stats.steps + stats.rejected >= crate::ode::MAX_STEPS {
                    return Err(SolveError::MaxSteps { t: at(t) });
                }
                if h.abs() < 16.0 * f64::EPSILON * t.abs().max(f64::MIN_POSITIVE) {
                    return Err(SolveError::StepTooSmall { t: at(t), h });
                }
                if dir * (t + h - t1) > 0.0 {
                    h = t1 - t;
                }
                dopri5_step(&mut rhs, t, &y, h, &mut w, &mut y_new, &mut err);
                if !finite(&y_new) || !finite(&err) {
                    // Treat as a failed step first; give up when the step cannot shrink.
                    stats.rejected += 1;
                    h *= 0.25;
                    if h.abs() < 16.0 * f64::EPSILON * t.abs().max(f64::MIN_POSITIVE) {
                        return Err(SolveError::NonFinite { t: at(t) });
                    }
                    rejected_last = true;
                    continue;
                }
                let e = rms_norm(&err, &y, &y_new, &atol, rtol);
                if judge(e) {
                    t += h;
                    core::mem::swap(&mut y, &mut y_new);
                    let (k0, k6) = (0, 6);
                    let last = w.k[k6].clone();
                    w.k[k0].copy_from_slice(&last);
                    stats.steps += 1;
                    on_step(at(t), &Y::from_components(|i| y[i]));
                    h = next_step(h, e, 4, rejected_last);
                    rejected_last = false;
                } else {
                    stats.rejected += 1;
                    h = next_step(h, e, 4, true);
                    rejected_last = true;
                }
            }
        }
    }
    stats.evaluations = evals.get();
    Ok(Solution {
        t: at(t),
        y: Y::from_components(|i| y[i]),
        stats,
    })
}

/// Hairer's initial step-size estimate (DD-13): two right-hand-side evaluations.
#[allow(clippy::too_many_arguments)]
fn initial_step(
    rhs: &mut Rhs,
    t: f64,
    y: &[f64],
    f0: &[f64],
    dir: f64,
    q: u32,
    atol: &[f64],
    rtol: f64,
) -> f64 {
    let n = y.len();
    let zeros = vec![0.0; n];
    let d0 = rms_norm(y, y, &zeros, atol, rtol);
    let d1 = rms_norm(f0, y, &zeros, atol, rtol);
    let h0 = if d0 < 1e-5 || d1 < 1e-5 {
        1e-6
    } else {
        0.01 * d0 / d1
    };
    let y1: Vec<f64> = (0..n).map(|i| y[i] + dir * h0 * f0[i]).collect();
    let mut f1 = vec![0.0; n];
    rhs(t + dir * h0, &y1, &mut f1);
    let df: Vec<f64> = (0..n).map(|i| f1[i] - f0[i]).collect();
    let d2 = rms_norm(&df, y, &zeros, atol, rtol) / h0;
    let h1 = if d1.max(d2) <= 1e-15 {
        (h0 * 1e-3).max(1e-6)
    } else {
        (0.01 / d1.max(d2)).powf(1.0 / (q as f64 + 1.0))
    };
    dir * (100.0 * h0).min(h1)
}
