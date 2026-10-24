//! The solver manager (CMP-M02): runs the numerical procedures and converts between
//! dimensioned states and the engines' plain arrays (2b §CMP-M02, CC-07, CC-08).

use core::cell::Cell;
use core::ops::Range;

use tensorcheck::{Quantity, TimesDim, VectorSpace};

use crate::engine::error_control::{
    judge, max_norm, next_step, pick_interval, rms_norm, valid_tolerance,
};
use crate::engine::ode_step::{
    MatrixError, RosenbrockWork, Work, dopri5_dense, dopri5_step, hermite, illinois, rk4,
    rodas_step,
};
use crate::engine::quadrature::{gauss_kronrod15, gauss_legendre, romberg, simpson, trapezoid};
use crate::ode::{Solution, SolveError, Stats};
use crate::quad::{Integral, QuadError};
use tensorcheck::internal::engine::{central, jacobian_fd, ridders};

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
    // Input checks, before the integrand is evaluated (DD-13, DD-17).
    if !valid_tolerance(rtol, &atol) {
        return Err(QuadError::InvalidTolerance);
    }
    if let Some(n) = invalid_count(kind) {
        return Err(QuadError::InvalidCount { n });
    }
    if !a.is_finite() || !b.is_finite() {
        return Err(QuadError::InvalidRange);
    }
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

/// The count of a fixed rule if it is unsupported: zero, or so large that the
/// refined rule's count (`2n` intervals, Simpson's even `n`, `n + 1` points)
/// overflows (DD-14).
fn invalid_count(kind: QuadKind) -> Option<usize> {
    let ok = match kind {
        QuadKind::Trapezoid(n) => n >= 1 && n.checked_mul(2).is_some(),
        QuadKind::Simpson(n) => {
            n >= 1
                && n.checked_add(n % 2)
                    .and_then(|e| e.checked_mul(2))
                    .is_some()
        }
        QuadKind::GaussLegendre(n) => n >= 1 && n.checked_add(1).is_some(),
        QuadKind::Romberg | QuadKind::GaussKronrod => true,
    };
    match kind {
        QuadKind::Trapezoid(n) | QuadKind::Simpson(n) | QuadKind::GaussLegendre(n) if !ok => {
            Some(n)
        }
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// ODEs (CC-07)

/// The ODE method, as seen by the manager.
#[doc(hidden)]
#[derive(Clone, Copy, Debug)]
pub enum OdeKind {
    Rk4(f64),
    Dopri5,
    Rodas,
}

/// The event condition of a solve, on plain arrays (DD-12, FR-31). `()` means no
/// events: `ACTIVE` is false, so the checks compile away.
pub(crate) trait Watch {
    const ACTIVE: bool;
    /// `g(t, y)` in base units.
    fn value(&mut self, t: f64, y: &[f64]) -> f64;
    /// Whether a sign change starting from `g_before` counts (the direction).
    fn wants(&self, g_before: f64) -> bool;
    /// Whether the solve stops at the first event.
    fn terminal(&self) -> bool;
    /// Records an event.
    fn record(&mut self, t: f64, y: &[f64]);
}
impl Watch for () {
    const ACTIVE: bool = false;
    fn value(&mut self, _: f64, _: &[f64]) -> f64 {
        0.0
    }
    fn wants(&self, _: f64) -> bool {
        false
    }
    fn terminal(&self) -> bool {
        false
    }
    fn record(&mut self, _: f64, _: &[f64]) {}
}

/// How the state inside a step is interpolated for events (DD-12).
enum Dense<'a> {
    /// Cubic Hermite from `f(t0, y0)` and `f(t1, y1)`; `None` evaluates `f(t1, y1)`
    /// only when the sign changed (RK4).
    Hermite(&'a [f64], Option<&'a [f64]>),
    /// The continuous extension of Dormand–Prince from the stages of the step.
    Dopri5(&'a [Vec<f64>]),
}

/// Checks one accepted step from `(t0, y0)` to `(t1, y1)` for an event (DD-12).
/// Updates `g_before` to `g(t1, y1)`. Returns the point where a terminal event stops
/// the solve, or `Err(t)` when `g` was not finite at `t`.
#[allow(clippy::too_many_arguments)]
fn watch<W: Watch + ?Sized, G: FnMut(f64, &[f64], &mut [f64]) + ?Sized>(
    events: &mut W,
    g_before: &mut f64,
    rhs: &mut G,
    t0: f64,
    y0: &[f64],
    t1: f64,
    y1: &[f64],
    dense: Dense,
) -> Result<Option<(f64, Vec<f64>)>, f64> {
    let (g0, g1) = (*g_before, events.value(t1, y1));
    if !g1.is_finite() {
        return Err(t1);
    }
    *g_before = g1;
    // Compare signs, not the product: g0 g1 underflows to −0 for tiny values (PR #41).
    let changed = (g0 < 0.0 && g1 >= 0.0) || (g0 > 0.0 && g1 <= 0.0);
    if !changed || !events.wants(g0) {
        return Ok(None);
    }
    let mut computed = Vec::new();
    let h = t1 - t0;
    // What the interpolant is built from must be finite, like the state (PR #41).
    let finite = |v: &[f64]| v.iter().all(|x| x.is_finite());
    let interpolate = match dense {
        Dense::Hermite(f0, f1) => {
            let f1 = match f1 {
                Some(f) => f,
                None => {
                    computed.resize(y1.len(), 0.0);
                    rhs(t1, y1, &mut computed);
                    &computed
                }
            };
            if !finite(f1) {
                return Err(t1);
            }
            Box::new(move |theta: f64, out: &mut [f64]| hermite(theta, h, y0, f0, y1, f1, out))
                as Box<dyn Fn(f64, &mut [f64]) + '_>
        }
        Dense::Dopri5(k) => {
            if !k.iter().all(|ki| finite(ki)) {
                return Err(t1);
            }
            Box::new(move |theta: f64, out: &mut [f64]| dopri5_dense(theta, h, y0, y1, k, out))
        }
    };
    let mut y = vec![0.0; y1.len()];
    let mut bad = None;
    // Every interpolated state is checked before g sees it or it is returned: finite
    // inputs can still overflow in the interpolant (PR #41).
    let mut g = |theta: f64| {
        interpolate(theta, &mut y);
        let v = if finite(&y) {
            events.value(t0 + theta * h, &y)
        } else {
            f64::NAN
        };
        if !v.is_finite() && bad.is_none() {
            bad = Some(t0 + theta * h);
        }
        v
    };
    let theta = illinois(&mut g, g0, g1);
    if let Some(t) = bad {
        return Err(t);
    }
    let t = if theta == 1.0 { t1 } else { t0 + theta * h };
    interpolate(theta, &mut y);
    if !finite(&y) {
        return Err(t);
    }
    events.record(t, &y);
    Ok(events.terminal().then_some((t, y)))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_ivp<T, Y, R, F, J, A, O, W>(
    mut f: F,
    span: Range<Quantity<T>>,
    y0: Y,
    kind: OdeKind,
    max_steps: usize,
    given_jacobian: &mut J,
    rtol: f64,
    atol: &A,
    on_step: &mut O,
    events: &mut W,
) -> Result<Solution<T, Y>, SolveError<T>>
where
    F: FnMut(Quantity<T>, &Y) -> R,
    Y: VectorSpace,
    R: VectorSpace,
    J: FnMut(f64, &[f64], &mut [f64]) -> bool + ?Sized,
    A: VectorSpace,
    O: FnMut(Quantity<T>, &Y) + ?Sized,
    W: Watch + ?Sized,
{
    let n = Y::LEN;
    let (t0, t1) = (span.start.base_value(), span.end.base_value());
    let dir = if t1 >= t0 { 1.0 } else { -1.0 };
    let atol: Vec<f64> = (0..n).map(|i| atol.component(i)).collect();
    // Input checks, before the right-hand side is evaluated (DD-12, DD-13, DD-17).
    if !t0.is_finite() || !t1.is_finite() {
        return Err(SolveError::InvalidSpan);
    }
    if !valid_tolerance(rtol, &atol) {
        return Err(SolveError::InvalidTolerance);
    }
    if let OdeKind::Rk4(h) = kind {
        // A positive magnitude: the interval gives the direction. Zero never advances
        // and a negative step moves away from the end, both until MaxSteps.
        if !(h.is_finite() && h > 0.0) {
            return Err(SolveError::InvalidStep { h });
        }
    }
    if (0..n).any(|i| !y0.component(i).is_finite()) {
        // Even for an empty interval: a non-finite state is never a solution.
        return Err(SolveError::NonFinite { t: span.start });
    }
    if max_steps == 0 && dir * (t1 - t0) > 0.0 {
        // No step allowed: stop before the adaptive methods evaluate f for the
        // initial step size.
        return Err(SolveError::MaxSteps { t: span.start });
    }
    let mut g_before = 0.0;
    if W::ACTIVE && dir * (t1 - t0) > 0.0 {
        let y: Vec<f64> = (0..n).map(|i| y0.component(i)).collect();
        g_before = events.value(t0, &y);
        if !g_before.is_finite() {
            return Err(SolveError::NonFinite { t: span.start });
        }
    }
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
    // The end of the solve at a terminal event: on_step sees the event point last.
    macro_rules! stop_at {
        ($watched:expr) => {
            match $watched {
                Err(tb) => return Err(SolveError::NonFinite { t: at(tb) }),
                Ok(Some((te, ye))) => {
                    stats.steps += 1;
                    let ye = Y::from_components(|i| ye[i]);
                    on_step(at(te), &ye);
                    stats.evaluations = evals.get();
                    return Ok(Solution {
                        t: at(te),
                        y: ye,
                        stats,
                    });
                }
                Ok(None) => {}
            }
        };
    }
    match kind {
        OdeKind::Rk4(h_abs) => {
            let mut w = Work::new(4, n);
            while dir * (t1 - t) > 0.0 {
                if stats.steps >= max_steps {
                    return Err(SolveError::MaxSteps { t: at(t) });
                }
                let h = dir * h_abs.min(dir * (t1 - t));
                rk4(&mut rhs, t, &y, h, &mut w, &mut y_new);
                if !finite(&y_new) {
                    return Err(SolveError::NonFinite { t: at(t) });
                }
                if W::ACTIVE {
                    let w0 = &w.k[0];
                    stop_at!(watch(
                        events,
                        &mut g_before,
                        &mut rhs,
                        t,
                        &y,
                        t + h,
                        &y_new,
                        Dense::Hermite(w0, None)
                    ));
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
                if stats.steps + stats.rejected >= max_steps {
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
                    if W::ACTIVE {
                        let k = &w.k;
                        stop_at!(watch(
                            events,
                            &mut g_before,
                            &mut rhs,
                            t,
                            &y,
                            t + h,
                            &y_new,
                            Dense::Dopri5(k)
                        ));
                    }
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
        OdeKind::Rodas => {
            let mut w = RosenbrockWork::new(n);
            let mut err = vec![0.0; n];
            let mut f0 = vec![0.0; n];
            let mut ft = vec![0.0; n];
            let mut f1 = vec![0.0; n];
            let mut jac = vec![0.0; n * n];
            let scale: Vec<f64> = atol.iter().map(|a| a.max(f64::MIN_POSITIVE)).collect();
            let t_scale = (t1 - t0).abs().max(f64::MIN_POSITIVE);
            rhs(t, &y, &mut f0);
            let mut h = initial_step(&mut rhs, t, &y, &f0, dir, 3, &atol, rtol);
            let mut rejected_last = false;
            // Singular matrices E so far (never reset, as `NSING` in `rodas.f`).
            let mut singular = 0;
            while dir * (t1 - t) > 0.0 {
                if stats.steps + stats.rejected >= max_steps {
                    return Err(SolveError::MaxSteps { t: at(t) });
                }
                if h.abs() < 16.0 * f64::EPSILON * t.abs().max(f64::MIN_POSITIVE) {
                    return Err(SolveError::StepTooSmall { t: at(t), h });
                }
                if dir * (t + h - t1) > 0.0 {
                    h = t1 - t;
                }
                // Jacobian ∂f/∂y at (t, y): the user's (FR-30), else central differences
                // (E08). The user's must be finite (DD-15).
                if given_jacobian(t, &y, &mut jac) {
                    if !finite(&jac) {
                        return Err(SolveError::NonFinite { t: at(t) });
                    }
                } else {
                    let tt = t;
                    jac = jacobian_fd(
                        &mut |yy: &[f64], out: &mut [f64]| rhs(tt, yy, out),
                        &y,
                        n,
                        &scale,
                    );
                }
                // Time derivative ∂f/∂t at (t, y) by central differences (E08).
                let yy = y.clone();
                let dfdt = jacobian_fd(
                    &mut |ts: &[f64], out: &mut [f64]| rhs(ts[0], &yy, out),
                    &[t],
                    n,
                    &[t_scale],
                );
                ft.copy_from_slice(&dfdt);
                match rodas_step(
                    &mut rhs, t, &y, h, &f0, &jac, &ft, &mut w, &mut y_new, &mut err,
                ) {
                    Ok(()) => {}
                    // E = I/(hγ) − J depends on h: a smaller step can make it regular.
                    // As in `ROSCOR` of `rodas.f` (label 80, `NSING`), halve it, and
                    // give up at the fifth singular E of the solve (PR #59 review).
                    Err(MatrixError::Singular) => {
                        singular += 1;
                        if singular >= 5 {
                            return Err(SolveError::Singular { t: at(t) });
                        }
                        stats.rejected += 1;
                        h *= 0.5;
                        rejected_last = true;
                        continue;
                    }
                    // A non-finite J does not depend on h, and a smaller step only
                    // enlarges 1/(hγ): not retried.
                    Err(MatrixError::NonFinite) => {
                        return Err(SolveError::NonFinite { t: at(t) });
                    }
                }
                if !finite(&y_new) || !finite(&err) {
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
                    // f at the new point: for the next step, and for the event check.
                    rhs(t + h, &y_new, &mut f1);
                    if W::ACTIVE {
                        stop_at!(watch(
                            events,
                            &mut g_before,
                            &mut rhs,
                            t,
                            &y,
                            t + h,
                            &y_new,
                            Dense::Hermite(&f0, Some(&f1))
                        ));
                    }
                    t += h;
                    core::mem::swap(&mut y, &mut y_new);
                    core::mem::swap(&mut f0, &mut f1);
                    stats.steps += 1;
                    on_step(at(t), &Y::from_components(|i| y[i]));
                    h = next_step(h, e, 3, rejected_last);
                    rejected_last = false;
                } else {
                    stats.rejected += 1;
                    h = next_step(h, e, 3, true);
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
fn initial_step<G: FnMut(f64, &[f64], &mut [f64]) + ?Sized>(
    rhs: &mut G,
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

// ---------------------------------------------------------------------------
// Numerical differentiation (2b §CMP-M02 `derivative`, `jacobian`)

/// `df/dx` per component, in base units, and the largest error estimate.
pub(crate) fn derivative<X, Y, F>(mut f: F, x: Quantity<X>, richardson: bool) -> (Vec<f64>, f64)
where
    F: FnMut(Quantity<X>) -> Y,
    Y: VectorSpace,
{
    let x0 = x.base_value();
    let mut values = vec![0.0; Y::LEN];
    let mut worst: f64 = 0.0;
    for (k, v) in values.iter_mut().enumerate() {
        let mut g = |xx: f64| f(Quantity::from_base_components(xx)).component(k);
        if richardson {
            let (d, e) = ridders(&mut g, x0, 0.1 * x0.abs().max(1.0));
            *v = d;
            worst = worst.max(e);
        } else {
            let d = central(&mut g, x0, 1.0);
            let h = f64::EPSILON.cbrt() * x0.abs().max(1.0) * 0.5;
            let d2 = (g(x0 + h) - g(x0 - h)) / (2.0 * h);
            *v = d;
            worst = worst.max((d - d2).abs());
        }
    }
    (values, worst)
}

/// The Jacobian of `f` at `y` by central differences, row-major, base units.
pub(crate) fn jacobian<Y, R, F>(mut f: F, y: &Y) -> Vec<f64>
where
    F: FnMut(&Y) -> R,
    Y: VectorSpace,
    R: VectorSpace,
{
    let y0: Vec<f64> = (0..Y::LEN).map(|i| y.component(i)).collect();
    let scale = vec![1.0; Y::LEN];
    jacobian_fd(
        &mut |yy: &[f64], out: &mut [f64]| {
            let r = f(&Y::from_components(|i| yy[i]));
            for (k, o) in out.iter_mut().enumerate() {
                *o = r.component(k);
            }
        },
        &y0,
        R::LEN,
        &scale,
    )
}
