//! The solver manager (CMP-M02): runs the numerical procedures and converts between
//! dimensioned states and the engines' plain arrays (2b §CMP-M02, CC-07, CC-08).

use core::cell::Cell;
use core::ops::Range;

use tensorcheck::{Quantity, TimesDim, VectorSpace};

use crate::engine::error_control::{judge, max_norm, pick_interval};
use crate::engine::quadrature::{gauss_kronrod15, gauss_legendre, romberg, simpson, trapezoid};
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
