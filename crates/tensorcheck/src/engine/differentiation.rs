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
