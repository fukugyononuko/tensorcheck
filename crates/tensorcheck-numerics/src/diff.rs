//! Numerical derivatives that keep dimensions (API-15, API-49).
//!
//! For metric tensors and tensor fields written generically over `N: Scalar`, use
//! the exact automatic differentiation of `tensorcheck::geometry` instead. These
//! finite-difference routines are for functions that are not generic.
//!
//! ```
//! use tensorcheck::Quantity;
//! use tensorcheck::units::si::{Length, Time, Velocity, m, s};
//! use tensorcheck_numerics::diff::{Richardson, derivative};
//!
//! let x = |t: Quantity<Time>| -> Quantity<Length> { Quantity::new(t.value_in(s).sin(), m) };
//! let d = derivative(x, 0.5 * s, Richardson);
//! let v: Quantity<Velocity> = d.value;
//! assert!((v.value_in(m / s) - 0.5f64.cos()).abs() < 1e-10);
//! ```

use tensorcheck::{PerDim, Quantity, VectorSpace};

use crate::manager::solver;

/// Central difference with the step balancing rounding and truncation errors.
#[derive(Clone, Copy, Debug)]
pub struct Central;
/// Ridders' method: Richardson extrapolation of central differences, with an
/// error estimate.
#[derive(Clone, Copy, Debug)]
pub struct Richardson;

/// A finite-difference method.
pub trait Method: Copy {
    #[doc(hidden)]
    fn richardson(self) -> bool;
}
impl Method for Central {
    fn richardson(self) -> bool {
        false
    }
}
impl Method for Richardson {
    fn richardson(self) -> bool {
        true
    }
}

/// A derivative and its estimated error.
#[derive(Clone, Copy, Debug)]
pub struct Derivative<V> {
    /// `df/dx`.
    pub value: V,
    /// Estimated absolute error, per component, in base units (largest). For
    /// [`Central`], the difference to the result with half the step.
    pub error_estimate: f64,
}

/// `df/dx` at `x`, with the dimension of `f` divided by that of `x` (API-49).
pub fn derivative<X, Y, F, M>(f: F, x: Quantity<X>, method: M) -> Derivative<<Y as PerDim<X>>::Out>
where
    F: FnMut(Quantity<X>) -> Y,
    Y: VectorSpace + PerDim<X>,
    M: Method,
{
    let (values, err) = solver::derivative(f, x, method.richardson());
    Derivative {
        value: <Y as PerDim<X>>::Out::from_components(|k| values[k]),
        error_estimate: err,
    }
}

/// The Jacobian `∂f_i/∂y_j` of `f` at `y`, by central differences (API-49).
/// Element `(i, j)` has dimension `[f_i]/[y_j]`; the matrix stores base units.
pub fn jacobian<Y, R, F>(f: F, y: &Y) -> JacobianMatrix
where
    F: FnMut(&Y) -> R,
    Y: VectorSpace,
    R: VectorSpace,
{
    let data = solver::jacobian(f, y);
    JacobianMatrix {
        rows: R::LEN,
        cols: Y::LEN,
        data,
    }
}

/// A Jacobian matrix in base units, row-major.
#[derive(Clone, Debug)]
pub struct JacobianMatrix {
    /// Number of rows (components of `f`).
    pub rows: usize,
    /// Number of columns (components of `y`).
    pub cols: usize,
    data: Vec<f64>,
}
impl JacobianMatrix {
    /// Element `(i, j)` in base units.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.cols + j]
    }
}
