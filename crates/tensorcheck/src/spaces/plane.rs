//! The Euclidean plane: Cartesian and polar charts, metric tensors and transforms.

use crate::access::component::Extent;
use crate::contracts::Point;
use crate::engine::index::ChartTransform;
use crate::engine::metric::{Metric, SymComponents};
use crate::units::si::{Dimensionless, Length};
use crate::utility::scalar::Scalar;

crate::space!(
    /// The two-dimensional Euclidean plane (reference dimension: length).
    pub Plane, 2, Length
);
crate::chart!(
    /// Cartesian coordinates `(x, y)` on the plane.
    pub Cartesian2 on Plane { X: Length, Y: Length }
);
crate::chart!(
    /// Polar coordinates `(r, θ)`; `θ` is dimensionless (radians).
    pub Polar on Plane { R: Length, Theta: Dimensionless }
);

/// The Euclidean metric tensor `δ_ab` in Cartesian coordinates.
#[derive(Clone, Copy, Debug)]
pub struct EuclidMetric2;
impl Metric for EuclidMetric2 {
    type Chart = Cartesian2;
    type Dim = Dimensionless;
    fn components<N: Scalar>(&self, _p: &Point<Cartesian2, N>) -> SymComponents<Cartesian2, N> {
        [N::one(), N::zero(), N::one()]
    }
}

/// The Euclidean metric tensor in polar coordinates, `diag(1, r²)`.
#[derive(Clone, Copy, Debug)]
pub struct PolarMetric;
impl Metric for PolarMetric {
    type Chart = Polar;
    type Dim = Dimensionless;
    fn components<N: Scalar>(&self, p: &Point<Polar, N>) -> SymComponents<Polar, N> {
        let r = p.coord(R);
        [N::one(), N::zero(), r * r]
    }
}

/// Polar to Cartesian: `x = r cos θ`, `y = r sin θ`.
#[derive(Clone, Copy, Debug)]
pub struct PolarToCartesian2;
impl ChartTransform for PolarToCartesian2 {
    type From = Polar;
    type To = Cartesian2;
    fn map<N: Scalar>(&self, p: &Point<Polar, N>) -> Point<Cartesian2, N> {
        let (r, t) = (p.coord(R), p.coord(Theta));
        Point::from_base([r * t.cos(), r * t.sin()])
    }
    fn jacobian<N: Scalar>(&self, p: &Point<Polar, N>) -> Vec<N> {
        let (r, t) = (p.coord(R), p.coord(Theta));
        let (c, s) = (t.cos(), t.sin());
        vec![c, -(r * s), s, r * c]
    }
    fn inverse_jacobian<N: Scalar>(&self, p: &Point<Polar, N>) -> Option<Vec<N>> {
        let (r, t) = (p.coord(R), p.coord(Theta));
        let (c, s) = (t.cos(), t.sin());
        Some(vec![c, s, -(s / r), c / r])
    }
}

/// Cartesian to polar: `r = √(x² + y²)`, `θ = atan2(y, x)`.
#[derive(Clone, Copy, Debug)]
pub struct Cartesian2ToPolar;
impl ChartTransform for Cartesian2ToPolar {
    type From = Cartesian2;
    type To = Polar;
    fn map<N: Scalar>(&self, p: &Point<Cartesian2, N>) -> Point<Polar, N> {
        let (x, y) = (p.coord(X), p.coord(Y));
        Point::from_base([(x * x + y * y).sqrt(), y.atan2(x)])
    }
    fn jacobian<N: Scalar>(&self, p: &Point<Cartesian2, N>) -> Vec<N> {
        let (x, y) = (p.coord(X), p.coord(Y));
        let r2 = x * x + y * y;
        let r = r2.sqrt();
        vec![x / r, y / r, -(y / r2), x / r2]
    }
}

const _: () = assert!(<Plane as Extent>::N == 2);
