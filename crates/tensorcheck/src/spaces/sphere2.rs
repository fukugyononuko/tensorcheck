//! The two-sphere of radius `a`, in angular and stereographic coordinates.

use crate::contracts::{Point, Quantity};
use crate::engine::index::ChartTransform;
use crate::engine::metric::{Metric, SymComponents};
use crate::units::si::{Dimensionless, Length};
use crate::utility::scalar::Scalar;

crate::space!(
    /// The two-sphere (reference dimension: length).
    pub Sphere2, 2, Length
);
crate::chart!(
    /// Angular coordinates `(θ, φ)`: polar angle and azimuth.
    pub Angles on Sphere2 { Theta: Dimensionless, Phi: Dimensionless }
);
crate::chart!(
    /// Stereographic coordinates `(x, y)` (projection from the south pole onto the
    /// plane tangent at the north pole, scaled so that `ρ = 2a tan(θ/2)`).
    pub Stereographic on Sphere2 { SX: Length, SY: Length }
);

/// The round metric tensor `a² (dθ² + sin²θ dφ²)`.
#[derive(Clone, Copy, Debug)]
pub struct SphereMetric {
    a2: f64,
}
impl SphereMetric {
    /// The sphere of radius `a`.
    pub fn new(radius: Quantity<Length>) -> Self {
        let a = radius.base_value();
        SphereMetric { a2: a * a }
    }
}
impl Metric for SphereMetric {
    type Chart = Angles;
    type Dim = Dimensionless;
    fn components<N: Scalar>(&self, p: &Point<Angles, N>) -> SymComponents<Angles, N> {
        let s = p.coord(Theta).sin();
        let a2 = N::from_f64(self.a2);
        [a2, N::zero(), a2 * s * s]
    }
}

/// The round metric tensor in stereographic coordinates,
/// `(1 + ρ²/4a²)^-2 (dx² + dy²)`.
#[derive(Clone, Copy, Debug)]
pub struct StereographicMetric {
    a2: f64,
}
impl StereographicMetric {
    /// The sphere of radius `a`.
    pub fn new(radius: Quantity<Length>) -> Self {
        let a = radius.base_value();
        StereographicMetric { a2: a * a }
    }
}
impl Metric for StereographicMetric {
    type Chart = Stereographic;
    type Dim = Dimensionless;
    fn components<N: Scalar>(
        &self,
        p: &Point<Stereographic, N>,
    ) -> SymComponents<Stereographic, N> {
        let (x, y) = (p.coord(SX), p.coord(SY));
        let q = N::one() + (x * x + y * y) / N::from_f64(4.0 * self.a2);
        let w = N::one() / (q * q);
        [w, N::zero(), w]
    }
}

/// Angles to stereographic: `x = 2a tan(θ/2) cos φ`, `y = 2a tan(θ/2) sin φ`.
#[derive(Clone, Copy, Debug)]
pub struct AnglesToStereographic {
    a: f64,
}
impl AnglesToStereographic {
    /// For the sphere of radius `a`.
    pub fn new(radius: Quantity<Length>) -> Self {
        AnglesToStereographic {
            a: radius.base_value(),
        }
    }
}
impl ChartTransform for AnglesToStereographic {
    type From = Angles;
    type To = Stereographic;
    fn map<N: Scalar>(&self, p: &Point<Angles, N>) -> Point<Stereographic, N> {
        let rho = N::from_f64(2.0 * self.a) * (p.coord(Theta) * N::from_f64(0.5)).tan();
        let f = p.coord(Phi);
        Point::from_base([rho * f.cos(), rho * f.sin()])
    }
    fn jacobian<N: Scalar>(&self, p: &Point<Angles, N>) -> Vec<N> {
        let half = p.coord(Theta) * N::from_f64(0.5);
        let t = half.tan();
        let rho = N::from_f64(2.0 * self.a) * t;
        let drho = N::from_f64(self.a) * (N::one() + t * t);
        let f = p.coord(Phi);
        let (s, c) = (f.sin(), f.cos());
        vec![drho * c, -(rho * s), drho * s, rho * c]
    }
}
