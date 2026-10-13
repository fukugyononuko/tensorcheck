//! Three-dimensional Euclidean space: Cartesian, spherical and cylindrical charts.

use crate::contracts::Point;
use crate::engine::metric::{Metric, SymComponents};
use crate::units::si::{Dimensionless, Length};
use crate::utility::scalar::Scalar;

crate::space!(
    /// Three-dimensional Euclidean space (reference dimension: length).
    pub Euclid3, 3, Length
);
crate::chart!(
    /// Cartesian coordinates `(x, y, z)`.
    pub Cartesian on Euclid3 { X: Length, Y: Length, Z: Length }
);

/// The Euclidean metric tensor `δ_ij` in Cartesian coordinates.
#[derive(Clone, Copy, Debug)]
pub struct EuclidMetric;
impl Metric for EuclidMetric {
    type Chart = Cartesian;
    type Dim = Dimensionless;
    fn components<N: Scalar>(&self, _p: &Point<Cartesian, N>) -> SymComponents<Cartesian, N> {
        let (o, z) = (N::one(), N::zero());
        [o, z, z, o, z, o]
    }
}

/// Spherical coordinates `(r, θ, φ)`.
pub mod spherical {
    use super::{Cartesian, Euclid3};
    use crate::contracts::Point;
    use crate::engine::index::ChartTransform;
    use crate::engine::metric::{Metric, SymComponents};
    use crate::units::si::{Dimensionless, Length};
    use crate::utility::scalar::Scalar;

    crate::chart!(
        /// Spherical coordinates `(r, θ, φ)`: radius, polar angle, azimuth.
        pub Spherical on Euclid3 { R: Length, Theta: Dimensionless, Phi: Dimensionless }
    );

    /// The Euclidean metric tensor in spherical coordinates, `diag(1, r², r² sin²θ)`.
    #[derive(Clone, Copy, Debug)]
    pub struct SphericalMetric;
    impl Metric for SphericalMetric {
        type Chart = Spherical;
        type Dim = Dimensionless;
        fn components<N: Scalar>(&self, p: &Point<Spherical, N>) -> SymComponents<Spherical, N> {
            let (r, t) = (p.coord(R), p.coord(Theta));
            let s = t.sin();
            let z = N::zero();
            [N::one(), z, z, r * r, z, r * r * s * s]
        }
    }

    /// Spherical to Cartesian.
    #[derive(Clone, Copy, Debug)]
    pub struct SphericalToCartesian;
    impl ChartTransform for SphericalToCartesian {
        type From = Spherical;
        type To = Cartesian;
        fn map<N: Scalar>(&self, p: &Point<Spherical, N>) -> Point<Cartesian, N> {
            let (r, t, f) = (p.coord(R), p.coord(Theta), p.coord(Phi));
            Point::from_base([r * t.sin() * f.cos(), r * t.sin() * f.sin(), r * t.cos()])
        }
        fn jacobian<N: Scalar>(&self, p: &Point<Spherical, N>) -> Vec<N> {
            let (r, t, f) = (p.coord(R), p.coord(Theta), p.coord(Phi));
            let (st, ct, sf, cf) = (t.sin(), t.cos(), f.sin(), f.cos());
            vec![
                st * cf,
                r * ct * cf,
                -(r * st * sf),
                st * sf,
                r * ct * sf,
                r * st * cf,
                ct,
                -(r * st),
                N::zero(),
            ]
        }
    }
}

/// Cylindrical coordinates `(ρ, φ, z)`.
pub mod cylindrical {
    use super::{Cartesian, Euclid3};
    use crate::contracts::Point;
    use crate::engine::index::ChartTransform;
    use crate::engine::metric::{Metric, SymComponents};
    use crate::units::si::{Dimensionless, Length};
    use crate::utility::scalar::Scalar;

    crate::chart!(
        /// Cylindrical coordinates `(ρ, φ, z)`.
        pub Cylindrical on Euclid3 { Rho: Length, Phi: Dimensionless, Z: Length }
    );

    /// The Euclidean metric tensor in cylindrical coordinates, `diag(1, ρ², 1)`.
    #[derive(Clone, Copy, Debug)]
    pub struct CylindricalMetric;
    impl Metric for CylindricalMetric {
        type Chart = Cylindrical;
        type Dim = Dimensionless;
        fn components<N: Scalar>(
            &self,
            p: &Point<Cylindrical, N>,
        ) -> SymComponents<Cylindrical, N> {
            let rho = p.coord(Rho);
            let (o, z) = (N::one(), N::zero());
            [o, z, z, rho * rho, z, o]
        }
    }

    /// Cylindrical to Cartesian.
    #[derive(Clone, Copy, Debug)]
    pub struct CylindricalToCartesian;
    impl ChartTransform for CylindricalToCartesian {
        type From = Cylindrical;
        type To = Cartesian;
        fn map<N: Scalar>(&self, p: &Point<Cylindrical, N>) -> Point<Cartesian, N> {
            let (rho, f, z) = (p.coord(Rho), p.coord(Phi), p.coord(Z));
            Point::from_base([rho * f.cos(), rho * f.sin(), z])
        }
        fn jacobian<N: Scalar>(&self, p: &Point<Cylindrical, N>) -> Vec<N> {
            let (rho, f) = (p.coord(Rho), p.coord(Phi));
            let (s, c) = (f.sin(), f.cos());
            let (o, z) = (N::one(), N::zero());
            vec![c, -(rho * s), z, s, rho * c, z, z, z, o]
        }
    }
}
