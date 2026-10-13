//! Minkowski spacetime.

use crate::contracts::Point;
use crate::engine::metric::{Metric, SymComponents};
use crate::units::si::{Dimensionless, Length};
use crate::utility::scalar::Scalar;

crate::space!(
    /// Four-dimensional spacetime (reference dimension: length).
    pub Spacetime, 4, Length
);
crate::chart!(
    /// Inertial coordinates `(ct, x, y, z)`; the time coordinate is `x^0 = ct` (a length).
    pub Inertial on Spacetime { CT: Length, X: Length, Y: Length, Z: Length }
);

/// The Minkowski metric tensor `η = diag(+1, -1, -1, -1)` in inertial coordinates.
///
/// For the opposite signature, define your own metric (see the volatility test
/// TC-R05).
#[derive(Clone, Copy, Debug)]
pub struct MinkowskiMetric;
impl Metric for MinkowskiMetric {
    type Chart = Inertial;
    type Dim = Dimensionless;
    fn components<N: Scalar>(&self, _p: &Point<Inertial, N>) -> SymComponents<Inertial, N> {
        let (o, z) = (N::one(), N::zero());
        [o, z, z, z, -o, z, z, -o, z, -o]
    }
    fn inverse<N: Scalar>(&self, p: &Point<Inertial, N>) -> Option<SymComponents<Inertial, N>> {
        Some(self.components(p))
    }
}
