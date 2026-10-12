//! Three-dimensional Euclidean space.

use crate::units::si::Length;

crate::space!(
    /// Three-dimensional Euclidean space.
    pub Euclid3, 3
);
crate::chart!(
    /// Cartesian coordinates `(x, y, z)`.
    pub Cartesian on Euclid3 { X: Length, Y: Length, Z: Length }
);
