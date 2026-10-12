//! The Euclidean plane.

use crate::units::si::{Dimensionless, Length};

crate::space!(
    /// The two-dimensional Euclidean plane.
    pub Plane, 2
);
crate::chart!(
    /// Cartesian coordinates `(x, y)` on the plane.
    pub Cartesian2 on Plane { X: Length, Y: Length }
);
crate::chart!(
    /// Polar coordinates `(r, θ)`; `θ` is dimensionless (radians).
    pub Polar on Plane { R: Length, Theta: Dimensionless }
);
