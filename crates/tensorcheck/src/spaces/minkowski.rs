//! Minkowski spacetime.

use crate::units::si::Length;

crate::space!(
    /// Four-dimensional spacetime.
    pub Spacetime, 4
);
crate::chart!(
    /// Inertial coordinates `(ct, x, y, z)`; the time coordinate is `x^0 = ct` (a length).
    pub Inertial on Spacetime { CT: Length, X: Length, Y: Length, Z: Length }
);
