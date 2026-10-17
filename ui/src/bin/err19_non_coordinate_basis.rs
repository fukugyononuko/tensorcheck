// TC-F20 / ERR-19: Christoffel symbols in a non-coordinate (orthonormal) basis.
// expect: error[E0277]: `Frame` is not a coordinate basis; Christoffel symbols and curvature need a coordinate basis
use tensorcheck::geometry::christoffel;
use tensorcheck::spaces::plane::Plane;
use tensorcheck::units::si::{Dimensionless, Length, m};
use tensorcheck::{Metric, Point, Scalar};

tensorcheck::chart!(Frame on Plane { E1: Length, E2: Length }; non_coordinate);

struct FrameMetric;
impl Metric for FrameMetric {
    type Chart = Frame;
    type Dim = Dimensionless;
    fn components<N: Scalar>(&self, _p: &Point<Frame, N>) -> [N; 3] {
        [N::one(), N::zero(), N::one()]
    }
}

fn main() {
    let p = Point::<Frame>::new((1.0 * m, 0.0 * m));
    let _ = christoffel(&FrameMetric, &p);
}
