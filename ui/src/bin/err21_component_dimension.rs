// ERR-21: the θ component of a polar velocity given in m/s instead of 1/s.
// expect: error[E0277]: component `_1` has the wrong dimension for this vector
use tensorcheck::spaces::plane::Polar;
use tensorcheck::units::si::{Velocity, m, s};
use tensorcheck::Vector;

fn main() {
    let _: Vector<Velocity, Polar> = Vector::from_parts((1.0 * (m / s), 0.25 * (m / s)));
}
