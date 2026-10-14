//! UC-01: Newton's second law with checked dimensions.
//!
//! `force + acceleration` would not compile: the dimensions differ.

use tensorcheck::Vector;
use tensorcheck::spaces::euclid3::Cartesian;
use tensorcheck::units::si::{Acceleration, Force, N, kg, m, s};

fn main() {
    let mass = 2.0 * kg;
    let acc: Vector<Acceleration, Cartesian> =
        Vector::from_components([0.0, 0.0, -9.8], m / s.pow::<2>());
    let force: Vector<Force, Cartesian> = mass * acc;
    let f = force.components_in(N);
    println!("F = {f:?} N");
    assert_eq!(f, [0.0, 0.0, -19.6]);
}
