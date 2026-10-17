// TC-F01 / ERR-01: force + acceleration.
// expect: error[E0277]: cannot add or subtract quantities with different dimensions
use tensorcheck::units::si::{kg, m, s};

fn main() {
    let acc = 9.8 * (m / s.pow::<2>());
    let force = (2.0 * kg) * acc;
    let _ = force + acc;
}
