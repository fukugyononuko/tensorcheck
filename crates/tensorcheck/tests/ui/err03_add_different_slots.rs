// TC-F03 / ERR-03: an upper-index vector plus a lower-index covector.
use tensorcheck::spaces::euclid3::Cartesian;
use tensorcheck::units::si::{Length, m};
use tensorcheck::{Covector, Vector};

fn main() {
    let v: Vector<Length, Cartesian> = Vector::from_components([1.0, 2.0, 3.0], m);
    let w: Covector<Length, Cartesian> = Covector::from_components([1.0, 2.0, 3.0], m);
    let _ = v + w;
}
