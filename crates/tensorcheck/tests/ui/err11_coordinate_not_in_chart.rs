// TC-F11 / ERR-11: asking a Cartesian vector for its theta component.
use tensorcheck::spaces::euclid3::Cartesian;
use tensorcheck::spaces::plane::Theta;
use tensorcheck::units::si::{Length, m};
use tensorcheck::Vector;

fn main() {
    let v: Vector<Length, Cartesian> = Vector::from_components([1.0, 2.0, 3.0], m);
    let _ = v.at((Theta,));
}
