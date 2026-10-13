// TC-F09 / ERR-09: lowering a Cartesian index with the polar metric tensor.
use tensorcheck::spaces::plane::{Cartesian2, Polar, PolarMetric};
use tensorcheck::units::si::{Length, m, rad};
use tensorcheck::{Point, Vector, slot};

fn main() {
    let p = Point::<Polar>::new((1.0 * m, 0.5 * rad));
    let v: Vector<Length, Cartesian2> = Vector::from_components([1.0, 2.0], m);
    let _ = v.lower(slot::_0, &PolarMetric, &p);
}
