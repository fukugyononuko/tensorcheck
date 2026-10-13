// TC-F10 / ERR-10: applying a polar-to-Cartesian transform to a Cartesian vector.
use tensorcheck::spaces::plane::{Cartesian2, Polar, PolarToCartesian2};
use tensorcheck::units::si::{Length, m, rad};
use tensorcheck::{Point, Vector};

fn main() {
    let p = Point::<Polar>::new((1.0 * m, 0.5 * rad));
    let v: Vector<Length, Cartesian2> = Vector::from_components([1.0, 2.0], m);
    let _ = v.transform(&PolarToCartesian2, &p);
}
