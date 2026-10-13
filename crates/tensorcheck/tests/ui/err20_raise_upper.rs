// ERR-20: raising an index that is already upper.
use tensorcheck::spaces::plane::{Polar, PolarMetric};
use tensorcheck::units::si::{Length, m, rad};
use tensorcheck::{Point, Vector, slot};

fn main() {
    let p = Point::<Polar>::new((1.0 * m, 0.5 * rad));
    let v: Vector<Length, Polar> = Vector::from_components([1.0, 2.0], m);
    let _ = v.raise(slot::_0, &PolarMetric, &p);
}
