// TC-F19: correct expressions from the use cases compile.
// expect: compiles
use tensorcheck::labels::{a, b};
use tensorcheck::spaces::euclid3::Cartesian;
use tensorcheck::spaces::plane::{Polar, PolarMetric, PolarToCartesian2};
use tensorcheck::units::si::{Acceleration, Force, Velocity, kg, m, rad, s};
use tensorcheck::{Point, Vector, slot};

fn main() {
    let acc: Vector<Acceleration, Cartesian> = Vector::from_components([0.0, 0.0, -9.8], m / s.pow::<2>());
    let _force: Vector<Force, Cartesian> = (2.0 * kg) * acc;
    let p = Point::<Polar>::new((2.0 * m, 0.5 * rad));
    let v: Vector<Velocity, Polar> = Vector::from_parts((1.0 * (m / s), 0.25 * (rad / s)));
    let _v2 = (v.lower(slot::_0, &PolarMetric, &p).idx(a) * v.idx(a)).into_tensor();
    let _vc = v.transform(&PolarToCartesian2, &p).unwrap();
    let _ = b;
}
