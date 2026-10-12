// TC-F06 / ERR-06: a label used three times.
use tensorcheck::labels::a;
use tensorcheck::spaces::euclid3::Cartesian;
use tensorcheck::units::si::{Length, m};
use tensorcheck::{Covector, Vector};

fn main() {
    let v: Vector<Length, Cartesian> = Vector::from_components([1.0, 2.0, 3.0], m);
    let w: Covector<Length, Cartesian> = Covector::from_components([1.0, 2.0, 3.0], m);
    let _ = (v.idx(a) * w.idx(a) * v.idx(a)).into_tensor();
}
