// TC-F08 / ERR-08: contracting a polar index with a Cartesian index.
// expect: error[E0277]: index `a` is in chart `Polar` in one place and `Cartesian2` in the other
use tensorcheck::labels::a;
use tensorcheck::spaces::plane::{Cartesian2, Polar};
use tensorcheck::units::si::{Dimensionless, one};
use tensorcheck::{Covector, Vector};

fn main() {
    let v: Vector<Dimensionless, Polar> = Vector::from_components([1.0, 0.5], one);
    let w: Covector<Dimensionless, Cartesian2> = Covector::from_components([1.0, 2.0], one);
    let _ = (v.idx(a) * w.idx(a)).into_tensor();
}
