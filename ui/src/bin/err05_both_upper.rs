// TC-F05 / ERR-05: v^a w^a (both upper).
// expect: error[E0277]: cannot contract index `a`: it is upper in both places
use tensorcheck::labels::a;
use tensorcheck::spaces::euclid3::Cartesian;
use tensorcheck::units::si::{Length, m};
use tensorcheck::Vector;

fn main() {
    let v: Vector<Length, Cartesian> = Vector::from_components([1.0, 2.0, 3.0], m);
    let _ = (v.idx(a) * v.idx(a)).into_tensor();
}
