// TC-F13 / ERR-13: adding a symmetric tensor and a tensor without symmetry.
// expect: error[E0277]: cannot add tensors with different symmetry
use tensorcheck::spaces::euclid3::Cartesian;
use tensorcheck::units::si::{Dimensionless, one};
use tensorcheck::{Down, Rank2, Symmetric};

fn main() {
    let s: Rank2<Dimensionless, Down<Cartesian>, Down<Cartesian>, Symmetric> =
        Rank2::from_independent([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], one);
    let t: Rank2<Dimensionless, Down<Cartesian>, Down<Cartesian>> = Rank2::from_components([[1.0; 3]; 3], one);
    let _ = s + t;
}
