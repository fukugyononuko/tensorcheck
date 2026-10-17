// TC-F12 / ERR-12: a symmetric tensor with one upper and one lower index.
// expect: error[E0277]: only tensors whose indices are all of the same kind can be symmetric or antisymmetric
use tensorcheck::spaces::euclid3::Cartesian;
use tensorcheck::units::si::{Dimensionless, one};
use tensorcheck::{Down, Rank2, Symmetric, Up};

fn main() {
    let _: Rank2<Dimensionless, Up<Cartesian>, Down<Cartesian>, Symmetric> =
        Rank2::from_independent([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], one);
}
