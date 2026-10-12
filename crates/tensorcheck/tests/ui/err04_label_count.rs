// TC-F04 / ERR-04: one label for a rank-2 tensor.
use tensorcheck::labels::a;
use tensorcheck::spaces::euclid3::Cartesian;
use tensorcheck::units::si::{Dimensionless, one};
use tensorcheck::utility::typelevel::{Cons, Nil};
use tensorcheck::{Down, Tensor, Up};

fn main() {
    let t: Tensor<Dimensionless, Cons<Up<Cartesian>, Cons<Down<Cartesian>, Nil>>> =
        Tensor::from_components([[1.0; 3]; 3], one);
    let _ = t.idx(a);
}
