// TC-F17 / ERR-17: adding a chart-bound quantity to a tensor.
// expect: error[E0277]: cannot add a tensor and a chart-bound (non-tensor) quantity
use tensorcheck::spaces::plane::Polar;
use tensorcheck::units::si::{Length, m};
use tensorcheck::utility::typelevel::{Cons, Nil};
use tensorcheck::{ChartBound, NoSym, Tensor, Up, Vector};

fn main() {
    let v: Vector<Length, Polar> = Vector::from_components([1.0, 2.0], m);
    let w: Tensor<Length, Cons<Up<Polar>, Nil>, NoSym, ChartBound> = Tensor::from_components([1.0, 2.0], m);
    let _ = v + w;
}
