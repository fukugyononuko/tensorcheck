// TC-F16 / ERR-16: transforming a chart-bound (non-tensor) quantity.
use tensorcheck::spaces::plane::{Polar, PolarToCartesian2};
use tensorcheck::units::si::{Length, m, rad};
use tensorcheck::utility::typelevel::{Cons, Nil};
use tensorcheck::{ChartBound, NoSym, Point, Tensor, Up};

fn main() {
    let p = Point::<Polar>::new((1.0 * m, 0.5 * rad));
    let v: Tensor<Length, Cons<Up<Polar>, Nil>, NoSym, ChartBound> = Tensor::from_components([1.0, 2.0], m);
    let _ = v.transform(&PolarToCartesian2, &p);
}
