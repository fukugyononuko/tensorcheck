//! UC-04: changing charts. A polar vector transformed to Cartesian coordinates can
//! be added to a Cartesian vector; adding it untransformed would not compile.

use tensorcheck::spaces::plane::{Cartesian2, Polar, PolarToCartesian2};
use tensorcheck::units::si::{Length, m, rad};
use tensorcheck::{Point, Vector};

fn main() {
    let p = Point::<Polar>::new((2.0 * m, std::f64::consts::FRAC_PI_2 * rad));
    let radial: Vector<Length, Polar> = Vector::from_components([1.0, 0.0], m);
    let as_cartesian = radial.transform(&PolarToCartesian2, &p).unwrap();
    let w: Vector<Length, Cartesian2> = Vector::from_components([1.0, 0.0], m);
    let sum = (as_cartesian + w).components_in(m);
    println!("∂_r at θ = π/2 plus (1, 0) m = {sum:?} m");
    assert!((sum[0] - 1.0).abs() < 1e-15 && (sum[1] - 1.0).abs() < 1e-15);
}
