// TC-F28: a given Jacobian on a vector state in polar coordinates, read with `at`
// (quantities). The θ component of J·dy must have the dimension rad/s; writing c dr
// with c in 1/s gives m/s, which `from_parts` rejects (PR #52, fourth review).
// expect: error[E0277]: component `_1` has the wrong dimension for this vector
use tensorcheck::spaces::plane::{Polar, R};
use tensorcheck::units::si::{Length, Time, Velocity, one, s};
use tensorcheck::{Quantity, Vector};
use tensorcheck_numerics::ode::Jacobian;

fn main() {
    let c = 1000.0 * (one / s);
    let jv = |_t: Quantity<Time>, _y: &Vector<Length, Polar>, dy: &Vector<Length, Polar>| -> Vector<Velocity, Polar> {
        let dr = dy.at((R,));
        Vector::from_parts((c * dr, c * dr))
    };
    let _ = Jacobian::given(jv);
}
