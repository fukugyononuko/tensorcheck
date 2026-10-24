// TC-F26 / ERR-01: an entry of a given Jacobian has the wrong dimension. For the state
// (x, v), ∂(dv/dt)/∂x must have the dimension 1/s²; writing c dx with c in 1/s
// makes a term of dimension m/s, which cannot be added to c dv (m/s²).
// The sum then has no valid type, so the return type is reported as well.
// expect: error[E0277]: cannot add or subtract quantities with different dimensions
// expect: error[E0277]: the Jacobian must return J·dy with the type of the right-hand side
use tensorcheck::Quantity;
use tensorcheck::units::si::{Length, Time, Velocity, one, s};
use tensorcheck_numerics::ode::Jacobian;

type Osc = (Quantity<Length>, Quantity<Velocity>);

fn main() {
    let c = 1000.0 * (one / s);
    let jv = |_t: Quantity<Time>, _y: &Osc, dy: &Osc| (dy.1, -(c * dy.0) - c * dy.1);
    let _ = Jacobian::given(jv);
}
