// TC-F24 / ERR-23: the Jacobian returns dy itself (no division by time).
// expect: error[E0277]: the Jacobian must return J·dy with the type of the right-hand side
use tensorcheck::Quantity;
use tensorcheck::units::si::{Length, Time};
use tensorcheck_numerics::ode::Jacobian;

fn main() {
    let jv = |_t: Quantity<Time>, _y: &Quantity<Length>, dy: &Quantity<Length>| *dy;
    let _ = Jacobian::given(jv);
}
