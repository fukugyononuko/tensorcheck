// TC-F14 / ERR-14: the right-hand side returns the state itself (no division by time).
// expect: error[E0277]: the right-hand side must return the state type divided by the time dimension
use tensorcheck::Quantity;
use tensorcheck::units::si::{Length, Time, m, s};
use tensorcheck_numerics::ode::{DormandPrince54, Tolerance, solve_ivp};

fn main() {
    let rhs = |_t: Quantity<Time>, y: &Quantity<Length>| *y;
    let _ = solve_ivp(rhs, (0.0 * s)..(1.0 * s), 1.0 * m, DormandPrince54, Tolerance { rtol: 1e-6, atol: 1e-9 * m });
}
