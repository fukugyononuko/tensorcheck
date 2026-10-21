// TC-F27 / ERR-24: a raw matrix of numbers cannot be passed as a Jacobian by
// implementing the solver's trait outside the crate; it would skip the dimension
// checks (PR #52 review).
// expect: error[E0277]: this is not a Jacobian for this solve
use tensorcheck::Quantity;
use tensorcheck::units::si::{Length, Time, Velocity, m, s};
use tensorcheck_numerics::ode::{JacobianFor, Rosenbrock, Tolerance, solve_ivp};

#[derive(Clone, Copy)]
struct Raw;
impl JacobianFor<Time, Quantity<Length>> for Raw {
    fn fill(&mut self, _t: f64, _y: &[f64], out: &mut [f64]) -> bool {
        out[0] = 123.0;
        true
    }
}

fn main() {
    let rhs = |_t: Quantity<Time>, y: &Quantity<Length>| -> Quantity<Velocity> { *y / (1.0 * s) };
    let method = Rosenbrock { jacobian: Raw };
    let _ = solve_ivp(rhs, (0.0 * s)..(1.0 * s), 1.0 * m, method, Tolerance { rtol: 1e-6, atol: 1e-9 * m });
}
