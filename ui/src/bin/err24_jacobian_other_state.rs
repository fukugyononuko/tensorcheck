// TC-F25 / ERR-24: the Jacobian is written for a velocity, the solve's state is a length.
// expect: error[E0277]: this Jacobian does not fit the solve: it takes other time or state types
use tensorcheck::Quantity;
use tensorcheck::units::si::{Acceleration, Length, Time, Velocity, m, s};
use tensorcheck_numerics::ode::{Jacobian, Rosenbrock, Tolerance, solve_ivp};

fn main() {
    let rhs = |_t: Quantity<Time>, y: &Quantity<Length>| -> Quantity<Velocity> { *y / (1.0 * s) };
    let jv = |_t: Quantity<Time>, _y: &Quantity<Velocity>, dy: &Quantity<Velocity>| -> Quantity<Acceleration> {
        *dy / (1.0 * s)
    };
    let method = Rosenbrock { jacobian: Jacobian::given(jv) };
    let _ = solve_ivp(rhs, (0.0 * s)..(1.0 * s), 1.0 * m, method, Tolerance { rtol: 1e-6, atol: 1e-9 * m });
}
