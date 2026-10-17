// TC-F15 / ERR-15: the absolute tolerance of a work integral given in metres.
// expect: error[E0277]: the absolute tolerance must have the same type as the state
use tensorcheck::Quantity;
use tensorcheck::units::si::{Force, Length, N, m};
use tensorcheck_numerics::quad::{GaussKronrod, Tolerance, integrate};

fn main() {
    let force = |x: Quantity<Length>| -> Quantity<Force> { Quantity::new(-200.0 * x.value_in(m), N) };
    let _ = integrate(force, (0.0 * m)..(0.1 * m), GaussKronrod, Tolerance { rtol: 1e-12, atol: 1e-15 * m });
}
