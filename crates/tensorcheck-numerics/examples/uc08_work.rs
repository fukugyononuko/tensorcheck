//! UC-08: work of a spring force, W = ∫ F dx. The result is typed as an energy.

use tensorcheck::Quantity;
use tensorcheck::units::si::{Energy, Force, J, Length, N, m};
use tensorcheck_numerics::quad::{GaussKronrod, Tolerance, integrate};

fn main() {
    let k = 200.0; // N/m
    let force = |x: Quantity<Length>| -> Quantity<Force> { Quantity::new(-k * x.value_in(m), N) };
    let w = integrate(
        force,
        (0.0 * m)..(0.1 * m),
        GaussKronrod,
        Tolerance {
            rtol: 1e-12,
            atol: 1e-15 * J,
        },
    )
    .unwrap();
    let work: Quantity<Energy> = w.value;
    println!(
        "W = {} J (expected −1 J), {} evaluations",
        work.value_in(J),
        w.evaluations
    );
    assert!((work.value_in(J) + 1.0).abs() < 1e-12);
}
