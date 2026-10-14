//! TC-N01 (numerics): RK4 on a dimensioned oscillator state versus the same loop on
//! `f64`. Run with `cargo bench -p tensorcheck-numerics --bench overhead`.

use std::hint::black_box;
use std::time::Instant;

use tensorcheck::Quantity;
use tensorcheck::units::si::{Acceleration, Length, Time, Velocity, m, s};
use tensorcheck_numerics::ode::{Rk4, Tolerance, solve_ivp};

fn best_of<F: FnMut() -> f64>(rounds: usize, mut f: F) -> f64 {
    (0..rounds).map(|_| f()).fold(f64::INFINITY, f64::min)
}

fn main() {
    let steps = 90_000; // below the step limit (ode::MAX_STEPS)
    let h = 1e-4;
    let checked = best_of(5, || {
        let t = Instant::now();
        let rhs = |_t: Quantity<Time>,
                   y: &(Quantity<Length>, Quantity<Velocity>)|
         -> (Quantity<Velocity>, Quantity<Acceleration>) {
            (y.1, Quantity::of(-y.0.value_in(m), m / (s * s)))
        };
        let sol = solve_ivp(
            rhs,
            (0.0 * s)..((steps as f64 * h) * s),
            black_box((1.0 * m, 0.0 * (m / s))),
            Rk4 { h },
            Tolerance {
                rtol: 1e-6,
                atol: (1e-9 * m, 1e-9 * (m / s)),
            },
        )
        .unwrap();
        black_box(sol.y.0.value_in(m));
        t.elapsed().as_secs_f64()
    });
    let hand = best_of(5, || {
        let t = Instant::now();
        let f = |y: [f64; 2]| [y[1], -y[0]];
        let mut y = black_box([1.0, 0.0]);
        for _ in 0..steps {
            let k1 = f(y);
            let k2 = f([y[0] + 0.5 * h * k1[0], y[1] + 0.5 * h * k1[1]]);
            let k3 = f([y[0] + 0.5 * h * k2[0], y[1] + 0.5 * h * k2[1]]);
            let k4 = f([y[0] + h * k3[0], y[1] + h * k3[1]]);
            for i in 0..2 {
                y[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            }
        }
        black_box(y);
        t.elapsed().as_secs_f64()
    });
    println!(
        "RK4, 2-component state, {steps} steps: ratio {:.2}",
        checked / hand
    );
}
