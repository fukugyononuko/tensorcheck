//! UC-07: a stiff system (Robertson's chemical kinetics) with the Rosenbrock method,
//! to t = 1e11 s, compared with the published reference solution (INdAM-Bari test
//! set, problem ROBER).

use tensorcheck::units::si::{Dimensionless, Frequency, Time, one, s};
use tensorcheck::{Quantity, Vector};
use tensorcheck_numerics::ode::{Jacobian, Rosenbrock, Tolerance, solve_ivp};

tensorcheck::space!(Species, 3, Dimensionless);
tensorcheck::chart!(Conc on Species { Y1: Dimensionless, Y2: Dimensionless, Y3: Dimensionless });

fn main() {
    let rhs = |_t: Quantity<Time>, y: &Vector<Dimensionless, Conc>| -> Vector<Frequency, Conc> {
        let [y1, y2, y3] = y.components_in(one);
        Vector::from_components(
            [
                -0.04 * y1 + 1e4 * y2 * y3,
                0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2 * y2,
                3e7 * y2 * y2,
            ],
            one / s,
        )
    };
    let y0 = Vector::from_components([1.0, 0.0, 0.0], one);
    let atol = Vector::from_components([1e-14, 1e-20, 1e-14], one);
    let sol = solve_ivp(
        rhs,
        (0.0 * s)..(1e11 * s),
        y0,
        Rosenbrock {
            jacobian: Jacobian::Numerical,
        },
        Tolerance { rtol: 1e-8, atol },
    )
    .unwrap();
    let y = sol.y.components_in(one);
    // Table II.10.3 of the ROBER report (digits as published).
    #[allow(clippy::excessive_precision)]
    let reference = [
        0.2083340149701255e-7,
        0.8333360770334713e-13,
        0.9999999791665050,
    ];
    println!("y(1e11 s) = {y:?} in {} steps", sol.stats.steps);
    for i in 0..3 {
        assert!(((y[i] - reference[i]) / reference[i]).abs() < 1e-7);
    }
}
