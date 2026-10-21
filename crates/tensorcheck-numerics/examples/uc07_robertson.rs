//! UC-07: a stiff system (Robertson's chemical kinetics) with the Rosenbrock method,
//! to t = 1e11 s, compared with the published reference solution (INdAM-Bari test
//! set, problem ROBER). Solved twice: with the Jacobian by central differences, and
//! with the Jacobian given as J·dy (FR-30).

use tensorcheck::units::si::{Dimensionless, Frequency, Time, one, s};
use tensorcheck::{Quantity, Vector};
use tensorcheck_numerics::ode::{Jacobian, Method, Rosenbrock, Tolerance, solve_ivp};

tensorcheck::space!(Species, 3, Dimensionless);
tensorcheck::chart!(Conc on Species { Y1: Dimensionless, Y2: Dimensionless, Y3: Dimensionless });

type State = Vector<Dimensionless, Conc>;

fn main() {
    let rhs = |_t: Quantity<Time>, y: &State| -> Vector<Frequency, Conc> {
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
    // J·dy: the derivative of the right-hand side in the direction dy.
    let jv = |_t: Quantity<Time>, y: &State, dy: &State| -> Vector<Frequency, Conc> {
        let [_, y2, y3] = y.components_in(one);
        let [d1, d2, d3] = dy.components_in(one);
        Vector::from_components(
            [
                -0.04 * d1 + 1e4 * (y3 * d2 + y2 * d3),
                0.04 * d1 - 1e4 * (y3 * d2 + y2 * d3) - 6e7 * y2 * d2,
                6e7 * y2 * d2,
            ],
            one / s,
        )
    };
    // Table II.10.3 of the ROBER report (digits as published).
    #[allow(clippy::excessive_precision)]
    let reference = [
        0.2083340149701255e-7,
        0.8333360770334713e-13,
        0.9999999791665050,
    ];
    run(
        "numerical Jacobian",
        rhs,
        Rosenbrock {
            jacobian: Jacobian::Numerical,
        },
        reference,
    );
    run(
        "given Jacobian",
        rhs,
        Rosenbrock {
            jacobian: Jacobian::given(jv),
        },
        reference,
    );
}

fn run<F, M>(name: &str, rhs: F, method: M, reference: [f64; 3])
where
    F: FnMut(Quantity<Time>, &State) -> Vector<Frequency, Conc>,
    M: Method,
    M::Jac: tensorcheck_numerics::ode::JacobianFor<Time, State>,
{
    let y0 = Vector::from_components([1.0, 0.0, 0.0], one);
    let atol = Vector::from_components([1e-14, 1e-20, 1e-14], one);
    let sol = solve_ivp(
        rhs,
        (0.0 * s)..(1e11 * s),
        y0,
        method,
        Tolerance { rtol: 1e-8, atol },
    )
    .unwrap();
    let y = sol.y.components_in(one);
    println!(
        "{name}: y(1e11 s) = {y:?} in {} steps, {} evaluations of f",
        sol.stats.steps, sol.stats.evaluations
    );
    for i in 0..3 {
        assert!(((y[i] - reference[i]) / reference[i]).abs() < 1e-7);
    }
}
