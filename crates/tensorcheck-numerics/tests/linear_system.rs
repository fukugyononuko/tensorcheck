//! The linear system `E = I/(hγ) − J` of the Rosenbrock method (issue #54): failures
//! of TC-U27. `Singular` only when `E`, as stored, is singular; overflow is not that.

use tensorcheck::units::si::{Dimensionless, Frequency, Time, one, s};
use tensorcheck::{Quantity, Vector};
use tensorcheck_numerics::ode::{Jacobian, Rosenbrock, SolveError, Tolerance, solve_ivp};

tensorcheck::space!(Plane, 2, Dimensionless);
tensorcheck::chart!(Xy on Plane { X: Dimensionless, Y: Dimensionless });
type Z = Vector<Dimensionless, Xy>;
type Rate = Vector<Frequency, Xy>;

const MAX: f64 = f64::MAX;

/// `A·v` for the row-major 2×2 matrix `a`.
fn product(a: [f64; 4], v: &Z) -> Rate {
    let [x, y] = v.components_in(one);
    Vector::from_components([a[0] * x + a[1] * y, a[2] * x + a[3] * y], one / s)
}

/// Solves `y' = A y` from `y0` over 1 s, with `A` given or by central differences.
fn solve(a: [f64; 4], y0: [f64; 2], given: bool) -> Result<Z, SolveError<Time>> {
    let rhs = move |_t: Quantity<Time>, y: &Z| product(a, y);
    let y0 = Vector::from_components(y0, one);
    let tol = Tolerance {
        rtol: 1e-6,
        atol: Vector::from_components([1e-310, 1e-310], one),
    };
    let span = (0.0 * s)..(1.0 * s);
    let r = if given {
        let jv = move |_t: Quantity<Time>, _y: &Z, dy: &Z| product(a, dy);
        solve_ivp(
            rhs,
            span,
            y0,
            Rosenbrock {
                jacobian: Jacobian::given(jv),
            },
            tol,
        )
    } else {
        solve_ivp(
            rhs,
            span,
            y0,
            Rosenbrock {
                jacobian: Jacobian::Numerical,
            },
            tol,
        )
    };
    r.map(|sol| sol.y)
}

#[test]
fn a_finite_regular_matrix_is_not_singular_whatever_its_scale() {
    // Eliminating these overflowed before the rows were scaled (issue #54): the
    // first pivot is ±MAX and the second becomes 2 MAX. Every matrix is regular:
    // det = 2 a² for the rotations, a² for the diagonal.
    for a in [MAX, 1e300, 1.0, 1e-300] {
        let matrices = [
            [-a, -a, a, -a],
            [-a, a, -a, -a],
            [-a, 0.0, 0.0, -a],
            [-a, -a / 2.0, a, -a],
        ];
        // E is factored at every step whatever the state. From a state other than
        // 0, a ≥ 1e300 stops at the initial step-size estimate instead (issue #58).
        let y0 = if a >= 1e300 { [0.0, 0.0] } else { [1.0, 2.0] };
        for m in matrices {
            for given in [true, false] {
                let r = solve(m, y0, given);
                let y = r.unwrap_or_else(|e| panic!("{m:?} given={given}: {e:?}"));
                // The solution decays: |y(1 s)| ≤ |y0| (the eigenvalues have
                // negative real parts).
                let [x, yy] = y.components_in(one);
                assert!(
                    x.is_finite() && yy.is_finite() && x.hypot(yy) <= y0[0].hypot(y0[1]),
                    "{m:?} given={given}: {x}, {yy}"
                );
            }
        }
    }
}

#[test]
fn a_matrix_singular_as_stored_is_singular() {
    // E = I/(hγ) − A rounds to −A when 1/(hγ) is below half an ulp of a, and
    // −A has equal rows.
    for a in [MAX, 1e100] {
        for m in [[a, a, a, a], [-a, a, -a, a]] {
            for given in [true, false] {
                let r = solve(m, [1e-300 / a, 1e-300 / a], given);
                assert_eq!(
                    r.unwrap_err(),
                    SolveError::Singular { t: 0.0 * s },
                    "{m:?} given={given}"
                );
            }
        }
    }
}

#[test]
fn a_matrix_that_is_not_finite_is_non_finite() {
    // f = (MAX tanh(1000 y), 0) is 0 at y = 0, so the first step is of order one,
    // but ∂f/∂y = 1000 MAX is not representable and the central differences
    // overflow. (A given matrix is checked before the factorization, DD-15.)
    let rhs = |_t: Quantity<Time>, z: &Z| -> Rate {
        let [_, y] = z.components_in(one);
        Vector::from_components([MAX * (1000.0 * y).tanh(), 0.0], one / s)
    };
    let r = solve_ivp(
        rhs,
        (0.0 * s)..(1.0 * s),
        Vector::from_components([0.0, 0.0], one),
        Rosenbrock {
            jacobian: Jacobian::Numerical,
        },
        Tolerance {
            rtol: 1e-6,
            atol: Vector::from_components([1e-6, 1e-6], one),
        },
    );
    assert_eq!(r.unwrap_err(), SolveError::NonFinite { t: 0.0 * s });
}
