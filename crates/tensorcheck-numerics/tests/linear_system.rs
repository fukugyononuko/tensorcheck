//! The linear system `E = I/(hγ) − J` of the Rosenbrock method (issue #54): failures
//! of TC-U27. Overflow in the elimination of a finite `E` is not `Singular`, and a
//! singular `E` is retried with a smaller step.

use tensorcheck::units::si::{Dimensionless, Frequency, Time, one, s};
use tensorcheck::utility::diagnostics::{JacobianReturns, RhsReturns};
use tensorcheck::{PerDim, Quantity, Vector, VectorSpace};
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
fn a_matrix_factored_without_scaling_is_factored_as_before() {
    // From y0 = 0 with f = 0 the first step is 1e-6 s, so 1/(hγ) = 4e6 and
    // E = [[MAX, 2^-1074], [MAX, 0]]: regular, and eliminated without overflow.
    // Scaling its rows by 2^-1024 rounds 2^-1074 to 0 and makes the rows equal
    // (PR #59 review). Only a given matrix: central differences would round the
    // tiny entry to 0 anyway.
    // No step is rejected: a smaller step would make E regular after scaling too,
    // which would hide a factorization that differs from before.
    let tiny = f64::MIN_POSITIVE * f64::EPSILON;
    let a = [-MAX, -tiny, -MAX, 4e6];
    let sol = solve_ivp(
        |_t: Quantity<Time>, _y: &Z| Vector::from_components([0.0, 0.0], one / s),
        (0.0 * s)..(1.0 * s),
        Vector::from_components([0.0, 0.0], one),
        Rosenbrock {
            jacobian: Jacobian::given(move |_t: Quantity<Time>, _y: &Z, dy: &Z| product(a, dy)),
        },
        Tolerance {
            rtol: 1e-6,
            atol: Vector::from_components([1e-310, 1e-310], one),
        },
    )
    .unwrap();
    assert_eq!(sol.y.components_in(one), [0.0, 0.0]);
    assert_eq!(sol.stats.rejected, 0);
}

#[test]
fn a_matrix_singular_for_one_step_is_retried_with_a_smaller_step() {
    // y' = 4e6 y from 0: at the first step 1/(hγ) = 4e6, so E = 0 exactly. Half
    // the step makes E = 4e6 I (PR #59 review; RODAS halves the step too).
    for given in [true, false] {
        let r = solve([4e6, 0.0, 0.0, 4e6], [0.0, 0.0], given);
        let y = r.unwrap_or_else(|e| panic!("given={given}: {e:?}"));
        assert_eq!(y.components_in(one), [0.0, 0.0], "given={given}");
    }
}

#[test]
fn the_fifth_singular_matrix_ends_the_solve() {
    // E rounds to −A, singular for every step tried: one product per unit vector
    // per attempt, and the fifth attempt gives up (`NSING` in `rodas.f`).
    let a = [1e100; 4];
    let calls = std::cell::Cell::new(0);
    let jv = |_t: Quantity<Time>, _y: &Z, dy: &Z| {
        calls.set(calls.get() + 1);
        product(a, dy)
    };
    let r = solve_ivp(
        move |_t: Quantity<Time>, y: &Z| product(a, y),
        (0.0 * s)..(1.0 * s),
        Vector::from_components([0.0, 0.0], one),
        Rosenbrock {
            jacobian: Jacobian::given(jv),
        },
        Tolerance {
            rtol: 1e-6,
            atol: Vector::from_components([1e-6, 1e-6], one),
        },
    );
    assert_eq!(r.unwrap_err(), SolveError::Singular { t: 0.0 * s });
    assert_eq!(calls.get(), 5 * 2);
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

/// A state of `N` components without dimensions (a user's type, EXT).
const N: usize = 1025;
#[derive(Clone, Copy)]
struct Big([f64; N]);
impl VectorSpace for Big {
    const LEN: usize = N;
    fn component(&self, i: usize) -> f64 {
        self.0[i]
    }
    fn from_components(mut f: impl FnMut(usize) -> f64) -> Self {
        Big(core::array::from_fn(|i| f(i)))
    }
}
impl PerDim<Time> for Big {
    type Out = Big;
}
impl RhsReturns<Big> for Big {}
impl JacobianReturns<Big> for Big {}

/// `W v` for Wilkinson's matrix of size N: 1 on the diagonal and in the last
/// column, −1 below the diagonal. Partial pivoting does not swap its rows, and the
/// last column doubles at each elimination, to 2^(N−1) = 2^1024.
fn wilkinson(v: &Big) -> Big {
    let mut below = 0.0; // Σ_{j<i} v_j
    Big(core::array::from_fn(|i| {
        let last = if i + 1 < N { v.0[N - 1] } else { 0.0 };
        let w = v.0[i] + last - below;
        below += v.0[i];
        w
    }))
}

#[test]
fn a_finite_matrix_whose_elimination_overflows_is_retried_with_a_smaller_step() {
    // y' = J y from 0 with J = 4e6 I − W: the first step is 1e-6 s, so E = W, whose
    // elimination overflows even with the rows scaled (PR #59 review). Half the
    // step gives E = 4e6 I + W, which is diagonally dominant. The interval is just
    // the first step, so that it is tried, as each factorization takes seconds
    // unoptimized.
    let jv = |_t: Quantity<Time>, _y: &Big, dy: &Big| {
        let w = wilkinson(dy);
        Big(core::array::from_fn(|i| 4e6 * dy.0[i] - w.0[i]))
    };
    let sol = solve_ivp(
        move |t: Quantity<Time>, y: &Big| jv(t, y, y),
        (0.0 * s)..(1e-6 * s),
        Big([0.0; N]),
        Rosenbrock {
            jacobian: Jacobian::given(jv),
        },
        Tolerance {
            rtol: 1e-6,
            atol: Big([1e-6; N]),
        },
    )
    .unwrap();
    assert!(sol.y.0.iter().all(|&x| x == 0.0));
    assert_eq!(
        (sol.stats.rejected, sol.stats.steps),
        (1, 2),
        "{:?}",
        sol.stats
    );
}
