//! One step of an ODE method (CMP-E05, DD-12). Plain arrays only:
//! `f(t, y, out)` writes `dy/dt` into `out`.

// Index loops mirror the formulas (several arrays share the index).
#![allow(clippy::needless_range_loop)]

/// Work arrays for a method with `stages` stages and `n` components.
pub struct Work {
    /// Stage derivatives `k_i`.
    pub k: Vec<Vec<f64>>,
    /// Temporary state.
    pub tmp: Vec<f64>,
}
impl Work {
    /// Allocates once per solve (RES-02).
    pub fn new(stages: usize, n: usize) -> Self {
        Work {
            k: vec![vec![0.0; n]; stages],
            tmp: vec![0.0; n],
        }
    }
}

/// The classical fourth-order Runge–Kutta step: `y_new = y + h Σ b_i k_i`.
pub fn rk4<F: FnMut(f64, &[f64], &mut [f64]) + ?Sized>(
    f: &mut F,
    t: f64,
    y: &[f64],
    h: f64,
    w: &mut Work,
    y_new: &mut [f64],
) {
    let n = y.len();
    f(t, y, &mut w.k[0]);
    for i in 0..n {
        w.tmp[i] = y[i] + 0.5 * h * w.k[0][i];
    }
    let (k01, rest) = w.k.split_at_mut(1);
    f(t + 0.5 * h, &w.tmp, &mut rest[0]);
    for i in 0..n {
        w.tmp[i] = y[i] + 0.5 * h * rest[0][i];
    }
    let (k12, rest2) = rest.split_at_mut(1);
    f(t + 0.5 * h, &w.tmp, &mut rest2[0]);
    for i in 0..n {
        w.tmp[i] = y[i] + h * rest2[0][i];
    }
    let (k2, k3) = rest2.split_at_mut(1);
    f(t + h, &w.tmp, &mut k3[0]);
    for i in 0..n {
        y_new[i] = y[i] + h / 6.0 * (k01[0][i] + 2.0 * k12[0][i] + 2.0 * k2[0][i] + k3[0][i]);
    }
}

/// Cubic Hermite interpolation on a step of size `h` (DD-12, FR-31): writes the state
/// at `t + θ h` from the values `y0`, `y1` and derivatives `f0`, `f1` at both ends.
pub fn hermite(
    theta: f64,
    h: f64,
    y0: &[f64],
    f0: &[f64],
    y1: &[f64],
    f1: &[f64],
    out: &mut [f64],
) {
    let s = 1.0 - theta;
    let (h00, h10) = ((1.0 + 2.0 * theta) * s * s, theta * s * s);
    let (h01, h11) = (theta * theta * (3.0 - 2.0 * theta), -theta * theta * s);
    for i in 0..out.len() {
        out[i] = h00 * y0[i] + h10 * h * f0[i] + h01 * y1[i] + h11 * h * f1[i];
    }
}

/// The continuous extension of order 4 of a Dormand–Prince step (DD-12, FR-31):
/// writes the state at `t + θ h` from the step's start `y0`, end `y1` and stages
/// `k[0..7]` (`k[6] = f(t + h, y1)`). Coefficients `D1`–`D7` and the form of
/// `CONTD5` from Hairer's `dopri5.f`
/// (<https://www.unige.ch/~hairer/prog/nonstiff/dopri5.f>).
pub fn dopri5_dense(theta: f64, h: f64, y0: &[f64], y1: &[f64], k: &[Vec<f64>], out: &mut [f64]) {
    const D1: f64 = -12715105075.0 / 11282082432.0;
    const D3: f64 = 87487479700.0 / 32700410799.0;
    const D4: f64 = -10690763975.0 / 1880347072.0;
    const D5: f64 = 701980252875.0 / 199316789632.0;
    const D6: f64 = -1453857185.0 / 822651844.0;
    const D7: f64 = 69997945.0 / 29380423.0;
    let s = 1.0 - theta;
    for i in 0..out.len() {
        let ydiff = y1[i] - y0[i];
        let bspl = h * k[0][i] - ydiff;
        let c3 = ydiff - h * k[6][i] - bspl;
        let c4 = h
            * (D1 * k[0][i]
                + D3 * k[2][i]
                + D4 * k[3][i]
                + D5 * k[4][i]
                + D6 * k[5][i]
                + D7 * k[6][i]);
        out[i] = y0[i] + theta * (ydiff + s * (bspl + theta * (c3 + s * c4)));
    }
}

/// A root of `g` on `[0, 1]` by the Illinois method (DD-12, FR-31), given
/// `g0 = g(0) ≠ 0` and `g1 = g(1)` of the other sign or zero. Returns a point where
/// `g` was evaluated to exactly zero, or the end of a bracket at most `4ε` wide on
/// the side of `g1`. For a continuous `g`, that bracket end is never before the root
/// it encloses; for a discontinuous one, see `Event` (issue #42). When two
/// steps in a row do not halve the bracket, the next step bisects: the bracket
/// halves at least every three steps, so it always closes within the 200
/// evaluations allowed (PR #41).
pub fn illinois<G: FnMut(f64) -> f64 + ?Sized>(g: &mut G, g0: f64, g1: f64) -> f64 {
    // The weights wa, wb start as g(a), g(b) and are halved by the Illinois rule, so
    // they may underflow to zero: only an evaluated zero ends the search (PR #41).
    let (mut a, mut wa, mut b, mut wb) = (0.0, g0, 1.0, g1);
    let mut zero = g1 == 0.0;
    // Which end was kept last: halve the other end's weight when the same end stays twice.
    let mut kept = 0;
    // The bracket width at the last halving, and the steps since then.
    let (mut halved_at, mut stalled) = (1.0, 0);
    for _ in 0..200 {
        if zero || b - a <= 4.0 * f64::EPSILON {
            break;
        }
        let mut c = (a * wb - b * wa) / (wb - wa);
        if stalled >= 2 || !(c > a && c < b) {
            c = 0.5 * (a + b);
        }
        let gc = g(c);
        if gc == 0.0 || gc.signum() == g1.signum() {
            (b, wb, zero) = (c, gc, gc == 0.0);
            if kept == -1 {
                wa *= 0.5;
            }
            kept = -1;
        } else {
            (a, wa) = (c, gc);
            if kept == 1 {
                wb *= 0.5;
            }
            kept = 1;
        }
        if b - a <= 0.5 * halved_at {
            (halved_at, stalled) = (b - a, 0);
        } else {
            stalled += 1;
        }
    }
    b
}

/// Dormand–Prince 5(4) coefficients (Dormand & Prince, J. Comput. Appl. Math. 6(1),
/// 19–26, 1980). `A[i]` holds row `i + 1` of the Butcher matrix.
pub mod dopri5 {
    /// Nodes `c_i`.
    pub const C: [f64; 7] = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];
    /// Butcher matrix rows 2..=7.
    pub const A: [[f64; 6]; 6] = [
        [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
        [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
        [
            19372.0 / 6561.0,
            -25360.0 / 2187.0,
            64448.0 / 6561.0,
            -212.0 / 729.0,
            0.0,
            0.0,
        ],
        [
            9017.0 / 3168.0,
            -355.0 / 33.0,
            46732.0 / 5247.0,
            49.0 / 176.0,
            -5103.0 / 18656.0,
            0.0,
        ],
        [
            35.0 / 384.0,
            0.0,
            500.0 / 1113.0,
            125.0 / 192.0,
            -2187.0 / 6784.0,
            11.0 / 84.0,
        ],
    ];
    /// Fifth-order weights (equal to the last row: first same as last).
    pub const B: [f64; 7] = [
        35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0,
        0.0,
    ];
    /// Fourth-order (embedded) weights.
    pub const BHAT: [f64; 7] = [
        5179.0 / 57600.0,
        0.0,
        7571.0 / 16695.0,
        393.0 / 640.0,
        -92097.0 / 339200.0,
        187.0 / 2100.0,
        1.0 / 40.0,
    ];
}

/// One Dormand–Prince 5(4) step. `w.k[0]` must hold `f(t, y)` on entry (first same
/// as last); on return `w.k[6]` holds `f(t + h, y_new)`. Writes the error estimate
/// `y5 − y4` into `err`.
pub fn dopri5_step<F: FnMut(f64, &[f64], &mut [f64]) + ?Sized>(
    f: &mut F,
    t: f64,
    y: &[f64],
    h: f64,
    w: &mut Work,
    y_new: &mut [f64],
    err: &mut [f64],
) {
    use dopri5::{A, B, BHAT, C};
    let n = y.len();
    for s in 1..7 {
        for i in 0..n {
            let mut acc = y[i];
            for (j, a) in A[s - 1].iter().enumerate().take(s) {
                acc += h * a * w.k[j][i];
            }
            w.tmp[i] = acc;
        }
        let (_, rest) = w.k.split_at_mut(s);
        f(t + C[s] * h, &w.tmp, &mut rest[0]);
    }
    for i in 0..n {
        let mut y5 = y[i];
        let mut e = 0.0;
        for s in 0..7 {
            y5 += h * B[s] * w.k[s][i];
            e += h * (B[s] - BHAT[s]) * w.k[s][i];
        }
        y_new[i] = y5;
        err[i] = e;
    }
}

#[cfg(test)]
mod tests {
    //! TC-U23 (order of Dormand–Prince with fixed steps; checks the coefficients).
    use super::*;

    /// Global error at t = 1 for y' = y cos t (exact y = exp(sin t)) with `steps` steps.
    fn dopri_error(steps: usize) -> f64 {
        let mut f = |t: f64, y: &[f64], out: &mut [f64]| out[0] = y[0] * t.cos();
        let h = 1.0 / steps as f64;
        let mut w = Work::new(7, 1);
        let (mut y, mut y_new, mut err) = (vec![1.0], vec![0.0], vec![0.0]);
        let mut t = 0.0;
        f(t, &y, &mut w.k[0]);
        for _ in 0..steps {
            dopri5_step(&mut f, t, &y, h, &mut w, &mut y_new, &mut err);
            let last = w.k[6].clone();
            w.k[0].copy_from_slice(&last);
            t += h;
            core::mem::swap(&mut y, &mut y_new);
        }
        (y[0] - 1f64.sin().exp()).abs()
    }

    #[test]
    fn illinois_always_returns_a_closed_bracket_or_a_zero() {
        // TC-U29: the condition itself, for g monotone near the root: g(b) = 0, or
        // g(b) has the sign of g(1) and g(b − 4ε) the other sign. PR #41 reviews 4 and
        // 5 (Codex): a ratio of 1e30 between the ends, and a weight that underflows.
        type Case = (&'static str, fn(f64) -> f64);
        let cases: [Case; 6] = [
            ("cubic", |x| x * x * x - 0.2),
            ("sqrt", |x| x.sqrt() - 0.5),
            ("x^100", |x| x.powi(100) - 0.5f64.powi(100)),
            ("exp", |x| (50.0 * x).exp() - 2.0),
            ("falling", |x| 0.3 - x),
            ("subnormal", |x| {
                (x - 0.9) * (-7420.0 * (x - 0.9).max(0.0)).exp()
            }),
        ];
        for (name, g) in cases {
            let (g0, g1) = (g(0.0), g(1.0));
            let b = illinois(&mut |x: f64| g(x), g0, g1);
            let gb = g(b);
            let closed = gb == 0.0
                || (gb.signum() == g1.signum()
                    && g(b - 4.0 * f64::EPSILON).signum() == g0.signum());
            assert!(closed, "{name}: b = {b}, g(b) = {gb:e}");
        }
    }

    #[test]
    fn dopri5_dense_output_is_fourth_order() {
        // TC-U29: the error at the middle of one step of y' = y cos t from t = 0 is
        // O(h⁵) (order 4), and the ends are the step's own values.
        let mut f = |t: f64, y: &[f64], out: &mut [f64]| out[0] = y[0] * t.cos();
        let mut mid_error = |h: f64| {
            let mut w = Work::new(7, 1);
            let (y, mut y_new, mut err, mut out) = ([1.0], [0.0], [0.0], [0.0]);
            f(0.0, &y, &mut w.k[0]);
            dopri5_step(&mut f, 0.0, &y, h, &mut w, &mut y_new, &mut err);
            dopri5_dense(0.0, h, &y, &y_new, &w.k, &mut out);
            assert_eq!(out[0], y[0]);
            dopri5_dense(1.0, h, &y, &y_new, &w.k, &mut out);
            assert!((out[0] - y_new[0]).abs() <= 2.0 * f64::EPSILON);
            dopri5_dense(0.5, h, &y, &y_new, &w.k, &mut out);
            (out[0] - (0.5 * h).sin().exp()).abs()
        };
        let order = (mid_error(0.2) / mid_error(0.1)).log2();
        assert!((order - 5.0).abs() <= 0.3, "local order {order}");
    }

    #[test]
    fn hermite_is_exact_for_cubics() {
        // TC-U29: y = t³ − 2t on [0.5, 2], with y' = 3t² − 2 at the ends.
        let (a, h) = (0.5, 1.5);
        let y = |t: f64| t * t * t - 2.0 * t;
        let dy = |t: f64| 3.0 * t * t - 2.0;
        let mut out = [0.0];
        for theta in [0.0, 0.1, 0.5, 0.9, 1.0] {
            hermite(
                theta,
                h,
                &[y(a)],
                &[dy(a)],
                &[y(a + h)],
                &[dy(a + h)],
                &mut out,
            );
            assert!((out[0] - y(a + theta * h)).abs() < 1e-14, "θ = {theta}");
        }
    }

    #[test]
    fn illinois_returns_the_root_from_the_far_side() {
        // TC-U29: root of θ³ − 0.2 at 0.2^(1/3); the result is never before it.
        let mut calls = 0;
        let mut g = |x: f64| {
            calls += 1;
            x * x * x - 0.2
        };
        let root = illinois(&mut g, -0.2, 0.8);
        let exact = 0.2f64.cbrt();
        assert!(
            root >= exact && root - exact <= 4.0 * f64::EPSILON,
            "{root} vs {exact}"
        );
        assert!(calls < 20, "{calls} evaluations");
        // A concave g keeps the other end; without halving, either case takes over 50.
        let mut calls = 0;
        let mut g = |x: f64| {
            calls += 1;
            x.sqrt() - 0.5
        };
        let root = illinois(&mut g, -0.5, 0.5);
        assert!(root >= 0.25 && root - 0.25 <= 4.0 * f64::EPSILON, "{root}");
        assert!(calls < 20, "{calls} evaluations");
        // PR #41 fourth review (Codex): g(1)/|g(0)| ≈ 1e30 stalls the Illinois weights;
        // the bracket must still close on the root.
        let c = 0.5f64.powi(100);
        let mut g = |x: f64| x.powi(100) - c;
        let root = illinois(&mut g, -c, 1.0 - c);
        assert!(root >= 0.5 && root - 0.5 <= 4.0 * f64::EPSILON, "{root}");
        // Falling, and a zero at the end of the step.
        let mut g = |x: f64| 0.3 - x;
        let root = illinois(&mut g, 0.3, -0.7);
        assert!(root >= 0.3 && root - 0.3 <= 4.0 * f64::EPSILON, "{root}");
        assert_eq!(illinois(&mut |x: f64| x - 1.0, -1.0, 0.0), 1.0);
    }

    #[test]
    fn dormand_prince_is_fifth_order() {
        // Halving h divides the global error by 2^5 = 32. Allowance ±0.3 in the order.
        let order = (dopri_error(10) / dopri_error(20)).log2();
        assert!((order - 5.0).abs() <= 0.3, "order {order}");
    }

    /// Global error at t = 1 for y' = y cos t with `steps` fixed RODAS steps and the
    /// exact Jacobian cos t and ∂f/∂t = −y sin t.
    fn rodas_error(steps: usize) -> f64 {
        let mut f = |t: f64, y: &[f64], out: &mut [f64]| out[0] = y[0] * t.cos();
        let h = 1.0 / steps as f64;
        let mut w = RosenbrockWork::new(1);
        let (mut y, mut y_new, mut err) = (vec![1.0], vec![0.0], vec![0.0]);
        let mut t = 0.0;
        for _ in 0..steps {
            let f0 = [y[0] * f64::cos(t)];
            let jac = [f64::cos(t)];
            let ft = [-y[0] * f64::sin(t)];
            rodas_step(
                &mut f, t, &y, h, &f0, &jac, &ft, &mut w, &mut y_new, &mut err,
            )
            .unwrap();
            t += h;
            core::mem::swap(&mut y, &mut y_new);
        }
        (y[0] - 1f64.sin().exp()).abs()
    }

    #[test]
    fn rodas_is_fourth_order() {
        // TC-U25: halving h divides the global error by 2^4. Allowance ±0.3 (04).
        let order = (rodas_error(10) / rodas_error(20)).log2();
        assert!((order - 4.0).abs() <= 0.3, "order {order}");
    }

    #[test]
    fn first_same_as_last_holds() {
        // The last stage is evaluated at (t + h, y_new): it equals f(t + h, y_new).
        let mut f = |t: f64, y: &[f64], out: &mut [f64]| out[0] = y[0] * t.cos();
        let mut w = Work::new(7, 1);
        let (y, mut y_new, mut err) = (vec![1.0], vec![0.0], vec![0.0]);
        f(0.0, &y, &mut w.k[0]);
        dopri5_step(&mut f, 0.0, &y, 0.1, &mut w, &mut y_new, &mut err);
        let mut direct = vec![0.0];
        f(0.1, &y_new, &mut direct);
        assert!((w.k[6][0] - direct[0]).abs() <= 2.0 * f64::EPSILON);
    }
}

/// RODAS coefficients, order 4(3), stiffly accurate, L-stable. Copied from
/// subroutine `ROCOE` (METH = 1) of Hairer's `rodas.f`
/// (<https://www.unige.ch/~hairer/prog/stiff/rodas.f>), the code accompanying
/// Hairer & Wanner, *Solving Ordinary Differential Equations II*, §IV.7.
#[allow(clippy::excessive_precision)] // digits as published in rodas.f, for comparison with the source
pub mod rodas {
    /// Nodes.
    pub const C2: f64 = 0.386;
    /// Nodes.
    pub const C3: f64 = 0.21;
    /// Nodes.
    pub const C4: f64 = 0.63;
    /// `γ_i` for the time derivative.
    pub const D: [f64; 4] = [0.25, -0.1043, 0.1035, -0.3620000000000023e-01];
    /// Stage weights `a_ij`.
    pub const A21: f64 = 0.1544e+01;
    /// Stage weights.
    pub const A31: f64 = 0.9466785280815826;
    /// Stage weights.
    pub const A32: f64 = 0.2557011698983284;
    /// Stage weights.
    pub const A41: f64 = 0.3314825187068521e+01;
    /// Stage weights.
    pub const A42: f64 = 0.2896124015972201e+01;
    /// Stage weights.
    pub const A43: f64 = 0.9986419139977817;
    /// Stage weights.
    pub const A51: f64 = 0.1221224509226641e+01;
    /// Stage weights.
    pub const A52: f64 = 0.6019134481288629e+01;
    /// Stage weights.
    pub const A53: f64 = 0.1253708332932087e+02;
    /// Stage weights.
    pub const A54: f64 = -0.6878860361058950;
    /// Coupling `c_ij` (divided by `h` in the stages).
    pub const C21: f64 = -0.56688e+01;
    /// Coupling.
    pub const C31: f64 = -0.2430093356833875e+01;
    /// Coupling.
    pub const C32: f64 = -0.2063599157091915;
    /// Coupling.
    pub const C41: f64 = -0.1073529058151375;
    /// Coupling.
    pub const C42: f64 = -0.9594562251023355e+01;
    /// Coupling.
    pub const C43: f64 = -0.2047028614809616e+02;
    /// Coupling.
    pub const C51: f64 = 0.7496443313967647e+01;
    /// Coupling.
    pub const C52: f64 = -0.1024680431464352e+02;
    /// Coupling.
    pub const C53: f64 = -0.3399990352819905e+02;
    /// Coupling.
    pub const C54: f64 = 0.1170890893206160e+02;
    /// Coupling.
    pub const C61: f64 = 0.8083246795921522e+01;
    /// Coupling.
    pub const C62: f64 = -0.7981132988064893e+01;
    /// Coupling.
    pub const C63: f64 = -0.3152159432874371e+02;
    /// Coupling.
    pub const C64: f64 = 0.1631930543123136e+02;
    /// Coupling.
    pub const C65: f64 = -0.6058818238834054e+01;
    /// Diagonal `γ`.
    pub const GAMMA: f64 = 0.25;
}

/// Why the matrix `E = I/(hγ) − J` of a Rosenbrock step could not be factored.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatrixError {
    /// `E` is finite, and both as it is and with its rows scaled, has a zero row or
    /// meets an exactly zero pivot.
    Singular,
    /// `E` is finite, but its elimination overflows both as it is and with its rows
    /// scaled.
    Overflow,
    /// `E` has an entry that is NaN or infinite (from `J`, or from `1/(hγ) − J_ii`).
    NonFinite,
}

/// Work arrays of a Rosenbrock step.
pub struct RosenbrockWork {
    /// `E = I/(hγ) − J` before scaling.
    pub raw: Vec<f64>,
    /// Row scales `D` (powers of two) of `E`.
    pub d: Vec<f64>,
    /// `D E`, LU-factored in place.
    pub e: Vec<f64>,
    /// Row permutation of the factorization.
    pub perm: Vec<usize>,
    /// Stage increments `k_1 … k_6`.
    pub k: Vec<Vec<f64>>,
    /// `f` at a stage point.
    pub dy: Vec<f64>,
    /// Stage point.
    pub ys: Vec<f64>,
}
impl RosenbrockWork {
    /// Allocates once per solve (RES-02).
    pub fn new(n: usize) -> Self {
        RosenbrockWork {
            raw: vec![0.0; n * n],
            d: vec![0.0; n],
            e: vec![0.0; n * n],
            perm: vec![0; n],
            k: vec![vec![0.0; n]; 6],
            dy: vec![0.0; n],
            ys: vec![0.0; n],
        }
    }
}

/// One RODAS step (DD-12), following `ROSCOR` in Hairer's `rodas.f`.
///
/// `f0 = f(t, y)`, `jac` the Jacobian `∂f/∂y` at `(t, y)` (row-major), `ft` the
/// time derivative `∂f/∂t`. Writes the solution into `y_new` and the error estimate
/// (the last stage increment) into `err`. Fails if `E` is not finite, or if LU
/// decomposition meets an exactly zero pivot both on `E` and on `E` with its rows
/// scaled (issue #54).
#[allow(clippy::too_many_arguments)]
pub fn rodas_step<F: FnMut(f64, &[f64], &mut [f64]) + ?Sized>(
    f: &mut F,
    t: f64,
    y: &[f64],
    h: f64,
    f0: &[f64],
    jac: &[f64],
    ft: &[f64],
    w: &mut RosenbrockWork,
    y_new: &mut [f64],
    err: &mut [f64],
) -> Result<(), MatrixError> {
    use rodas::*;
    use tensorcheck::internal::utility::{lu_decompose_with, lu_solve, scale_rows};
    let n = y.len();
    let fac = 1.0 / (h * GAMMA);
    for i in 0..n {
        for j in 0..n {
            w.raw[i * n + j] = -jac[i * n + j] + if i == j { fac } else { 0.0 };
        }
    }
    if !w.raw.iter().all(|x| x.is_finite()) {
        return Err(MatrixError::NonFinite);
    }
    // E as it is first, as in `rodas.f`: every E factored so before issue #54 gives
    // the same factors. Only an exactly zero or non-finite pivot fails here (see
    // `lu_decompose_with`).
    w.e.copy_from_slice(&w.raw);
    w.d.fill(1.0);
    if lu_decompose_with(&mut w.e, n, &mut w.perm, false).is_err() {
        // Then with the rows scaled by powers of two (DD-08), so that the largest
        // entry of each is of order one: rows near f64::MAX no longer overflow in
        // the elimination (issue #54), though growth by partial pivoting still can
        // (2^(n−1) at most, see below). The scaling rounds
        // an entry only if it underflows, by at most 2^-1075 / d_i, far below the
        // rounding of the row (ε max_j |E_ij| ≈ ε / d_i; PR #59 review).
        scale_rows(&w.raw, n, &mut w.d, &mut w.e).map_err(|_| MatrixError::Singular)?;
        if lu_decompose_with(&mut w.e, n, &mut w.perm, false).is_err() {
            // The scaled entries are below about √2, but partial pivoting may still
            // double them at each elimination: Wilkinson's matrix of size 1025
            // overflows (PR #59 review).
            return Err(if w.e.iter().all(|x| x.is_finite()) {
                MatrixError::Singular
            } else {
                MatrixError::Overflow
            });
        }
    }
    let (hd1, hd2, hd3, hd4) = (h * D[0], h * D[1], h * D[2], h * D[3]);
    // E x = b is (D E) x = D b (D = I if E was factored as it is). Scaling b by a
    // power of two is exact unless D b is subnormal; the error is then below
    // 2^-1075 / d_i, far below the rounding of E's row i (about ε / d_i).
    let solve = |w: &mut RosenbrockWork, s: usize| {
        let mut b = core::mem::take(&mut w.k[s]);
        for (x, d) in b.iter_mut().zip(&w.d) {
            *x *= d;
        }
        lu_solve(&w.e, n, &w.perm, &mut b);
        w.k[s] = b;
    };
    // Stage 1.
    for i in 0..n {
        w.k[0][i] = f0[i] + hd1 * ft[i];
    }
    solve(w, 0);
    // Stage 2.
    for i in 0..n {
        w.ys[i] = y[i] + A21 * w.k[0][i];
    }
    f(t + C2 * h, &w.ys, &mut w.dy);
    for i in 0..n {
        w.k[1][i] = w.dy[i] + hd2 * ft[i] + C21 / h * w.k[0][i];
    }
    solve(w, 1);
    // Stage 3.
    for i in 0..n {
        w.ys[i] = y[i] + A31 * w.k[0][i] + A32 * w.k[1][i];
    }
    f(t + C3 * h, &w.ys, &mut w.dy);
    for i in 0..n {
        w.k[2][i] = w.dy[i] + hd3 * ft[i] + (C31 * w.k[0][i] + C32 * w.k[1][i]) / h;
    }
    solve(w, 2);
    // Stage 4.
    for i in 0..n {
        w.ys[i] = y[i] + A41 * w.k[0][i] + A42 * w.k[1][i] + A43 * w.k[2][i];
    }
    f(t + C4 * h, &w.ys, &mut w.dy);
    for i in 0..n {
        w.k[3][i] =
            w.dy[i] + hd4 * ft[i] + (C41 * w.k[0][i] + C42 * w.k[1][i] + C43 * w.k[2][i]) / h;
    }
    solve(w, 3);
    // Stage 5 (autonomous part only, as in ROSCOR).
    for i in 0..n {
        w.ys[i] = y[i] + A51 * w.k[0][i] + A52 * w.k[1][i] + A53 * w.k[2][i] + A54 * w.k[3][i];
    }
    f(t + h, &w.ys, &mut w.dy);
    for i in 0..n {
        w.k[4][i] =
            w.dy[i] + (C51 * w.k[0][i] + C52 * w.k[1][i] + C53 * w.k[2][i] + C54 * w.k[3][i]) / h;
    }
    solve(w, 4);
    // Embedded solution, then stage 6.
    for i in 0..n {
        w.ys[i] += w.k[4][i];
    }
    f(t + h, &w.ys, &mut w.dy);
    for i in 0..n {
        w.k[5][i] = w.dy[i]
            + (C61 * w.k[0][i]
                + C62 * w.k[1][i]
                + C63 * w.k[2][i]
                + C64 * w.k[3][i]
                + C65 * w.k[4][i])
                / h;
    }
    solve(w, 5);
    for i in 0..n {
        y_new[i] = w.ys[i] + w.k[5][i];
        err[i] = w.k[5][i];
    }
    Ok(())
}
