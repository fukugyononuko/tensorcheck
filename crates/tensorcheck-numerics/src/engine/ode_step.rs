//! One step of an ODE method (CMP-E05, DD-12). Plain arrays only:
//! `f(t, y, out)` writes `dy/dt` into `out`.

// Index loops mirror the formulas (several arrays share the index).
#![allow(clippy::needless_range_loop)]

/// The right-hand side on plain arrays.
pub type Rhs<'a> = dyn FnMut(f64, &[f64], &mut [f64]) + 'a;

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
pub fn rk4(f: &mut Rhs, t: f64, y: &[f64], h: f64, w: &mut Work, y_new: &mut [f64]) {
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
pub fn dopri5_step(
    f: &mut Rhs,
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
    fn dormand_prince_is_fifth_order() {
        // Halving h divides the global error by 2^5 = 32. Allowance ±0.3 in the order.
        let order = (dopri_error(10) / dopri_error(20)).log2();
        assert!((order - 5.0).abs() <= 0.3, "order {order}");
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
