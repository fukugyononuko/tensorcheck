//! Error norms and step-size control (CMP-E06, DD-13).

// The step-size functions are used by the ODE solver (implementation phase 7).
#![allow(dead_code)]

/// Hairer's weighted RMS norm of an error estimate:
/// `sqrt(mean((e_i / (atol_i + rtol · max(|y_i|, |y'_i|)))²))`.
pub fn rms_norm(err: &[f64], y_old: &[f64], y_new: &[f64], atol: &[f64], rtol: f64) -> f64 {
    let n = err.len().max(1);
    let s: f64 = err
        .iter()
        .zip(y_old)
        .zip(y_new)
        .zip(atol)
        .map(|(((e, a), b), t)| {
            let sc = t + rtol * a.abs().max(b.abs());
            (e / sc).powi(2)
        })
        .sum();
    (s / n as f64).sqrt()
}

/// The largest component of `e_i / (atol_i + rtol · |v_i|)`.
pub fn max_norm(err: &[f64], value: &[f64], atol: &[f64], rtol: f64) -> f64 {
    err.iter()
        .zip(value)
        .zip(atol)
        .map(|((e, v), t)| e / (t + rtol * v.abs()))
        .fold(0.0, f64::max)
}

/// Accept the step if the scaled error is at most one.
pub fn judge(err: f64) -> bool {
    err <= 1.0
}

/// Bounds on the step-size factor (DD-13).
pub const FAC_MIN: f64 = 0.2;
/// Upper bound on the step-size factor.
pub const FAC_MAX: f64 = 5.0;
/// Safety factor.
pub const SAFETY: f64 = 0.9;

/// The next step size: `h · min(f_max, max(f_min, 0.9 · err^(-1/(q+1))))`, with
/// `f_max = 1` right after a rejection.
pub fn next_step(h: f64, err: f64, q: u32, just_rejected: bool) -> f64 {
    let fmax = if just_rejected { 1.0 } else { FAC_MAX };
    let fac = if err == 0.0 {
        fmax
    } else {
        SAFETY * err.powf(-1.0 / (q as f64 + 1.0))
    };
    h * fac.clamp(FAC_MIN, fmax)
}

/// The interval with the largest scaled error (for adaptive quadrature).
pub fn pick_interval(scaled_errors: &[f64]) -> usize {
    scaled_errors
        .iter()
        .enumerate()
        .fold(
            (0, f64::NEG_INFINITY),
            |acc, (i, &e)| if e > acc.1 { (i, e) } else { acc },
        )
        .0
}

#[cfg(test)]
mod tests {
    //! TC-U20.
    use super::*;

    #[test]
    fn rms_norm_matches_hand_calculation() {
        // e = (1e-6, 2e-6), y = (1, -4) → (2, -3), atol = (1e-6, 1e-6), rtol = 1e-3:
        // scales 1e-6 + 2e-3 and 1e-6 + 4e-3.
        let n = rms_norm(
            &[1e-6, 2e-6],
            &[1.0, -4.0],
            &[2.0, -3.0],
            &[1e-6, 1e-6],
            1e-3,
        );
        let s1: f64 = 1e-6 / (1e-6 + 2e-3);
        let s2 = 2e-6 / (1e-6 + 4e-3);
        assert_eq!(n, ((s1 * s1 + s2 * s2) / 2.0).sqrt());
    }

    #[test]
    fn step_factor_is_bounded() {
        for &err in &[0.0, 1e-12, 0.1, 1.0, 10.0, 1e12] {
            for q in [1, 3, 4] {
                let f = next_step(1.0, err, q, false);
                assert!((FAC_MIN..=FAC_MAX).contains(&f), "{err} {q} {f}");
                // After a rejection the step never grows.
                assert!(next_step(1.0, err, q, true) <= 1.0);
            }
        }
        assert_eq!(next_step(2.0, 1.0, 4, false), 2.0 * SAFETY);
    }

    #[test]
    fn picks_the_interval_with_the_largest_error() {
        assert_eq!(pick_interval(&[0.1, 3.0, 2.0]), 1);
        assert!(judge(1.0) && !judge(1.0000001));
    }
}
