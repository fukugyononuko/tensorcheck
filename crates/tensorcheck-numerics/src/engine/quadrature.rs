//! Quadrature rules (CMP-E07, DD-14). Integrands are vector-valued: `f(x, out)`
//! writes the `m` components of `f(x)` into `out`.

/// Composite trapezoid rule with `n` intervals.
pub fn trapezoid(
    f: &mut dyn FnMut(f64, &mut [f64]),
    a: f64,
    b: f64,
    n: usize,
    m: usize,
) -> Vec<f64> {
    let h = (b - a) / n as f64;
    let mut sum = vec![0.0; m];
    let mut v = vec![0.0; m];
    for i in 0..=n {
        f(a + h * i as f64, &mut v);
        let w = if i == 0 || i == n { 0.5 } else { 1.0 };
        for k in 0..m {
            sum[k] += w * v[k];
        }
    }
    sum.iter().map(|s| s * h).collect()
}

/// Composite Simpson rule with `n` (even) intervals.
pub fn simpson(f: &mut dyn FnMut(f64, &mut [f64]), a: f64, b: f64, n: usize, m: usize) -> Vec<f64> {
    let n = n + n % 2;
    let h = (b - a) / n as f64;
    let mut sum = vec![0.0; m];
    let mut v = vec![0.0; m];
    for i in 0..=n {
        f(a + h * i as f64, &mut v);
        let w = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        for k in 0..m {
            sum[k] += w * v[k];
        }
    }
    sum.iter().map(|s| s * h / 3.0).collect()
}

/// Romberg integration: trapezoid with halving step and Richardson extrapolation.
/// Returns `(estimate, |R(k,k) − R(k−1,k−1)|, levels used)`; stops when `done`
/// accepts the pair or after `max_levels`.
pub fn romberg(
    f: &mut dyn FnMut(f64, &mut [f64]),
    a: f64,
    b: f64,
    m: usize,
    max_levels: usize,
    done: &mut dyn FnMut(&[f64], &[f64]) -> bool,
) -> (Vec<f64>, Vec<f64>, usize) {
    let mut v = vec![0.0; m];
    let mut prev: Vec<Vec<f64>> = Vec::new();
    let mut h = b - a;
    // R(0,0)
    let mut t = vec![0.0; m];
    f(a, &mut v);
    t.iter_mut().zip(&v).for_each(|(s, x)| *s += 0.5 * x);
    f(b, &mut v);
    t.iter_mut().zip(&v).for_each(|(s, x)| *s += 0.5 * x);
    let mut row = vec![t.iter().map(|s| s * h).collect::<Vec<f64>>()];
    let mut err = vec![f64::INFINITY; m];
    for level in 1..max_levels {
        // Trapezoid with half the step: add the midpoints.
        let pts = 1usize << (level - 1);
        let mut mid = vec![0.0; m];
        for i in 0..pts {
            f(a + h * (i as f64 + 0.5), &mut v);
            mid.iter_mut().zip(&v).for_each(|(s, x)| *s += x);
        }
        h *= 0.5;
        prev.push(row[row.len() - 1].clone());
        let mut new_row = vec![
            row[0]
                .iter()
                .zip(&mid)
                .map(|(r, s)| 0.5 * r + h * s)
                .collect::<Vec<f64>>(),
        ];
        let mut factor = 1.0;
        for j in 1..=level {
            factor *= 4.0;
            let r: Vec<f64> = new_row[j - 1]
                .iter()
                .zip(&row[j - 1])
                .map(|(x, y)| x + (x - y) / (factor - 1.0))
                .collect();
            new_row.push(r);
        }
        let best = new_row[level].clone();
        err = best
            .iter()
            .zip(&prev[level - 1])
            .map(|(x, y)| (x - y).abs())
            .collect();
        row = new_row;
        if level >= 2 && done(&best, &err) {
            return (best, err, level + 1);
        }
    }
    (row[row.len() - 1].clone(), err, max_levels)
}

/// Nodes and weights of the `n`-point Gauss–Legendre rule on `[-1, 1]`, by Newton's
/// method on the Legendre polynomial `P_n` (DD-14).
pub fn gauss_legendre_nodes(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    for i in 0..n.div_ceil(2) {
        // P_n(z) and P_n'(z) by the three-term recurrence.
        let legendre = |z: f64| {
            let (mut p1, mut p2) = (1.0, 0.0);
            for j in 1..=n {
                let p3 = p2;
                p2 = p1;
                p1 = ((2 * j - 1) as f64 * z * p2 - (j - 1) as f64 * p3) / j as f64;
            }
            (p1, n as f64 * (z * p1 - p2) / (z * z - 1.0))
        };
        let mut z = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        for _ in 0..100 {
            let (p, dp) = legendre(z);
            let z1 = z;
            z = z1 - p / dp;
            if (z - z1).abs() <= 3.0 * f64::EPSILON {
                break;
            }
        }
        // The weight needs P_n' at the converged root, not at the previous iterate.
        let (_, dp) = legendre(z);
        x[i] = -z;
        x[n - 1 - i] = z;
        let wi = 2.0 / ((1.0 - z * z) * dp * dp);
        w[i] = wi;
        w[n - 1 - i] = wi;
    }
    (x, w)
}

/// The `n`-point Gauss–Legendre rule on `[a, b]`.
pub fn gauss_legendre(
    f: &mut dyn FnMut(f64, &mut [f64]),
    a: f64,
    b: f64,
    n: usize,
    m: usize,
) -> Vec<f64> {
    let (x, w) = gauss_legendre_nodes(n);
    let (c, r) = (0.5 * (a + b), 0.5 * (b - a));
    let mut sum = vec![0.0; m];
    let mut v = vec![0.0; m];
    for (xi, wi) in x.iter().zip(&w) {
        f(c + r * xi, &mut v);
        for k in 0..m {
            sum[k] += wi * v[k];
        }
    }
    sum.iter().map(|s| s * r).collect()
}

/// Kronrod nodes of the 15-point rule on `[0, 1]` (QUADPACK `qk15`); odd positions
/// are the 7-point Gauss nodes. Verified by exactness tests (TC-U21).
#[allow(clippy::excessive_precision)] // digits as published, for comparison with the source
pub const XGK: [f64; 8] = [
    0.991_455_371_120_812_639_206_854_697_526_329,
    0.949_107_912_342_758_524_526_189_684_047_851,
    0.864_864_423_359_769_072_789_712_788_640_926,
    0.741_531_185_599_394_439_863_864_773_280_788,
    0.586_087_235_467_691_130_294_144_845_693_013,
    0.405_845_151_377_397_166_906_606_412_076_961,
    0.207_784_955_007_898_467_600_689_403_773_245,
    0.0,
];
/// Kronrod weights (QUADPACK `qk15`).
#[allow(clippy::excessive_precision)] // digits as published, for comparison with the source
pub const WGK: [f64; 8] = [
    0.022_935_322_010_529_224_963_732_008_058_970,
    0.063_092_092_629_978_553_290_700_663_189_204,
    0.104_790_010_322_250_183_839_876_322_541_518,
    0.140_653_259_715_525_918_745_189_590_510_238,
    0.169_004_726_639_267_902_826_583_426_598_550,
    0.190_350_578_064_785_409_913_256_402_421_014,
    0.204_432_940_075_298_892_414_161_999_234_649,
    0.209_482_141_084_727_828_012_999_174_891_714,
];
/// Weights of the 7-point Gauss rule at `XGK[1], XGK[3], XGK[5], XGK[7]`.
#[allow(clippy::excessive_precision)] // digits as published, for comparison with the source
pub const WG: [f64; 4] = [
    0.129_484_966_168_869_693_270_611_432_679_082,
    0.279_705_391_489_276_667_901_467_771_423_780,
    0.381_830_050_505_118_944_950_369_775_488_975,
    0.417_959_183_673_469_387_755_102_040_816_327,
];

/// One 15-point Gauss–Kronrod estimate on `[a, b]`, with QUADPACK's error estimate
/// per component. Returns `(kronrod result, error estimate)`.
pub fn gauss_kronrod15(
    f: &mut dyn FnMut(f64, &mut [f64]),
    a: f64,
    b: f64,
    m: usize,
) -> (Vec<f64>, Vec<f64>) {
    let (c, h) = (0.5 * (a + b), 0.5 * (b - a));
    let mut fv = vec![vec![0.0; m]; 15];
    // fv[0] = center, then pairs (c − h x, c + h x) for j = 0..7
    f(c, &mut fv[0]);
    for j in 0..7 {
        let (lo, hi) = fv.split_at_mut(2 * j + 2);
        f(c - h * XGK[j], &mut lo[2 * j + 1]);
        f(c + h * XGK[j], &mut hi[0]);
    }
    let mut out = vec![0.0; m];
    let mut err = vec![0.0; m];
    for k in 0..m {
        let fc = fv[0][k];
        let mut resk = WGK[7] * fc;
        let mut resg = WG[3] * fc;
        let mut resabs = WGK[7] * fc.abs();
        for j in 0..7 {
            let (f1, f2) = (fv[2 * j + 1][k], fv[2 * j + 2][k]);
            resk += WGK[j] * (f1 + f2);
            resabs += WGK[j] * (f1.abs() + f2.abs());
            if j % 2 == 1 {
                resg += WG[j / 2] * (f1 + f2);
            }
        }
        let reskh = 0.5 * resk;
        let mut resasc = WGK[7] * (fc - reskh).abs();
        for j in 0..7 {
            resasc +=
                WGK[j] * ((fv[2 * j + 1][k] - reskh).abs() + (fv[2 * j + 2][k] - reskh).abs());
        }
        let result = resk * h;
        let (resabs, resasc) = (resabs * h.abs(), resasc * h.abs());
        let mut e = ((resk - resg) * h).abs();
        if resasc != 0.0 && e != 0.0 {
            e = resasc * (200.0 * e / resasc).powf(1.5).min(1.0);
        }
        if resabs > f64::MIN_POSITIVE / (50.0 * f64::EPSILON) {
            e = e.max(50.0 * f64::EPSILON * resabs);
        }
        out[k] = result;
        err[k] = e;
    }
    (out, err)
}

#[cfg(test)]
mod tests {
    //! TC-U21: exactness and convergence orders; verification of the QUADPACK constants.
    use super::*;

    fn poly(k: i32) -> impl FnMut(f64, &mut [f64]) {
        move |x, out| out[0] = x.powi(k)
    }
    fn exact(k: i32, a: f64, b: f64) -> f64 {
        (b.powi(k + 1) - a.powi(k + 1)) / (k + 1) as f64
    }
    fn rel(x: f64, y: f64) -> f64 {
        (x - y).abs() / y.abs().max(1e-300)
    }

    #[test]
    fn rules_are_exact_for_their_polynomial_degree() {
        let (a, b) = (-0.7, 1.9);
        // Trapezoid: degree 1. Simpson: degree 3. Allowance: 10 n ε relative (n sums
        // of rounded terms).
        for k in 0..=1 {
            assert!(
                rel(trapezoid(&mut poly(k), a, b, 8, 1)[0], exact(k, a, b)) <= 80.0 * f64::EPSILON
            );
        }
        for k in 0..=3 {
            assert!(
                rel(simpson(&mut poly(k), a, b, 8, 1)[0], exact(k, a, b)) <= 80.0 * f64::EPSILON
            );
        }
        // Gauss–Legendre, degree 2n − 1: 10 n ε relative (n rounded terms).
        // Measured (2026-10-14, `measure::gauss_legendre_errors`): worst 15.5 ε at
        // n = 16. The first measurement gave 208 ε; the cause was the weights, which
        // used P_n' at the iterate before the last Newton step (fixed).
        for n in 1..=20 {
            for k in 0..(2 * n as i32) {
                let got = gauss_legendre(&mut poly(k), a, b, n, 1)[0];
                let tol = 10.0 * n as f64 * f64::EPSILON;
                assert!(rel(got, exact(k, a, b)) <= tol, "n={n} k={k}");
            }
        }
    }

    #[test]
    fn convergence_orders_of_composite_rules() {
        // ∫_0^1 e^x: halving h divides the error by 4 (trapezoid) and 16 (Simpson).
        let mut f = |x: f64, out: &mut [f64]| out[0] = x.exp();
        let exact = 1f64.exp() - 1.0;
        let order = |e1: f64, e2: f64| (e1 / e2).log2();
        let t = (trapezoid(&mut f, 0.0, 1.0, 16, 1)[0] - exact).abs();
        let t2 = (trapezoid(&mut f, 0.0, 1.0, 32, 1)[0] - exact).abs();
        assert!((order(t, t2) - 2.0).abs() <= 0.2);
        let s = (simpson(&mut f, 0.0, 1.0, 8, 1)[0] - exact).abs();
        let s2 = (simpson(&mut f, 0.0, 1.0, 16, 1)[0] - exact).abs();
        assert!((order(s, s2) - 4.0).abs() <= 0.2);
    }

    #[test]
    fn kronrod_constants_are_correct() {
        // The 7-point Gauss part matches nodes and weights computed independently.
        let (x, w) = gauss_legendre_nodes(7);
        for (j, g) in [(1usize, 0usize), (3, 1), (5, 2), (7, 3)] {
            let xi = XGK[j];
            let k = x.iter().position(|v| (v - xi).abs() < 1e-12).expect("node");
            assert!((x[k] - xi).abs() <= 4.0 * f64::EPSILON);
            assert!(rel(w[k], WG[g]) <= 8.0 * f64::EPSILON);
        }
        // Weights sum to 2 on [-1, 1]; the 15-point rule is exact up to degree 22.
        let sum: f64 = 2.0 * WGK[..7].iter().sum::<f64>() + WGK[7];
        assert!((sum - 2.0).abs() <= 4.0 * f64::EPSILON);
        for k in 0..=22 {
            let (r, _) = gauss_kronrod15(&mut poly(k), -1.0, 1.0, 1);
            let want = exact(k, -1.0, 1.0);
            assert!(
                (r[0] - want).abs() <= 20.0 * f64::EPSILON * want.abs().max(1.0),
                "k={k}: {} vs {want}",
                r[0]
            );
        }
        // Degree 24 is no longer exact (sanity check that the test can fail).
        let (r, _) = gauss_kronrod15(&mut poly(24), -1.0, 1.0, 1);
        assert!((r[0] - exact(24, -1.0, 1.0)).abs() > 1e-12);
    }
}

#[cfg(test)]
mod measure {
    use super::*;
    #[test]
    #[ignore]
    fn gauss_legendre_errors() {
        let (a, b) = (-0.7f64, 1.9f64);
        for n in [4usize, 8, 10, 12, 16, 20] {
            let mut worst: f64 = 0.0;
            let mut worst_abs_scaled: f64 = 0.0;
            for k in 0..(2 * n as i32) {
                let got =
                    gauss_legendre(&mut |x: f64, o: &mut [f64]| o[0] = x.powi(k), a, b, n, 1)[0];
                let want = (b.powi(k + 1) - a.powi(k + 1)) / (k + 1) as f64;
                // Scale by ∫|x^k| so cancellation between the two ends does not inflate the ratio.
                let absint = (b.powi(k + 1) + a.abs().powi(k + 1)) / (k + 1) as f64;
                worst = worst.max((got - want).abs() / want.abs());
                worst_abs_scaled = worst_abs_scaled.max((got - want).abs() / absint);
            }
            let (x, _) = gauss_legendre_nodes(n);
            println!(
                "n={n:2}: rel {:8.1} ε   rel to ∫|f| {:6.2} ε   x0 = {:.17}",
                worst / f64::EPSILON,
                worst_abs_scaled / f64::EPSILON,
                x[0]
            );
        }
    }
}

#[cfg(test)]
mod measure_weights {
    use super::*;
    #[test]
    #[ignore]
    fn weight_sums() {
        for n in 1..=8 {
            let (x, w) = gauss_legendre_nodes(n);
            let s: f64 = w.iter().sum();
            println!(
                "n={n}: Σw − 2 = {:+.2} ε  nodes {:?} weights {:?}",
                (s - 2.0) / f64::EPSILON,
                x,
                w
            );
        }
    }
}
