//! Measurement for TC-U07 (run with --ignored --nocapture); thresholds are set from this.
use tensorcheck::utility::linear::inverse;

fn mat(n: usize, seed: u64) -> Vec<f64> {
    let mut x = seed;
    (0..n * n)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((x >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
        })
        .collect()
}
fn norm_inf(a: &[f64], n: usize) -> f64 {
    (0..n)
        .map(|i| (0..n).map(|j| a[i * n + j].abs()).sum::<f64>())
        .fold(0.0, f64::max)
}
fn residual(a: &[f64], inv: &[f64], n: usize) -> f64 {
    let mut r: f64 = 0.0;
    for i in 0..n {
        let mut row = 0.0;
        for j in 0..n {
            let mut s = 0.0;
            for k in 0..n {
                s += a[i * n + k] * inv[k * n + j];
            }
            row += (s - if i == j { 1.0 } else { 0.0 }).abs();
        }
        r = r.max(row);
    }
    r
}

#[test]
#[ignore]
fn measure() {
    for n in 1..=6 {
        let mut worst: f64 = 0.0;
        for seed in 0..2000 {
            let a = mat(n, seed + 1);
            let mut inv = vec![0.0; n * n];
            if inverse(&a, n, &mut inv).is_err() {
                continue;
            }
            let kappa = norm_inf(&a, n) * norm_inf(&inv, n);
            let ratio = residual(&a, &inv, n) / (n as f64 * kappa * f64::EPSILON);
            worst = worst.max(ratio);
        }
        println!("n={n}: max residual / (n κ ε) = {worst:.3}");
    }
}
