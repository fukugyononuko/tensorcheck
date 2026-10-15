//! Phase 3: metric tensors, raising/lowering, chart transforms, inverses.
//! TC-U07, TC-U08, TC-I02, TC-I04, TC-R04, TC-R05.

use tensorcheck::labels::{a, b};
use tensorcheck::spaces::euclid3::spherical::{SphericalMetric, SphericalToCartesian};
use tensorcheck::spaces::euclid3::{Cartesian, EuclidMetric};
use tensorcheck::spaces::minkowski::{Inertial, MinkowskiMetric};
use tensorcheck::spaces::plane::{
    Cartesian2, Cartesian2ToPolar, EuclidMetric2, Polar, PolarMetric, PolarToCartesian2,
};
use tensorcheck::units::si::{Area, Dimensionless, Length, Velocity, m, one, rad, s};
use tensorcheck::utility::linear::{LinearError, inverse};
use tensorcheck::{ChartTransform, Covector, Point, Quantity, Vector, slot};

fn rel(x: f64, y: f64) -> f64 {
    (x - y).abs() / y.abs().max(f64::MIN_POSITIVE)
}

/// Deterministic pseudo-random matrix entries in [-1, 1).
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
fn norm_inf(mx: &[f64], n: usize) -> f64 {
    (0..n)
        .map(|i| (0..n).map(|j| mx[i * n + j].abs()).sum::<f64>())
        .fold(0.0, f64::max)
}

#[test]
fn inverse_residual_is_bounded_by_condition_number() {
    // TC-U07. Theory: ‖A A⁻¹ − I‖∞ ≲ c n κ(A) ε. Measured over 2000 random matrices
    // per size (tests/measure_inverse.rs): the worst ratio was 0.50 (n = 1), 0.41,
    // 0.25, 0.31 (n = 4, adjugate), 0.13, 0.11 (n = 6, LU). We allow c = 10.
    for n in 1..=6 {
        for seed in 1..=500 {
            let am = mat(n, seed);
            let mut inv = vec![0.0; n * n];
            if inverse(&am, n, &mut inv).is_err() {
                continue;
            }
            let kappa = norm_inf(&am, n) * norm_inf(&inv, n);
            let mut worst: f64 = 0.0;
            for i in 0..n {
                let mut row = 0.0;
                for j in 0..n {
                    let sum: f64 = (0..n).map(|k| am[i * n + k] * inv[k * n + j]).sum();
                    row += (sum - if i == j { 1.0 } else { 0.0 }).abs();
                }
                worst = worst.max(row);
            }
            assert!(
                worst <= 10.0 * n as f64 * kappa * f64::EPSILON,
                "n={n} seed={seed}"
            );
        }
    }
}

#[test]
fn singular_matrices_are_rejected() {
    // TC-U07: exactly singular matrices (two equal rows) for every size.
    for n in 2..=6 {
        let mut am = mat(n, 7);
        for j in 0..n {
            am[(n - 1) * n + j] = am[j];
        }
        let mut inv = vec![0.0; n * n];
        assert_eq!(
            inverse(&am, n, &mut inv),
            Err(LinearError::Singular),
            "n={n}"
        );
    }
    // Nearly singular but well inside f64 (condition number ~1e8) is accepted.
    for n in 2..=6 {
        // Last row = first row + 1e-8 × an independent random row.
        let mut am = mat(n, 11);
        let noise = mat(n, 99);
        for j in 0..n {
            am[(n - 1) * n + j] = am[j] + 1e-8 * noise[j];
        }
        let mut inv = vec![0.0; n * n];
        assert!(inverse(&am, n, &mut inv).is_ok(), "n={n}");
    }
}

#[test]
fn inverse_does_not_depend_on_the_scale() {
    // TC-U07: the singularity test and the inverse work at any scale (Codex review
    // 2026-10-14, finding 1). With the old test the determinant and the product of
    // the row norms overflowed (1e200) or underflowed (1e-200), and every one of
    // these well-conditioned matrices was rejected.
    let residual = |am: &[f64], inv: &[f64], n: usize| {
        let mut worst: f64 = 0.0;
        for i in 0..n {
            let mut row = 0.0;
            for j in 0..n {
                // A A⁻¹ is of order one whatever the scale of A, so this cannot overflow.
                let sum: f64 = (0..n).map(|k| am[i * n + k] * inv[k * n + j]).sum();
                row += (sum - if i == j { 1.0 } else { 0.0 }).abs();
            }
            worst = worst.max(row);
        }
        worst
    };
    for n in 1..=6 {
        for &scale in &[1e200, 1e-200, 1e300, 1e-300] {
            // A diagonal matrix: condition number 1, so only a few roundings (det, a/det).
            let mut am = vec![0.0; n * n];
            for i in 0..n {
                am[i * n + i] = scale;
            }
            let mut inv = vec![0.0; n * n];
            inverse(&am, n, &mut inv).unwrap_or_else(|e| panic!("n={n} scale={scale}: {e}"));
            for i in 0..n {
                let r = rel(inv[i * n + i], 1.0 / scale);
                assert!(r <= 4.0 * f64::EPSILON, "n={n} scale={scale} rel={r}");
            }
            // A general well-conditioned matrix, scaled.
            let base = mat(n, 3);
            let am: Vec<f64> = base.iter().map(|x| x * scale).collect();
            let mut inv = vec![0.0; n * n];
            inverse(&am, n, &mut inv).unwrap_or_else(|e| panic!("n={n} scale={scale}: {e}"));
            let kappa = norm_inf(&am, n) * norm_inf(&inv, n);
            let r = residual(&am, &inv, n);
            assert!(
                r <= 10.0 * n as f64 * kappa * f64::EPSILON,
                "n={n} scale={scale} residual={r}"
            );
            // Exactly singular stays singular at any scale.
            if n >= 2 {
                let mut am = am.clone();
                for j in 0..n {
                    am[(n - 1) * n + j] = am[j];
                }
                assert_eq!(
                    inverse(&am, n, &mut inv),
                    Err(LinearError::Singular),
                    "n={n} scale={scale}"
                );
            }
        }
        // Rows of very different scales (the Hadamard ratio does not see row scales).
        let row_scale = |i: usize| if i.is_multiple_of(2) { 1e150 } else { 1e-150 };
        let base = mat(n, 5);
        let am: Vec<f64> = base
            .iter()
            .enumerate()
            .map(|(k, x)| x * row_scale(k / n))
            .collect();
        let mut inv = vec![0.0; n * n];
        inverse(&am, n, &mut inv).unwrap_or_else(|e| panic!("n={n} mixed rows: {e}"));
        // A A⁻¹ − I is not scale free here (off-diagonal terms are 1e300 times their
        // difference), so check the equilibrated problem: (D A)(A⁻¹ D⁻¹) − I with
        // D A = base, whose condition number is that of `base`.
        let inv_eq: Vec<f64> = inv
            .iter()
            .enumerate()
            .map(|(k, x)| x * row_scale(k % n))
            .collect();
        let mut base_inv = vec![0.0; n * n];
        inverse(&base, n, &mut base_inv).unwrap();
        let kappa = norm_inf(&base, n) * norm_inf(&base_inv, n);
        let r = residual(&base, &inv_eq, n);
        assert!(
            r <= 10.0 * n as f64 * kappa * f64::EPSILON,
            "n={n} mixed rows residual={r}"
        );
    }
}

#[test]
fn polar_cartesian_round_trip() {
    // TC-U08. Error budget: a few trigonometric evaluations and products,
    // each within 1 ulp; we allow 10 ε relative.
    let p = Point::<Polar>::new((2.0 * m, 0.7 * rad));
    let q = PolarToCartesian2.map(&p);
    let back = Cartesian2ToPolar.map(&q);
    assert!(rel(back.get(0), 2.0) <= 10.0 * f64::EPSILON);
    assert!(rel(back.get(1), 0.7) <= 10.0 * f64::EPSILON);
    // J · J⁻¹ = I using the closed-form inverse Jacobian.
    let j = PolarToCartesian2.jacobian(&p);
    let ji = PolarToCartesian2.inverse_jacobian(&p).unwrap();
    for r in 0..2 {
        for c in 0..2 {
            let sum: f64 = (0..2).map(|k| j[r * 2 + k] * ji[k * 2 + c]).sum();
            assert!((sum - if r == c { 1.0 } else { 0.0 }).abs() <= 10.0 * f64::EPSILON);
        }
    }
}

#[test]
fn speed_squared_agrees_between_polar_and_cartesian() {
    // TC-I02 / CC-02: g_ab V^a V^b in polar coordinates equals |v|² in Cartesian.
    let p = Point::<Polar>::new((2.0 * m, 0.5 * rad));
    let v: Vector<Velocity, Polar> = Vector::from_components([1.5, 0.25], m / s);
    let v_low = v.lower(slot::_0, &PolarMetric, &p);
    let v2: Quantity<tensorcheck::Prod<Velocity, Velocity>> =
        (v_low.idx(a) * v.idx(a)).into_tensor();
    // 1.5² + (r θ')² = 2.25 + (2 · 0.25)² = 2.5
    assert!(rel(v2.value_in(m * m / (s * s)), 2.5) <= 4.0 * f64::EPSILON);

    let vc = v.transform(&PolarToCartesian2, &p).unwrap();
    let q = PolarToCartesian2.map(&p);
    let vc_low = vc.lower(slot::_0, &EuclidMetric2, &q);
    let v2c = (vc_low.idx(b) * vc.idx(b)).into_tensor();
    assert!(rel(v2c.value_in(m * m / (s * s)), 2.5) <= 10.0 * f64::EPSILON);
}

#[test]
fn raising_undoes_lowering() {
    let p = Point::<Polar>::new((3.0 * m, 1.1 * rad));
    let v: Vector<Velocity, Polar> = Vector::from_components([0.3, -0.7], m / s);
    let back = v
        .lower(slot::_0, &PolarMetric, &p)
        .raise(slot::_0, &PolarMetric, &p)
        .unwrap();
    let got = back.components_in(m / s);
    assert!(rel(got[0], 0.3) <= 4.0 * f64::EPSILON);
    assert!(rel(got[1], -0.7) <= 4.0 * f64::EPSILON);
}

#[test]
fn transformed_vector_adds_to_cartesian_vector() {
    // TC-I04 / CC-04: (r, θ) = (2, π/2), V = ∂_r → (0, 1) in Cartesian.
    let p = Point::<Polar>::new((2.0 * m, std::f64::consts::FRAC_PI_2 * rad));
    let v: Vector<Length, Polar> = Vector::from_components([1.0, 0.0], m);
    let vc = v.transform(&PolarToCartesian2, &p).unwrap();
    let w: Vector<Length, Cartesian2> = Vector::from_components([1.0, 0.0], m);
    let sum = (vc + w).components_in(m);
    assert!((sum[0] - 1.0).abs() <= 10.0 * f64::EPSILON);
    assert!((sum[1] - 1.0).abs() <= 10.0 * f64::EPSILON);
    // A covector transforms with the inverse Jacobian: df = ∂_θ f dθ with f = θ
    // becomes (-y/r², x/r²) in Cartesian.
    let dtheta: Covector<tensorcheck::Quot<Dimensionless, Length>, Polar> =
        Covector::from_components([0.0, 1.0], one / m);
    let dc = dtheta
        .transform(&PolarToCartesian2, &p)
        .unwrap()
        .components_in(one / m);
    assert!((dc[0] - (-1.0 / 2.0)).abs() <= 10.0 * f64::EPSILON);
    assert!(dc[1].abs() <= 10.0 * f64::EPSILON);
}

#[test]
fn spherical_metric_matches_cartesian_after_transform() {
    // g_spherical = Jᵀ δ J: lowering in spherical equals transforming, lowering in
    // Cartesian, and transforming back.
    let p = Point::<tensorcheck::spaces::euclid3::spherical::Spherical>::new((
        2.0 * m,
        0.9 * rad,
        0.4 * rad,
    ));
    let v: Vector<Velocity, tensorcheck::spaces::euclid3::spherical::Spherical> =
        Vector::from_components([0.5, 0.1, -0.2], m / s);
    let low_sph: [f64; 3] = v.lower(slot::_0, &SphericalMetric, &p).components_in(m / s);
    let vc = v.transform(&SphericalToCartesian, &p).unwrap();
    let q = SphericalToCartesian.map(&p);
    let low_cart = vc.lower(slot::_0, &EuclidMetric, &q);
    // Pull back: T_a(spherical) = J^i_a T_i(Cartesian).
    let jac = SphericalToCartesian.jacobian(&p);
    let lc = low_cart.components_in(m / s);
    for (ai, want) in low_sph.iter().enumerate() {
        let got: f64 = (0..3).map(|i| jac[i * 3 + ai] * lc[i]).sum();
        assert!(
            (got - want).abs() <= 20.0 * f64::EPSILON * want.abs().max(1.0),
            "{ai}"
        );
    }
    let _ = Cartesian;
}

// TC-R04 (VOL-04): a user-defined chart and transform (spherical, again, from scratch).
mod user_chart {
    use tensorcheck::spaces::euclid3::{Cartesian, Euclid3};
    use tensorcheck::units::si::{Dimensionless, Length};
    use tensorcheck::{ChartTransform, Point, Scalar};
    tensorcheck::chart!(pub MySpherical on Euclid3 { Rr: Length, Th: Dimensionless, Ph: Dimensionless });
    pub struct MyToCartesian;
    impl ChartTransform for MyToCartesian {
        type From = MySpherical;
        type To = Cartesian;
        fn map<N: Scalar>(&self, p: &Point<MySpherical, N>) -> Point<Cartesian, N> {
            let (r, t, f) = (p.get(0), p.get(1), p.get(2));
            Point::from_base([r * t.sin() * f.cos(), r * t.sin() * f.sin(), r * t.cos()])
        }
        fn jacobian<N: Scalar>(&self, p: &Point<MySpherical, N>) -> Vec<N> {
            let (r, t, f) = (p.get(0), p.get(1), p.get(2));
            let (st, ct, sf, cf) = (t.sin(), t.cos(), f.sin(), f.cos());
            vec![
                st * cf,
                r * ct * cf,
                -(r * st * sf),
                st * sf,
                r * ct * sf,
                r * st * cf,
                ct,
                -(r * st),
                N::zero(),
            ]
        }
    }
}

#[test]
fn user_defined_chart_transform() {
    use user_chart::{MySpherical, MyToCartesian};
    let p = Point::<MySpherical>::new((1.0 * m, 0.3 * rad, 1.2 * rad));
    let v: Vector<Length, MySpherical> = Vector::from_components([1.0, 0.0, 0.0], m);
    let vc = v.transform(&MyToCartesian, &p).unwrap().components_in(m);
    // ∂_r is the unit radial vector.
    let want = [
        0.3f64.sin() * 1.2f64.cos(),
        0.3f64.sin() * 1.2f64.sin(),
        0.3f64.cos(),
    ];
    for k in 0..3 {
        assert!((vc[k] - want[k]).abs() <= 10.0 * f64::EPSILON);
    }
}

// TC-R05 (VOL-05): a user-defined Minkowski metric with signature (−, +, +, +).
struct MostlyPlus;
impl tensorcheck::Metric for MostlyPlus {
    type Chart = Inertial;
    type Dim = Dimensionless;
    fn components<N: tensorcheck::Scalar>(&self, _p: &Point<Inertial, N>) -> [N; 10] {
        let (o, z) = (N::one(), N::zero());
        [-o, z, z, z, o, z, z, o, z, o]
    }
}

#[test]
fn user_defined_metric_signature() {
    let x = Point::<Inertial>::new((0.0 * m, 0.0 * m, 0.0 * m, 0.0 * m));
    let u: Vector<Length, Inertial> = Vector::from_components([2.0, 1.0, 0.0, 0.0], m);
    let plus = (u.lower(slot::_0, &MostlyPlus, &x).idx(a) * u.idx(a)).into_tensor();
    let minus = (u.lower(slot::_0, &MinkowskiMetric, &x).idx(a) * u.idx(a)).into_tensor();
    let p: Quantity<Area> = plus;
    assert_eq!(p.value_in(m * m), -3.0);
    assert_eq!(minus.value_in(m * m), 3.0);
}
