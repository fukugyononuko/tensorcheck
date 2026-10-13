//! Phase 5: connection, curvature, covariant derivatives.
//! TC-U10, TC-U11, TC-U12, TC-U14, TC-U16.
#![allow(clippy::needless_range_loop)]

use tensorcheck::geometry::{TensorField, christoffel, covariant_derivative, curvature};
use tensorcheck::spaces::minkowski::Spacetime;
use tensorcheck::spaces::plane::{Polar, PolarMetric};
use tensorcheck::spaces::sphere2::{
    Angles, AnglesToStereographic, SphereMetric, StereographicMetric,
};
use tensorcheck::units::si::{Dimensionless, Length, m, rad};
use tensorcheck::utility::typelevel::{Cons, Nil};
use tensorcheck::{ChartTransform, Down, Metric, NoSym, Point, Scalar, Tensor, Tensorial, Up};

const EPS: f64 = f64::EPSILON;

fn flat2<const N: usize>(x: [[f64; N]; N]) -> Vec<f64> {
    x.iter().flatten().copied().collect()
}
fn flat3<const N: usize>(x: [[[f64; N]; N]; N]) -> Vec<f64> {
    x.iter().flatten().flatten().copied().collect()
}
fn flat4<const N: usize>(x: [[[[f64; N]; N]; N]; N]) -> Vec<f64> {
    x.iter().flatten().flatten().flatten().copied().collect()
}
fn max_abs(v: &[f64]) -> f64 {
    v.iter().fold(0.0, |acc, x| acc.max(x.abs()))
}

// ---------------------------------------------------------------------------
// Metric tensors used by the tests (all user-defined, signature −+++).

tensorcheck::chart!(SchwarzschildChart on Spacetime { CT: Length, Rr: Length, Th: Dimensionless, Ph: Dimensionless });
tensorcheck::chart!(Comoving on Spacetime { T0: Length, X1: Length, X2: Length, X3: Length });

/// Schwarzschild metric tensor with Schwarzschild radius `rs` (metres).
struct Schwarzschild {
    rs: f64,
}
impl Metric for Schwarzschild {
    type Chart = SchwarzschildChart;
    type Dim = Dimensionless;
    fn components<N: Scalar>(&self, p: &Point<SchwarzschildChart, N>) -> [N; 10] {
        let (r, th) = (p.get(1), p.get(2));
        let f = N::one() - N::from_f64(self.rs) / r;
        let z = N::zero();
        let s = th.sin();
        [-f, z, z, z, N::one() / f, z, z, r * r, z, r * r * s * s]
    }
}

/// Kerr metric tensor in Boyer–Lindquist coordinates (mass `mm` and spin `aa` in metres).
struct Kerr {
    mm: f64,
    aa: f64,
}
impl Metric for Kerr {
    type Chart = SchwarzschildChart;
    type Dim = Dimensionless;
    fn components<N: Scalar>(&self, p: &Point<SchwarzschildChart, N>) -> [N; 10] {
        let (r, th) = (p.get(1), p.get(2));
        let (mm, aa) = (N::from_f64(self.mm), N::from_f64(self.aa));
        let (s, c) = (th.sin(), th.cos());
        let sigma = r * r + aa * aa * c * c;
        let delta = r * r - N::from_f64(2.0) * mm * r + aa * aa;
        let two = N::from_f64(2.0);
        let z = N::zero();
        let gtt = -(N::one() - two * mm * r / sigma);
        let gtp = -(two * mm * aa * r * s * s / sigma);
        let gpp = (r * r + aa * aa + two * mm * aa * aa * r * s * s / sigma) * s * s;
        [gtt, z, z, gtp, sigma / delta, z, z, sigma, z, gpp]
    }
}

/// A spatially flat FLRW metric tensor with scale factor `(1 + 0.3 ct/m)^0.7`.
struct Flrw;
impl Metric for Flrw {
    type Chart = Comoving;
    type Dim = Dimensionless;
    fn components<N: Scalar>(&self, p: &Point<Comoving, N>) -> [N; 10] {
        let a = (N::one() + N::from_f64(0.3) * p.get(0)).powf(0.7);
        let a2 = a * a;
        let z = N::zero();
        [-N::one(), z, z, z, a2, z, z, a2, z, a2]
    }
}

/// A generic smooth positive-definite metric tensor `δ_ab + 0.2 h_ab(x)` in `n`
/// dimensions, with every component depending on every coordinate.
macro_rules! generic_metric {
    ($space:ident, $chart:ident, $metric:ident, $n:literal, $($x:ident),+) => {
        tensorcheck::space!($space, $n, Length);
        tensorcheck::chart!($chart on $space { $($x: Length),+ });
        struct $metric;
        impl Metric for $metric {
            type Chart = $chart;
            type Dim = Dimensionless;
            fn components<N: Scalar>(&self, p: &Point<$chart, N>) -> [N; $n * ($n + 1) / 2] {
                let n = $n;
                core::array::from_fn(|k| {
                    // invert k -> (i, j) with i <= j
                    let (mut i, mut rem) = (0, k);
                    while rem >= n - i {
                        rem -= n - i;
                        i += 1;
                    }
                    let j = i + rem;
                    let mut h = N::zero();
                    for q in 0..n {
                        let w = N::from_f64(0.3 + 0.1 * ((i + 2 * j + 3 * q) % 5) as f64);
                        h = h + (p.get(q) * w + N::from_f64((i + j) as f64)).sin() * N::from_f64(1.0 / n as f64);
                    }
                    let base = if i == j { N::one() } else { N::zero() };
                    base + h * N::from_f64(0.2)
                })
            }
        }
    };
}
generic_metric!(G2, C2, M2, 2, A1, A2);
generic_metric!(G3, C3, M3, 3, B1, B2, B3);
generic_metric!(G4, C4, M4, 4, D1, D2, D3, D4);
generic_metric!(G5, C5, M5, 5, E1, E2, E3, E4, E5);

// ---------------------------------------------------------------------------
// Identity checks on flat arrays: g (n×n), Γ (n³), R^a_bcd (n⁴).

struct Deviations {
    gamma_sym: f64,
    antisym_cd: f64,
    antisym_ab: f64,
    pair: f64,
    bianchi1: f64,
    ricci_sym: f64,
    einstein_trace: f64,
    low_scale: f64,
}

fn identities(
    n: usize,
    g: &[f64],
    gamma: &[f64],
    r_up: &[f64],
    ricci: &[f64],
    scalar: f64,
    ein: &[f64],
) -> Deviations {
    let ri = |a: usize, b: usize, c: usize, d: usize| r_up[((a * n + b) * n + c) * n + d];
    // R_abcd = g_ae R^e_bcd
    let mut low = vec![0.0; n * n * n * n];
    for a in 0..n {
        for b in 0..n {
            for c in 0..n {
                for d in 0..n {
                    low[((a * n + b) * n + c) * n + d] =
                        (0..n).map(|e| g[a * n + e] * ri(e, b, c, d)).sum();
                }
            }
        }
    }
    let rl = |a: usize, b: usize, c: usize, d: usize| low[((a * n + b) * n + c) * n + d];
    let mut dv = Deviations {
        gamma_sym: 0.0,
        antisym_cd: 0.0,
        antisym_ab: 0.0,
        pair: 0.0,
        bianchi1: 0.0,
        ricci_sym: 0.0,
        einstein_trace: 0.0,
        low_scale: max_abs(&low),
    };
    for a in 0..n {
        for b in 0..n {
            for c in 0..n {
                dv.gamma_sym = dv
                    .gamma_sym
                    .max((gamma[(a * n + b) * n + c] - gamma[(a * n + c) * n + b]).abs());
                for d in 0..n {
                    dv.antisym_cd = dv.antisym_cd.max((rl(a, b, c, d) + rl(a, b, d, c)).abs());
                    dv.antisym_ab = dv.antisym_ab.max((rl(a, b, c, d) + rl(b, a, c, d)).abs());
                    dv.pair = dv.pair.max((rl(a, b, c, d) - rl(c, d, a, b)).abs());
                    dv.bianchi1 = dv
                        .bianchi1
                        .max((ri(a, b, c, d) + ri(a, c, d, b) + ri(a, d, b, c)).abs());
                }
            }
        }
    }
    // Ricci straight from the Riemann tensor (before the library symmetrizes it).
    for b in 0..n {
        for d in 0..n {
            let rbd: f64 = (0..n).map(|a| ri(a, b, a, d)).sum();
            let rdb: f64 = (0..n).map(|a| ri(a, d, a, b)).sum();
            dv.ricci_sym = dv.ricci_sym.max((rbd - rdb).abs());
            let _ = ricci;
        }
    }
    // g^ab G_ab = (1 − n/2) R
    let mut ginv = vec![0.0; n * n];
    tensorcheck::utility::linear::inverse(g, n, &mut ginv).unwrap();
    let tr: f64 = (0..n * n).map(|k| ginv[k] * ein[k]).sum();
    dv.einstein_trace = (tr - (1.0 - n as f64 / 2.0) * scalar).abs();
    dv
}

macro_rules! check_metric {
    ($name:expr, $metric:expr, $point:expr, $n:literal, $report:expr) => {{
        let g = $metric;
        let p = $point;
        let cv = curvature(&g, &p).unwrap();
        let gm = flat2(tensorcheck::metric_tensor(&g, &p).base_components());
        let gamma = flat3(cv.christoffel.base_components());
        let r = flat4(cv.riemann.base_components());
        let ric = flat2(cv.ricci.base_components());
        let ein = flat2(cv.einstein.base_components());
        let scale = max_abs(&r).max(1e-300);
        let dv = identities($n, &gm, &gamma, &r, &ric, cv.ricci_scalar.base_value(), &ein);
        if std::env::var("TC_DETAIL").is_ok() {
            let low_scale = dv.low_scale.max(1e-300);
            println!(
                "{:14} cd {:7.2} ab {:7.2} pair {:7.2} (÷max|R_abcd|)  bianchi {:7.2} ricci {:7.2} trace {:7.2} (÷max|R^a_bcd|)",
                $name,
                dv.antisym_cd / low_scale / EPS,
                dv.antisym_ab / low_scale / EPS,
                dv.pair / low_scale / EPS,
                dv.bianchi1 / scale / EPS,
                dv.ricci_sym / scale / EPS,
                dv.einstein_trace / scale / EPS
            );
        }
        // Identities of R_abcd are measured against max |R_abcd|, those of R^a_bcd
        // against max |R^a_bcd| (lowering multiplies components by g, e.g. r² ≈ 49).
        let low_scale = dv.low_scale.max(1e-300);
        let worst = [
            dv.antisym_cd / low_scale,
            dv.antisym_ab / low_scale,
            dv.pair / low_scale,
            dv.bianchi1 / scale,
            dv.ricci_sym / scale,
            dv.einstein_trace / scale,
        ];
        let gscale = max_abs(&gamma).max(1e-300);
        $report($name, dv.gamma_sym / gscale / EPS, worst.iter().fold(0.0f64, |acc, x| acc.max(*x)) / EPS);
    }};
}

fn all_identity_cases(report: &mut dyn FnMut(&str, f64, f64)) {
    let r = |x: f64| x * m;
    check_metric!(
        "polar",
        PolarMetric,
        Point::<Polar>::new((r(1.7), 0.4 * rad)),
        2,
        report
    );
    check_metric!(
        "sphere",
        SphereMetric::new(r(2.0)),
        Point::<Angles>::new((0.9 * rad, 0.3 * rad)),
        2,
        report
    );
    check_metric!(
        "schwarzschild",
        Schwarzschild { rs: 2.0 },
        Point::<SchwarzschildChart>::new((r(0.0), r(7.0), 1.1 * rad, 0.2 * rad)),
        4,
        report
    );
    check_metric!(
        "kerr",
        Kerr { mm: 1.0, aa: 0.6 },
        Point::<SchwarzschildChart>::new((r(0.0), r(6.0), 0.8 * rad, 0.3 * rad)),
        4,
        report
    );
    check_metric!(
        "flrw",
        Flrw,
        Point::<Comoving>::new((r(1.5), r(0.2), r(-0.4), r(0.7))),
        4,
        report
    );
    check_metric!(
        "generic2",
        M2,
        Point::<C2>::new((r(0.3), r(-0.7))),
        2,
        report
    );
    check_metric!(
        "generic3",
        M3,
        Point::<C3>::new((r(0.3), r(-0.7), r(1.1))),
        3,
        report
    );
    check_metric!(
        "generic4",
        M4,
        Point::<C4>::new((r(0.3), r(-0.7), r(1.1), r(0.2))),
        4,
        report
    );
    check_metric!(
        "generic5",
        M5,
        Point::<C5>::new((r(0.3), r(-0.7), r(1.1), r(0.2), r(-1.3))),
        5,
        report
    );
}

#[test]
#[ignore]
fn measure_identities() {
    all_identity_cases(&mut |name, gamma_sym, worst| {
        println!(
            "{name:14} Γ symmetry {gamma_sym:8.2} ε   Riemann identities {worst:8.2} ε (relative to max |R|)"
        );
    });
}

/// TC-U14 thresholds. Premise: the Levi-Civita connection in a coordinate basis.
/// Design value (04): 100 ε relative to the largest component of the same tensor.
/// Measured (2026-10-14, run with TC_DETAIL=1): worst 5.86 ε (generic 4D, trace
/// of the Einstein tensor); all others ≤ 3.2 ε. The Christoffel symmetry is exact
/// because the formula is symmetric in b, c (measured 0).
const IDENTITY_TOL_EPS: f64 = 100.0;

#[test]
fn riemann_identities_hold_for_special_and_generic_metrics() {
    let mut failures = Vec::new();
    all_identity_cases(&mut |name, gamma_sym, worst| {
        if gamma_sym > 1.0 || worst > IDENTITY_TOL_EPS {
            failures.push(format!("{name}: Γ {gamma_sym} ε, Riemann {worst} ε"));
        }
    });
    assert!(failures.is_empty(), "{failures:?}");
}

// ---------------------------------------------------------------------------
// TC-U10, TC-U11: values against analytic results.

#[test]
fn christoffel_symbols_of_polar_coordinates_and_the_sphere() {
    let p = Point::<Polar>::new((1.7 * m, 0.4 * rad));
    let g = flat3(christoffel(&PolarMetric, &p).unwrap().base_components());
    let n = 2;
    let at = |a: usize, b: usize, c: usize| g[(a * n + b) * n + c];
    // Γ^r_θθ = −r, Γ^θ_rθ = Γ^θ_θr = 1/r, others 0.
    assert!((at(0, 1, 1) + 1.7).abs() <= 10.0 * EPS * 1.7);
    assert!((at(1, 0, 1) - 1.0 / 1.7).abs() <= 10.0 * EPS);
    assert!((at(1, 1, 0) - 1.0 / 1.7).abs() <= 10.0 * EPS);
    for (a, b, c) in [(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0), (1, 1, 1)] {
        assert!(at(a, b, c).abs() <= 10.0 * EPS * 1.7, "{a}{b}{c}");
    }
    let th: f64 = 0.9;
    let q = Point::<Angles>::new((th * rad, 0.3 * rad));
    let gs = flat3(
        christoffel(&SphereMetric::new(2.0 * m), &q)
            .unwrap()
            .base_components(),
    );
    let at = |a: usize, b: usize, c: usize| gs[(a * n + b) * n + c];
    assert!((at(0, 1, 1) + th.sin() * th.cos()).abs() <= 10.0 * EPS);
    assert!((at(1, 0, 1) - th.cos() / th.sin()).abs() <= 10.0 * EPS);
}

#[test]
fn curvature_of_flat_polar_coordinates_and_the_sphere() {
    // Flat plane in polar coordinates: Riemann = 0.
    let p = Point::<Polar>::new((1.7 * m, 0.4 * rad));
    let cv = curvature(&PolarMetric, &p).unwrap();
    assert!(max_abs(&flat4(cv.riemann.base_components())) <= 100.0 * EPS);
    // Sphere of radius a: R = 2/a², R_ab = g_ab/a², G_ab = 0 in two dimensions.
    let a = 2.0;
    let q = Point::<Angles>::new((0.9 * rad, 0.3 * rad));
    let g = SphereMetric::new(a * m);
    let cv = curvature(&g, &q).unwrap();
    let want = 2.0 / (a * a);
    assert!((cv.ricci_scalar.value_in(m.pow::<-2>()) - want).abs() <= 100.0 * EPS * want);
    let gm = flat2(tensorcheck::metric_tensor(&g, &q).base_components());
    let ric = flat2(cv.ricci.base_components());
    for k in 0..4 {
        assert!((ric[k] - gm[k] / (a * a)).abs() <= 100.0 * EPS * gm[0].abs().max(1.0) / (a * a));
    }
    assert!(max_abs(&flat2(cv.einstein.base_components())) <= 100.0 * EPS);
}

#[test]
fn schwarzschild_and_kerr_are_ricci_flat() {
    for (name, ric, scale) in [
        {
            let p = Point::<SchwarzschildChart>::new((0.0 * m, 7.0 * m, 1.1 * rad, 0.2 * rad));
            let cv = curvature(&Schwarzschild { rs: 2.0 }, &p).unwrap();
            (
                "schwarzschild",
                flat2(cv.ricci.base_components()),
                max_abs(&flat4(cv.riemann.base_components())),
            )
        },
        {
            let p = Point::<SchwarzschildChart>::new((0.0 * m, 6.0 * m, 0.8 * rad, 0.3 * rad));
            let cv = curvature(&Kerr { mm: 1.0, aa: 0.6 }, &p).unwrap();
            (
                "kerr",
                flat2(cv.ricci.base_components()),
                max_abs(&flat4(cv.riemann.base_components())),
            )
        },
    ] {
        assert!(
            max_abs(&ric) <= 100.0 * EPS * scale,
            "{name}: {} vs scale {scale}",
            max_abs(&ric)
        );
    }
}

// ---------------------------------------------------------------------------
// TC-U12: covariant derivatives.

/// The metric tensor of the plane in polar coordinates, as a tensor field.
struct PolarMetricField;
impl TensorField for PolarMetricField {
    type Chart = Polar;
    type Dim = Dimensionless;
    type Slots = Cons<Down<Polar>, Cons<Down<Polar>, Nil>>;
    fn at<N: Scalar>(
        &self,
        p: &Point<Polar, N>,
    ) -> Tensor<Dimensionless, Self::Slots, NoSym, Tensorial, N> {
        let r = p.get(0);
        Tensor::from_base_components([[N::one(), N::zero()], [N::zero(), r * r]])
    }
}

/// The constant Cartesian field ∂_x written in polar coordinates:
/// V^r = cos θ, V^θ = −sin θ / r.
struct UnitX;
impl TensorField for UnitX {
    type Chart = Polar;
    type Dim = Length;
    type Slots = Cons<Up<Polar>, Nil>;
    fn at<N: Scalar>(
        &self,
        p: &Point<Polar, N>,
    ) -> Tensor<Length, Self::Slots, NoSym, Tensorial, N> {
        let (r, t) = (p.get(0), p.get(1));
        Tensor::from_base_components([t.cos(), -(t.sin() / r)])
    }
}

#[test]
fn metric_compatibility_and_parallel_fields() {
    let p = Point::<Polar>::new((1.7 * m, 0.4 * rad));
    let dg = covariant_derivative(&PolarMetricField, &PolarMetric, &p).unwrap();
    assert!(max_abs(&flat3(dg.base_components())) <= 100.0 * EPS * 1.7 * 1.7);
    let dv = covariant_derivative(&UnitX, &PolarMetric, &p).unwrap();
    assert!(max_abs(&flat2(dv.base_components())) <= 100.0 * EPS);
}

// ---------------------------------------------------------------------------
// TC-U16: independence of the chart.

#[test]
fn curvature_does_not_depend_on_the_chart() {
    let a = 2.0;
    let q = Point::<Angles>::new((0.9 * rad, 0.3 * rad));
    let xf = AnglesToStereographic::new(a * m);
    let s = xf.map(&q);
    let cv_angles = curvature(&SphereMetric::new(a * m), &q).unwrap();
    let cv_stereo = curvature(&StereographicMetric::new(a * m), &s).unwrap();
    let ra = cv_angles.ricci_scalar.value_in(m.pow::<-2>());
    let rs = cv_stereo.ricci_scalar.value_in(m.pow::<-2>());
    // Measured 0.0 and 2 ε; design value 1000 ε (04 TC-U16).
    assert!((ra - rs).abs() <= 1000.0 * EPS * ra.abs(), "{ra} vs {rs}");
    // The Ricci tensor transformed from angles equals the one computed directly.
    let moved = flat2(
        cv_angles
            .ricci
            .transform(&xf, &q)
            .unwrap()
            .base_components(),
    );
    let direct = flat2(cv_stereo.ricci.base_components());
    let scale = max_abs(&direct);
    for k in 0..4 {
        assert!(
            (moved[k] - direct[k]).abs() <= 1000.0 * EPS * scale,
            "{k}: {} vs {}",
            moved[k],
            direct[k]
        );
    }
}
