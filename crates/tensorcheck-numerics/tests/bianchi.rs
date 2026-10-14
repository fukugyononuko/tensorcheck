//! TC-U15: differential identities of the curvature (second Bianchi identity and
//! the contracted Bianchi identity), with the Riemann and Einstein tensors
//! differentiated as functions of the point by Ridders' method.
#![allow(clippy::needless_range_loop)]

use tensorcheck::geometry::{Einstein, Riemann, curvature};
use tensorcheck::units::si::{Dimensionless, Length, m};
use tensorcheck::{Metric, Point, Quantity, Scalar, VectorSpace};
use tensorcheck_numerics::diff::{Richardson, derivative};

tensorcheck::space!(G3, 3, Length);
tensorcheck::chart!(C3 on G3 { B1: Length, B2: Length, B3: Length });

/// A generic smooth positive-definite metric tensor in three dimensions.
struct M3;
impl Metric for M3 {
    type Chart = C3;
    type Dim = Dimensionless;
    fn components<N: Scalar>(&self, p: &Point<C3, N>) -> [N; 6] {
        let n = 3;
        core::array::from_fn(|k| {
            let (mut i, mut rem) = (0, k);
            while rem >= n - i {
                rem -= n - i;
                i += 1;
            }
            let j = i + rem;
            let mut h = N::zero();
            for q in 0..n {
                let w = N::from_f64(0.3 + 0.1 * ((i + 2 * j + 3 * q) % 5) as f64);
                h = h + (p.get(q) * w + N::from_f64((i + j) as f64)).sin() * N::from_f64(1.0 / 3.0);
            }
            (if i == j { N::one() } else { N::zero() }) + h * N::from_f64(0.2)
        })
    }
}

const N: usize = 3;
const X0: [f64; 3] = [0.3, -0.7, 1.1];

fn at(x: [f64; 3]) -> Point<C3> {
    Point::<C3>::new((x[0] * m, x[1] * m, x[2] * m))
}

/// ∂_e of every component of the Riemann and Einstein tensors at X0, with the
/// largest error estimates.
fn partials() -> (Vec<Vec<f64>>, Vec<Vec<f64>>, f64, f64) {
    let (mut dr, mut dg) = (Vec::new(), Vec::new());
    let (mut er, mut eg) = (0.0f64, 0.0f64);
    for e in 0..N {
        let shifted = |q: Quantity<Length>| {
            let mut x = X0;
            x[e] = q.value_in(m);
            x
        };
        let r = derivative(
            |q: Quantity<Length>| -> Riemann<C3> {
                curvature(&M3, &at(shifted(q))).unwrap().riemann
            },
            X0[e] * m,
            Richardson,
        );
        let g = derivative(
            |q: Quantity<Length>| -> Einstein<C3> {
                curvature(&M3, &at(shifted(q))).unwrap().einstein
            },
            X0[e] * m,
            Richardson,
        );
        dr.push((0..N.pow(4)).map(|k| r.value.component(k)).collect());
        dg.push((0..N * N).map(|k| g.value.component(k)).collect());
        er = er.max(r.error_estimate);
        eg = eg.max(g.error_estimate);
    }
    (dr, dg, er, eg)
}

/// (largest |second Bianchi residual|, its scale, largest |∇^a G_ab|, its scale,
/// derivative error estimates)
fn residuals() -> (f64, f64, f64, f64, f64, f64) {
    let cv = curvature(&M3, &at(X0)).unwrap();
    let gm: Vec<f64> = (0..9)
        .map(|k| tensorcheck::metric_tensor(&M3, &at(X0)).component(k))
        .collect();
    let mut ginv = vec![0.0; 9];
    tensorcheck::utility::linear::inverse(&gm, 3, &mut ginv).unwrap();
    let gam = |a: usize, b: usize, c: usize| cv.christoffel.component((a * N + b) * N + c);
    let rie = |a: usize, b: usize, c: usize, d: usize| {
        cv.riemann.component(((a * N + b) * N + c) * N + d)
    };
    let ein = |a: usize, b: usize| cv.einstein.component(a * N + b);
    let (dr, dg, er, eg) = partials();
    let d_r =
        |e: usize, a: usize, b: usize, c: usize, d: usize| dr[e][((a * N + b) * N + c) * N + d];
    // ∇_e R^a_bcd
    let nabla_r = |e: usize, a: usize, b: usize, c: usize, d: usize| {
        let mut v = d_r(e, a, b, c, d);
        for f in 0..N {
            v += gam(a, e, f) * rie(f, b, c, d)
                - gam(f, e, b) * rie(a, f, c, d)
                - gam(f, e, c) * rie(a, b, f, d)
                - gam(f, e, d) * rie(a, b, c, f);
        }
        v
    };
    let (mut bianchi, mut scale_r): (f64, f64) = (0.0, 0.0);
    for a in 0..N {
        for b in 0..N {
            for c in 0..N {
                for d in 0..N {
                    for e in 0..N {
                        let s = nabla_r(e, a, b, c, d)
                            + nabla_r(c, a, b, d, e)
                            + nabla_r(d, a, b, e, c);
                        bianchi = bianchi.max(s.abs());
                        scale_r = scale_r.max(nabla_r(e, a, b, c, d).abs());
                    }
                }
            }
        }
    }
    // ∇^a G_ab = g^{ac} (∂_c G_ab − Γ^f_ca G_fb − Γ^f_cb G_af)
    let (mut div, mut scale_g): (f64, f64) = (0.0, 0.0);
    for b in 0..N {
        let mut s = 0.0;
        for a in 0..N {
            for c in 0..N {
                let mut ng = dg[c][a * N + b];
                for f in 0..N {
                    ng -= gam(f, c, a) * ein(f, b) + gam(f, c, b) * ein(a, f);
                }
                s += ginv[a * N + c] * ng;
                scale_g = scale_g.max(ng.abs());
            }
        }
        div = div.max(s.abs());
    }
    (bianchi, scale_r, div, scale_g, er, eg)
}

#[test]
#[ignore]
fn measure() {
    let (b, sr, d, sg, er, eg) = residuals();
    println!(
        "second Bianchi: residual {b:e} (scale {sr:e}), derivative estimate {er:e}, ratio {:.2}",
        b / er
    );
    println!(
        "contracted Bianchi: residual {d:e} (scale {sg:e}), derivative estimate {eg:e}, ratio {:.2}",
        d / eg
    );
}

#[test]
fn second_and_contracted_bianchi_identities() {
    // Measured (2026-10-14): residuals 0.03 and 0.07 × the derivative error
    // estimate (Ridders). Design value (04, TC-U15): 10 × the estimate.
    let (b, _sr, d, _sg, er, eg) = residuals();
    assert!(
        b <= 10.0 * er,
        "second Bianchi residual {b} vs estimate {er}"
    );
    assert!(
        d <= 10.0 * eg,
        "contracted Bianchi residual {d} vs estimate {eg}"
    );
}
