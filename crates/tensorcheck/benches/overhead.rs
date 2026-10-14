//! TC-N01: run time of checked expressions versus the same computation written by
//! hand on `f64` arrays. Run with `cargo bench -p tensorcheck --bench overhead`.
//! Prints the ratio (checked / hand-written); the best of several rounds is used.

use std::hint::black_box;
use std::time::Instant;

use tensorcheck::labels::{a, b, c, d};
use tensorcheck::spaces::minkowski::{Inertial, MinkowskiMetric};
use tensorcheck::units::si::{Dimensionless, Length, m, one};
use tensorcheck::{Down, Point, Rank2, Rank4, Up, Vector, slot};

fn best_of<F: FnMut() -> f64>(rounds: usize, mut f: F) -> f64 {
    (0..rounds).map(|_| f()).fold(f64::INFINITY, f64::min)
}

fn time<F: FnMut()>(n: usize, mut f: F) -> f64 {
    let t = Instant::now();
    for _ in 0..n {
        f();
    }
    t.elapsed().as_secs_f64()
}

fn main() {
    let n = 2_000_000;
    let g4 = [
        [1.0, 0.1, 0.2, 0.3],
        [0.1, 2.0, 0.4, 0.5],
        [0.2, 0.4, 3.0, 0.6],
        [0.3, 0.5, 0.6, 4.0],
    ];
    let v4 = [1.0, 2.0, 3.0, 4.0];
    let g: Rank2<Dimensionless, Down<Inertial>, Down<Inertial>> = Rank2::from_components(g4, one);
    let v: Vector<Length, Inertial> = Vector::from_components(v4, m);

    // 1. g_ab V^b (rank 2 × rank 1).
    let checked = best_of(5, || {
        time(n, || {
            let w = (black_box(g).idx((a, b)) * black_box(v).idx(b)).into_tensor();
            black_box(w);
        })
    });
    let hand = best_of(5, || {
        time(n, || {
            let (gg, vv) = (black_box(g4), black_box(v4));
            let mut w = [0.0; 4];
            for i in 0..4 {
                for j in 0..4 {
                    w[i] += gg[i][j] * vv[j];
                }
            }
            black_box(w);
        })
    });
    println!(
        "contraction g_ab V^b (4D):       ratio {:.2}",
        checked / hand
    );

    // 2. R^a_bcd T^cd (rank 4 × rank 2).
    let r4: [[[[f64; 4]; 4]; 4]; 4] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            core::array::from_fn(|k| {
                core::array::from_fn(|l| ((i * 64 + j * 16 + k * 4 + l) % 7) as f64 - 3.0)
            })
        })
    });
    let t2 = [[0.5; 4]; 4];
    let r: Rank4<Dimensionless, Up<Inertial>, Down<Inertial>, Down<Inertial>, Down<Inertial>> =
        Rank4::from_components(r4, one);
    let t: Rank2<Dimensionless, Up<Inertial>, Up<Inertial>> = Rank2::from_components(t2, one);
    let m2 = n / 20;
    let checked = best_of(5, || {
        time(m2, || {
            let w = (black_box(r).idx((a, b, c, d)) * black_box(t).idx((c, d))).into_tensor();
            black_box(w);
        })
    });
    let hand = best_of(5, || {
        time(m2, || {
            let (rr, tt) = (black_box(r4), black_box(t2));
            let mut w = [[0.0; 4]; 4];
            for i in 0..4 {
                for j in 0..4 {
                    for k in 0..4 {
                        for l in 0..4 {
                            w[i][j] += rr[i][j][k][l] * tt[k][l];
                        }
                    }
                }
            }
            black_box(w);
        })
    });
    println!(
        "contraction R^a_bcd T^cd (4D):   ratio {:.2}",
        checked / hand
    );

    // 3. Lowering an index with the Minkowski metric tensor.
    let x = Point::<Inertial>::new((0.0 * m, 0.0 * m, 0.0 * m, 0.0 * m));
    let checked = best_of(5, || {
        time(n / 10, || {
            black_box(black_box(v).lower(slot::_0, &MinkowskiMetric, &x));
        })
    });
    let hand = best_of(5, || {
        time(n / 10, || {
            // The same work as the library: build η as a full matrix, then multiply.
            let vv = black_box(v4);
            let eta = black_box([
                [1.0, 0.0, 0.0, 0.0],
                [0.0, -1.0, 0.0, 0.0],
                [0.0, 0.0, -1.0, 0.0],
                [0.0, 0.0, 0.0, -1.0],
            ]);
            let w: [f64; 4] = core::array::from_fn(|i| (0..4).map(|k| eta[i][k] * vv[k]).sum());
            black_box(w);
        })
    });
    println!(
        "lower with η (4D):               ratio {:.2}",
        checked / hand
    );
}
