//! Gravity 1: the curvature of Schwarzschild spacetime, from the metric tensor alone.
//!
//! Outside a spherical mass the Ricci tensor vanishes (Einstein's equations in
//! vacuum, `R_ab = 0`), yet spacetime is curved: the Kretschmann scalar
//! `K = R_abcd R^abcd` equals `12 r_s² / r⁶`.

mod common;

use common::{Schw, SchwarzschildMetric, gm_sun, schwarzschild_radius, write_csv};
use tensorcheck::geometry::curvature;
use tensorcheck::labels::{a, b, c, d};
use tensorcheck::units::si::{m, rad};
use tensorcheck::{Point, slot};

fn main() {
    let rs = schwarzschild_radius(gm_sun()); // about 2.95 km
    let g = SchwarzschildMetric::new(rs);
    let rs_m = rs.value_in(m);
    println!("Schwarzschild radius of the Sun: {rs_m:.3} m");

    let mut rows = Vec::new();
    let mut worst_ricci: f64 = 0.0;
    let mut worst_k: f64 = 0.0;
    for i in 0..=40 {
        let r = rs_m * 1.2 * 10f64.powf(i as f64 / 20.0); // 1.2 r_s to 120 r_s
        let p = Point::<Schw>::new((0.0 * m, r * m, 1.1 * rad, 0.4 * rad));
        let cv = curvature(&g, &p).unwrap();

        // K = R_abcd R^abcd: lower the first index, raise the other three.
        let riem = cv.riemann;
        let low = riem.lower(slot::_0, &g, &p);
        let up = riem
            .raise(slot::_1, &g, &p)
            .and_then(|t| t.raise(slot::_2, &g, &p))
            .and_then(|t| t.raise(slot::_3, &g, &p))
            .unwrap();
        let k = (low.idx((a, b, c, d)) * up.idx((a, b, c, d))).into_tensor();
        let k = k.value_in(m.pow::<-4>());
        let want = 12.0 * rs_m * rs_m / r.powi(6);
        worst_k = worst_k.max((k - want).abs() / want);

        // R_ab R^ab, compared with K: both are invariants, so the ratio does not
        // depend on the chart.
        let ric_up = cv
            .ricci
            .raise(slot::_0, &g, &p)
            .and_then(|t| t.raise(slot::_1, &g, &p))
            .unwrap();
        let ric2 = (cv.ricci.idx((a, b)) * ric_up.idx((a, b))).into_tensor();
        worst_ricci = worst_ricci.max((ric2.value_in(m.pow::<-4>()).abs() / k).sqrt());
        rows.push(vec![r / rs_m, k * rs_m.powi(4), want * rs_m.powi(4)]);
    }
    println!("max sqrt(|R_ab R^ab| / K): {worst_ricci:e}");
    println!("max relative error of K against 12 r_s²/r⁶: {worst_k:e}");
    write_csv("r_over_rs,k_rs4,k_exact_rs4", &rows);
    assert!(worst_ricci < 1e-10);
    assert!(worst_k < 1e-10);
}
