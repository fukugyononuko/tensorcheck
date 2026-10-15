//! Gravity 4: the bending of light by the Sun. A light ray is a null geodesic of
//! Schwarzschild spacetime; we solve the same geodesic equation as for Mercury and
//! read the deflection from the direction of the ray far from the Sun. The result
//! is compared with `4GM / (c² b)`, which is 1.75″ for a ray grazing the Sun, and
//! with the next term of the series, `(15π/4) (GM / c² b)²`.

mod common;

use common::{GeoState, Phi, R, SchwarzschildMetric, c, gm_sun, schwarzschild_radius, write_csv};
use std::f64::consts::PI;
use tensorcheck::units::si::{m, rad, s};
use tensorcheck::{Point, Vector};
use tensorcheck_numerics::ode::{DormandPrince54, Tolerance, solve_ivp_with};

/// Sends a ray from `x = −x0` along +x with impact parameter `b` and follows it
/// until it is `x0` away again. Returns the deflection angle and the path as
/// (x/r_s, y/r_s).
fn shoot(rs: f64, b: f64, x0: f64) -> (f64, Vec<Vec<f64>>) {
    let g = SchwarzschildMetric::new(rs * m);
    let cc = c().value_in(m / s);
    let (r0, phi0) = (x0.hypot(b), PI - (b / x0).atan());
    // Spatial velocity c along +x, written in (r, φ); u^t from g(u, u) = 0.
    let (u_r, u_phi) = (cc * phi0.cos(), -cc * phi0.sin() / r0);
    let f = 1.0 - rs / r0;
    let u_t = ((u_r * u_r / f + r0 * r0 * u_phi * u_phi) / f).sqrt();
    let y0: GeoState = (
        Point::new((0.0 * m, r0 * m, 0.5 * PI * rad, phi0 * rad)),
        Vector::from_parts((
            u_t * (m / s),
            u_r * (m / s),
            0.0 * (rad / s),
            u_phi * (rad / s),
        )),
    );
    let tol = Tolerance {
        rtol: 1e-12,
        atol: (
            Point::new((1e-6 * b * m, 1e-9 * b * m, 1e-15 * rad, 1e-15 * rad)),
            Vector::from_parts((
                1e-6 * (m / s),
                1e-6 * (m / s),
                1e-24 * (rad / s),
                1e-24 * (rad / s),
            )),
        ),
    };
    let mut path = Vec::new();
    let sol = solve_ivp_with(
        geodesic_rhs!(&g),
        (0.0 * s)..((2.0 * r0 / cc) * s),
        y0,
        DormandPrince54,
        tol,
        |_, y| {
            let (r, f) = (y.0.coord(R), y.0.coord(Phi));
            path.push(vec![r * f.cos() / rs, r * f.sin() / rs]);
        },
    )
    .unwrap();
    // Direction of the outgoing ray in the plane: d(x, y) = (ṙ cos φ − r φ̇ sin φ, ṙ sin φ + r φ̇ cos φ).
    let (r, f) = (sol.y.0.coord(R), sol.y.0.coord(Phi));
    let u = sol.y.1.base_components();
    let (dx, dy) = (
        u[1] * f.cos() - r * u[3] * f.sin(),
        u[1] * f.sin() + r * u[3] * f.cos(),
    );
    (-dy.atan2(dx), path)
}

fn main() {
    let gm = gm_sun().base_value();
    let rs = schwarzschild_radius(gm_sun()).value_in(m);
    let r_sun = 6.957e8; // m (IAU nominal solar radius)
    let arcsec = 180.0 / PI * 3600.0;

    let mut worst: f64 = 0.0;
    for k in [1.0, 2.0, 5.0, 10.0] {
        let b = k * r_sun;
        let (defl, _) = shoot(rs, b, 1e4 * b);
        let x = gm / (c().value_in(m / s).powi(2) * b);
        let (first, second) = (4.0 * x, 4.0 * x + 3.75 * PI * x * x);
        println!(
            "b = {k:>4} R☉: deflection {:.5}″ (4GM/c²b = {:.5}″, relative difference {:.2e})",
            defl * arcsec,
            first * arcsec,
            (defl - first) / first
        );
        worst = worst.max((defl - second).abs() / second);
    }
    println!("largest relative difference with the second-order term: {worst:.1e}");
    // The rest is the third-order term and starting the ray at a finite distance.
    assert!(worst < 1e-6);

    // Near a black hole (impact parameters 2.7 r_s to 8 r_s) the same equation
    // bends rays by large angles; below b = (3√3/2) r_s ≈ 2.6 r_s they fall in.
    let mut rows = Vec::new();
    for (i, k) in [2.7, 3.0, 3.5, 4.5, 6.0, 8.0].into_iter().enumerate() {
        let (defl, path) = shoot(rs, k * rs, 25.0 * rs);
        println!("b = {k} r_s: deflection {:.1}°", defl.to_degrees());
        for p in path {
            rows.push(vec![i as f64, k, p[0], p[1]]);
        }
    }
    write_csv("ray,b_rs,x_rs,y_rs", &rows);
}
