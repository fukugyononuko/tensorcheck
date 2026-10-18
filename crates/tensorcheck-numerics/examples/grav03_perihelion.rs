//! Gravity 3: the perihelion precession of Mercury from the geodesic equation of
//! Schwarzschild spacetime, `d²x^a/dτ² = −Γ^a_bc u^b u^c`.
//!
//! The Christoffel symbols come from the metric tensor by automatic
//! differentiation at every step; nothing about orbits is put in by hand. The
//! precession per orbit is compared with `6πGM / (c² a (1 − e²))`, which gives the
//! famous 43″ per century.

mod common;

use common::{GeoState, Phi, R, SchwarzschildMetric, c, gm_sun, schwarzschild_radius, write_csv};
use std::f64::consts::PI;
use tensorcheck::units::si::{Dimensionless, Length, Time, Velocity, m, one, rad, s};
use tensorcheck::{Point, Quantity, Vector};
use tensorcheck_numerics::ode::{DormandPrince54, Event, Tolerance, solve_ivp_events_with};
use tensorcheck_numerics::quad::{self, GaussKronrod, integrate};

/// The state at perihelion `r_p` of the bound orbit with aphelion `r_a`, in the
/// equatorial plane.
fn at_perihelion(rs: f64, r_p: f64, r_a: f64) -> GeoState {
    // Turning points of (dr/dτ)² = c²[E² − (1 − r_s/r)(1 + L²/r²)] fix L (here a
    // length: L = r² dφ/(c dτ)).
    let (up, ua) = (1.0 / r_p, 1.0 / r_a);
    let l2 = rs * (ua - up) / ((ua * ua - rs * ua.powi(3)) - (up * up - rs * up.powi(3)));
    let cc = c().value_in(m / s);
    let u_phi = cc * l2.sqrt() / (r_p * r_p);
    // g(u, u) = c² fixes u^t.
    let u_t = ((cc * cc + r_p * r_p * u_phi * u_phi) / (1.0 - rs / r_p)).sqrt();
    (
        Point::new((0.0 * m, r_p * m, 0.5 * PI * rad, 0.0 * rad)),
        Vector::from_parts((
            u_t * (m / s),
            0.0 * (m / s),
            0.0 * (rad / s),
            u_phi * (rad / s),
        )),
    )
}

/// The exact advance per orbit, `2 ∫ (dφ/dr) dr − 2π` from `r_p` to `r_a`, by
/// quadrature. The substitution `r = r_p + (r_a − r_p)(1 − cos ψ)/2` removes the
/// inverse square roots at the turning points.
fn advance_by_quadrature(rs: f64, r_p: f64, r_a: f64) -> f64 {
    let (up, ua) = (1.0 / r_p, 1.0 / r_a);
    let l2 = rs * (ua - up) / ((ua * ua - rs * ua.powi(3)) - (up * up - rs * up.powi(3)));
    let e2 = (1.0 - rs * up) * (1.0 + l2 * up * up);
    let integrand = |psi: Quantity<Dimensionless>| -> Quantity<Dimensionless> {
        let psi = psi.value_in(one);
        let r = r_p + 0.5 * (r_a - r_p) * (1.0 - psi.cos());
        let dr = 0.5 * (r_a - r_p) * psi.sin();
        let v = e2 - (1.0 - rs / r) * (1.0 + l2 / (r * r));
        // Gauss–Kronrod nodes are inside the interval, so v > 0 at every node.
        let x = dr / v.sqrt();
        (l2.sqrt() / (r * r) * x) * one
    };
    let half = integrate(
        integrand,
        (0.0 * one)..(PI * one),
        GaussKronrod,
        quad::Tolerance {
            rtol: 1e-12,
            atol: 1e-14 * one,
        },
    )
    .unwrap()
    .value
    .value_in(one);
    2.0 * half - 2.0 * PI
}

fn tolerance(r_scale: f64) -> Tolerance<GeoState> {
    Tolerance {
        rtol: 1e-12,
        atol: (
            Point::new((
                1e-3 * r_scale * m,
                1e-9 * r_scale * m,
                1e-14 * rad,
                1e-14 * rad,
            )),
            Vector::from_parts((
                1e-6 * (m / s),
                1e-9 * (m / s),
                1e-24 * (rad / s),
                1e-24 * (rad / s),
            )),
        ),
    }
}

/// Follows the orbit for `tau_end` seconds of proper time. Returns the angle φ of
/// every perihelion passage after the start, and the path as (x/r_s, y/r_s).
fn perihelia(rs: f64, y0: GeoState, tau_end: f64) -> (Vec<f64>, Vec<Vec<f64>>) {
    let g = SchwarzschildMetric::new(rs * m);
    let tol = tolerance(y0.0.coord(R));
    let mut path = Vec::new();
    // dr/dτ goes from negative to positive at a perihelion.
    let perihelion = Event::new(|_tau: Quantity<Time>, y: &GeoState| -> Quantity<Velocity> {
        y.1.base_components()[1] * (m / s)
    })
    .rising();
    let sol = solve_ivp_events_with(
        geodesic_rhs!(&g),
        (0.0 * s)..(tau_end * s),
        y0,
        DormandPrince54,
        tol,
        perihelion,
        |_, y| {
            let (r, f) = (y.0.coord(R), y.0.coord(Phi));
            path.push(vec![r * f.cos() / rs, r * f.sin() / rs]);
        },
    )
    .unwrap();
    let found = sol.events.iter().map(|e| e.y.0.coord(Phi)).collect();
    (found, path)
}

fn main() {
    // Mercury: semi-major axis, eccentricity and orbital period.
    let a: Quantity<Length> = 5.790_905e10 * m;
    let e = 0.205_630;
    let period_days = 87.969_1;
    let rs = schwarzschild_radius(gm_sun()).value_in(m);
    let (r_p, r_a) = ((1.0 - e) * a.value_in(m), (1.0 + e) * a.value_in(m));

    // One orbit and a bit more, in proper time.
    let (phis, _) = perihelia(
        rs,
        at_perihelion(rs, r_p, r_a),
        1.05 * period_days * 86_400.0,
    );
    let shift = phis[0] - 2.0 * PI;
    let want = 6.0 * PI * gm_sun().base_value()
        / (c().value_in(m / s).powi(2) * a.value_in(m) * (1.0 - e * e));
    let arcsec = 180.0 / PI * 3600.0;
    let per_century = |x: f64| x * arcsec * 36_525.0 / period_days;
    println!(
        "Mercury: perihelion moves {:.4e} rad per orbit (formula {:.4e}), {:.2}″ per century",
        shift,
        want,
        per_century(shift)
    );
    assert!((shift - want).abs() < 1e-3 * want);
    assert!((per_century(shift) - 42.98).abs() < 0.01);

    // The same equation near a compact object (r_p = 10 r_s, r_a = 30 r_s): the
    // ellipse turns into a rosette.
    let (r_p, r_a) = (10.0 * rs, 30.0 * rs);
    let y0 = at_perihelion(rs, r_p, r_a);
    let orbit_time = 2.0 * PI * (0.5 * (r_p + r_a)).powf(1.5) / gm_sun().base_value().sqrt();
    let (phis, path) = perihelia(rs, y0, 6.0 * orbit_time);
    let exact = advance_by_quadrature(rs, r_p, r_a);
    let mut prev = 0.0;
    for (i, p) in phis.iter().enumerate() {
        let adv = p - prev - 2.0 * PI;
        println!(
            "compact object: perihelion {} at φ = {:.6} rad, advance {:.9} rad",
            i + 1,
            p,
            adv
        );
        assert!((adv - exact).abs() < 1e-8);
        prev = *p;
    }
    println!("advance per orbit by quadrature: {exact:.9} rad");
    assert!(phis.len() >= 5);
    write_csv("x_rs,y_rs", &path);
}
