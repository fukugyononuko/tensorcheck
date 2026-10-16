//! Gravity 2: Halley's comet around the Sun with Newton's law of gravitation,
//! `d²x/dt² = −GM x / |x|³`, solved in SI units.
//!
//! The acceleration is built from typed quantities: `GM / r³` has dimension
//! `1/s²`, so writing `GM / r²` by mistake would not compile. After one period
//! from Kepler's third law, `T = 2π √(a³/GM)`, the comet is back at perihelion.

mod common;

use common::{gm_sun, write_csv};
use tensorcheck::spaces::plane::{Cartesian2, X, Y};
use tensorcheck::units::si::{Acceleration, Length, Time, Velocity, m, s};
use tensorcheck::{Point, Quantity, Vector};
use tensorcheck_numerics::ode::{DormandPrince54, Tolerance, solve_ivp_with};

type State = (Point<Cartesian2>, Vector<Velocity, Cartesian2>);

const AU: f64 = 1.495_978_707e11; // m
const YEAR: f64 = 365.25 * 86_400.0; // s (Julian year)

fn main() {
    let gm = gm_sun();
    let a: Quantity<Length> = 17.834 * AU * m; // semi-major axis
    let e = 0.967_14; // eccentricity
    let q = (1.0 - e) * a; // perihelion distance

    let rhs = |_t: Quantity<Time>, y: &State| {
        let x = y.0.coord(X) * m;
        let yy = y.0.coord(Y) * m;
        let r = y.0.coord(X).hypot(y.0.coord(Y)) * m;
        let k = gm / r.powi::<3>(); // 1/s²
        let acc: Vector<Acceleration, Cartesian2> = Vector::from_parts((-(k * x), -(k * yy)));
        (y.1, acc)
    };

    // Start at perihelion, moving with the vis-viva speed v² = GM (2/q − 1/a).
    let v_q: Quantity<Velocity> = (gm * (2.0 / q - 1.0 / a)).sqrt();
    let y0: State = (
        Point::new((q, 0.0 * m)),
        Vector::from_parts((0.0 * (m / s), v_q)),
    );
    let period: Quantity<Time> = 2.0 * std::f64::consts::PI * (a.powi::<3>() / gm).sqrt();
    println!(
        "period from Kepler's third law: {:.3} years",
        period.value_in(s) / YEAR
    );

    let energy = |y: &State| {
        let (vx, vy) = (y.1.base_components()[0], y.1.base_components()[1]);
        let r = y.0.coord(X).hypot(y.0.coord(Y));
        0.5 * (vx * vx + vy * vy) - gm.base_value() / r
    };
    let e0 = energy(&y0);
    let mut worst_e: f64 = 0.0;
    let mut rows = Vec::new();
    let tol = Tolerance {
        rtol: 1e-12,
        atol: (
            Point::new((1.0 * m, 1.0 * m)),
            Vector::from_parts((1e-9 * (m / s), 1e-9 * (m / s))),
        ),
    };
    let sol = solve_ivp_with(rhs, (0.0 * s)..period, y0, DormandPrince54, tol, |t, y| {
        worst_e = worst_e.max(((energy(y) - e0) / e0).abs());
        rows.push(vec![
            t.value_in(s) / YEAR,
            y.0.coord(X) / AU,
            y.0.coord(Y) / AU,
        ]);
    })
    .unwrap();

    let miss = (sol.y.0.coord(X) - q.value_in(m)).hypot(sol.y.0.coord(Y));
    println!(
        "after one period: {:.3e} m from the starting point ({:.1e} of the perihelion distance), {} steps",
        miss,
        miss / q.value_in(m),
        sol.stats.steps
    );
    println!("energy kept to {worst_e:e}");
    write_csv("t_years,x_au,y_au", &rows);
    assert!(miss / q.value_in(m) < 1e-6);
    assert!(worst_e < 1e-9);
}
