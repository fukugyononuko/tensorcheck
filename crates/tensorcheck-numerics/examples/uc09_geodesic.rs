//! UC-09 (second part): geodesics on a sphere. The right-hand side computes the
//! Christoffel symbols at each point; the path stays on a great circle and keeps
//! its speed.

use tensorcheck::geometry::christoffel;
use tensorcheck::labels::{a, b, c};
use tensorcheck::spaces::sphere2::{Angles, SphereMetric};
use tensorcheck::units::si::{Time, Velocity, m, rad, s};
use tensorcheck::{Point, Quantity, Vector, metric_tensor};
use tensorcheck_numerics::ode::{DormandPrince54, Tolerance, solve_ivp_with};

type U = Vector<Velocity, Angles>;

fn main() {
    let g = SphereMetric::new(2.0 * m);
    let rhs = |_l: Quantity<Time>, y: &(Point<Angles>, U)| {
        let gamma = christoffel(&g, &y.0).unwrap();
        (
            y.1,
            -(gamma.idx((a, b, c)) * y.1.idx(b) * y.1.idx(c)).into_tensor(),
        )
    };
    let p0 = Point::<Angles>::new((0.7 * rad, 0.2 * rad));
    let u0: U = Vector::from_parts((0.3 * (rad / s), 0.5 * (rad / s)));
    let speed2 = |p: &Point<Angles>, u: &U| {
        (metric_tensor(&g, p).idx((a, b)) * u.idx(a) * u.idx(b))
            .into_tensor()
            .value_in(m * m / (s * s))
    };
    let v0 = speed2(&p0, &u0);
    let mut worst: f64 = 0.0;
    let tol = Tolerance {
        rtol: 1e-11,
        atol: (
            Point::<Angles>::new((1e-12 * rad, 1e-12 * rad)),
            Vector::from_parts((1e-12 * (rad / s), 1e-12 * (rad / s))),
        ),
    };
    let sol = solve_ivp_with(
        rhs,
        (0.0 * s)..(4.0 * s),
        (p0, u0),
        DormandPrince54,
        tol,
        |_, y| {
            worst = worst.max((speed2(&y.0, &y.1) - v0).abs() / v0);
        },
    )
    .unwrap();
    println!(
        "end point (θ, φ) = ({:.6}, {:.6}) rad after {} steps; speed kept to {worst:e}",
        sol.y.0.get(0),
        sol.y.0.get(1),
        sol.stats.steps
    );
    assert!(worst < 1e-10);
}
