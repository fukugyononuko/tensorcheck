//! UC-09 (first part): connection and curvature of a sphere of radius a, from the
//! metric tensor alone (derivatives by automatic differentiation).

use tensorcheck::geometry::curvature;
use tensorcheck::spaces::sphere2::{
    Angles, AnglesToStereographic, SphereMetric, StereographicMetric,
};
use tensorcheck::units::si::{m, rad};
use tensorcheck::{ChartTransform, Point};

fn main() {
    let a = 2.0;
    let g = SphereMetric::new(a * m);
    let p = Point::<Angles>::new((0.9 * rad, 0.3 * rad));
    let cv = curvature(&g, &p).unwrap();
    let r = cv.ricci_scalar.value_in(m.pow::<-2>());
    println!(
        "scalar curvature R = {r} m⁻² (expected 2/a² = {})",
        2.0 / (a * a)
    );
    let gamma = cv.christoffel.base_components();
    println!(
        "Γ^θ_φφ = {} (expected −sin θ cos θ = {})",
        gamma[0][1][1],
        -(0.9f64.sin() * 0.9f64.cos())
    );

    // The same scalar in stereographic coordinates.
    let xf = AnglesToStereographic::new(a * m);
    let r2 = curvature(&StereographicMetric::new(a * m), &xf.map(&p))
        .unwrap()
        .ricci_scalar
        .value_in(m.pow::<-2>());
    println!("R in stereographic coordinates = {r2} m⁻²");
    assert!((r - 0.5).abs() < 1e-14 && (r2 - 0.5).abs() < 1e-14);
}
