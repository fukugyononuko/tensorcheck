//! UC-02: speed squared in polar coordinates, g_ab V^a V^b.
//!
//! The components of a velocity in the coordinate basis have different dimensions:
//! V^r in m/s and V^θ in 1/s.

use tensorcheck::labels::{a, b};
use tensorcheck::spaces::plane::{EuclidMetric2, Polar, PolarMetric, PolarToCartesian2, Theta};
use tensorcheck::units::si::{Frequency, Velocity, m, rad, s};
use tensorcheck::{ChartTransform, Point, Quantity, Vector, slot};

fn main() {
    let p = Point::<Polar>::new((2.0 * m, 0.5 * rad));
    let v: Vector<Velocity, Polar> = Vector::from_parts((1.5 * (m / s), 0.25 * (rad / s)));
    let v_theta: Quantity<Frequency> = v.at((Theta,));

    let v2 = (v.lower(slot::_0, &PolarMetric, &p).idx(a) * v.idx(a)).into_tensor();
    let v2 = v2.value_in(m * m / (s * s));
    println!(
        "V^θ = {} rad/s, |v|² = {v2} m²/s² (polar)",
        v_theta.value_in(rad / s)
    );

    // The same speed from Cartesian components.
    let vc = v.transform(&PolarToCartesian2, &p).unwrap();
    let q = PolarToCartesian2.map(&p);
    let v2c = (vc.lower(slot::_0, &EuclidMetric2, &q).idx(b) * vc.idx(b))
        .into_tensor()
        .value_in(m * m / (s * s));
    println!("|v|² = {v2c} m²/s² (Cartesian)");
    assert!((v2 - 2.5).abs() < 1e-14 && (v2c - 2.5).abs() < 1e-14);
}
