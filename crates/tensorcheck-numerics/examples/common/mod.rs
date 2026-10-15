//! Definitions shared by the gravity examples: constants, Schwarzschild spacetime
//! (defined here, outside the library) and CSV output for the figures.

#![allow(dead_code)]

use std::io::Write;

use tensorcheck::dim::Quot;
use tensorcheck::units::si::{Dimensionless, Length, Time, Velocity, m, s};
use tensorcheck::{Metric, Point, Quantity, Scalar};

/// `GM` has dimension m³/s².
pub type GravParam = Quot<
    tensorcheck::dim::Prod<Length, tensorcheck::dim::Prod<Length, Length>>,
    tensorcheck::dim::Prod<Time, Time>,
>;

/// Speed of light.
pub fn c() -> Quantity<Velocity> {
    299_792_458.0 * (m / s)
}
/// Heliocentric gravitational constant `GM☉` (IAU 2015 nominal value).
pub fn gm_sun() -> Quantity<GravParam> {
    1.327_124_4e20 * (m * m * m / (s * s))
}
/// Schwarzschild radius `r_s = 2GM/c²`.
pub fn schwarzschild_radius(gm: Quantity<GravParam>) -> Quantity<Length> {
    2.0 * gm / (c() * c())
}

tensorcheck::space!(
    /// The spacetime outside a spherical mass (reference dimension: length).
    pub Exterior, 4, Length
);
tensorcheck::chart!(
    /// Schwarzschild coordinates `(ct, r, θ, φ)`.
    pub Schw on Exterior { CT: Length, R: Length, Theta: Dimensionless, Phi: Dimensionless }
);

/// `ds² = (1 − r_s/r) c²dt² − dr²/(1 − r_s/r) − r²(dθ² + sin²θ dφ²)`.
#[derive(Clone, Copy, Debug)]
pub struct SchwarzschildMetric {
    rs: f64,
}
impl SchwarzschildMetric {
    /// The metric tensor for Schwarzschild radius `r_s`.
    pub fn new(rs: Quantity<Length>) -> Self {
        SchwarzschildMetric {
            rs: rs.base_value(),
        }
    }
}
impl Metric for SchwarzschildMetric {
    type Chart = Schw;
    type Dim = Dimensionless;
    fn components<N: Scalar>(&self, p: &Point<Schw, N>) -> [N; 10] {
        let r = p.coord(R);
        let f = N::one() - N::from_f64(self.rs) / r;
        let st = p.coord(Theta).sin();
        let z = N::zero();
        [
            f,
            z,
            z,
            z,
            -(N::one() / f),
            z,
            z,
            -(r * r),
            z,
            -(r * r * st * st),
        ]
    }
}

/// The path after `--csv` on the command line, if any.
pub fn csv_path() -> Option<String> {
    let args: Vec<String> = std::env::args().collect();
    args.iter()
        .position(|a| a == "--csv")
        .and_then(|i| args.get(i + 1).cloned())
}

/// Writes `rows` under `header` to the `--csv` path, if one was given.
pub fn write_csv(header: &str, rows: &[Vec<f64>]) {
    let Some(path) = csv_path() else { return };
    let mut f = std::fs::File::create(&path).expect("cannot create the CSV file");
    writeln!(f, "{header}").unwrap();
    for row in rows {
        let line: Vec<String> = row.iter().map(|x| format!("{x:e}")).collect();
        writeln!(f, "{}", line.join(",")).unwrap();
    }
}

/// The state of a geodesic: the point and the velocity `dx^a/dτ`.
pub type GeoState = (Point<Schw>, tensorcheck::Vector<Velocity, Schw>);

/// The geodesic equation `d²x^a/dτ² = −Γ^a_bc u^b u^c` as the right-hand side of a
/// first-order system in `(x, u)`.
#[macro_export]
macro_rules! geodesic_rhs {
    ($g:expr) => {
        |_tau: tensorcheck::Quantity<tensorcheck::units::si::Time>, y: &common::GeoState| {
            use tensorcheck::labels::{a, b, c};
            let gamma = tensorcheck::geometry::christoffel($g, &y.0).unwrap();
            (
                y.1,
                -(gamma.idx((a, b, c)) * y.1.idx(b) * y.1.idx(c)).into_tensor(),
            )
        }
    };
}
