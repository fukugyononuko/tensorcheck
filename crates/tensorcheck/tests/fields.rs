//! Phase 4: symmetry, field space with spacetime, user definitions.
//! TC-U06, TC-I03, TC-I05, API-34.

use tensorcheck::labels::{I, J, a, b, mu};
use tensorcheck::spaces::euclid3::Cartesian;
use tensorcheck::spaces::minkowski::{Inertial, MinkowskiMetric};
use tensorcheck::spaces::plane::{Polar, R, Theta};
use tensorcheck::units::si::{Dimensionless, Frequency, Length, Velocity, m, one, rad, s};
use tensorcheck::{
    Antisymmetric, Down, Point, Quantity, Rank2, Symmetric, Up, Vector, metric_tensor, slot,
};

fn rel(x: f64, y: f64) -> f64 {
    (x - y).abs() / y.abs().max(f64::MIN_POSITIVE)
}

#[test]
fn symmetric_and_antisymmetric_construction() {
    // TC-U06.
    let s3: Rank2<Dimensionless, Down<Cartesian>, Down<Cartesian>, Symmetric> =
        Rank2::from_independent([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], one);
    let c = s3.components_in(one);
    assert_eq!(c, [[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]]);

    let f: Rank2<Dimensionless, Up<Cartesian>, Up<Cartesian>, Antisymmetric> =
        Rank2::from_independent([1.0, 2.0, 3.0], one);
    let c = f.components_in(one);
    assert_eq!(c, [[0.0, 1.0, 2.0], [-1.0, 0.0, 3.0], [-2.0, -3.0, 0.0]]);

    let t: Rank2<Dimensionless, Down<Cartesian>, Down<Cartesian>> =
        Rank2::from_components([[1.0, 2.0, 0.0], [4.0, 5.0, 6.0], [0.0, 0.0, 9.0]], one);
    let sym = t.symmetrize().components_in(one);
    let anti = t.antisymmetrize().components_in(one);
    for i in 0..3 {
        for j in 0..3 {
            assert_eq!(sym[i][j], sym[j][i]);
            assert_eq!(anti[i][j], -anti[j][i]);
            assert_eq!(
                sym[i][j] + anti[i][j],
                [[1.0, 2.0, 0.0], [4.0, 5.0, 6.0], [0.0, 0.0, 9.0]][i][j]
            );
        }
    }
    // Symmetric + symmetric stays symmetric; scaling keeps the marker.
    let doubled: Rank2<Dimensionless, Down<Cartesian>, Down<Cartesian>, Symmetric> = s3 + s3 * 1.0;
    assert_eq!(doubled.components_in(one)[0][1], 4.0);
}

#[test]
fn vector_from_parts_in_polar_coordinates() {
    // API-34: components with their own dimensions.
    let v: Vector<Velocity, Polar> = Vector::from_parts((1.0 * (m / s), 0.25 * (rad / s)));
    let vr: Quantity<Velocity> = v.at((R,));
    let vt: Quantity<Frequency> = v.at((Theta,));
    assert_eq!(vr.value_in(m / s), 1.0);
    assert_eq!(vt.value_in(rad / s), 0.25);
}

// A two-component field space with metric G = diag(1, exp(2 φ¹)) (a standard
// nonlinear sigma-model example).
tensorcheck::space!(Fields, 2, Dimensionless);
tensorcheck::chart!(Phi on Fields { P1: Dimensionless, P2: Dimensionless });
struct FieldMetric;
impl tensorcheck::Metric for FieldMetric {
    type Chart = Phi;
    type Dim = Dimensionless;
    fn components<N: tensorcheck::Scalar>(&self, p: &Point<Phi, N>) -> [N; 3] {
        let e = (p.get(0) * N::from_f64(2.0)).exp();
        [N::one(), N::zero(), e]
    }
}

#[test]
fn kinetic_term_mixes_spacetime_and_field_indices() {
    // TC-I03 / CC-03: G_IJ ∂_μφ^I ∂^μφ^J with η = diag(+,−,−,−).
    let phi = Point::<Phi>::new((0.3 * one, -1.0 * one));
    let x = Point::<Inertial>::new((0.0 * m, 0.0 * m, 0.0 * m, 0.0 * m));
    // ∂_μ φ^I: rows μ = 0..3, columns I = 1, 2 (dimension 1/m).
    let d: [[f64; 2]; 4] = [[0.5, -0.2], [0.1, 0.3], [0.0, 0.4], [-0.6, 0.2]];
    let dphi: Rank2<tensorcheck::Quot<Dimensionless, Length>, Down<Inertial>, Up<Phi>> =
        Rank2::from_components(d, one / m);
    let dphi_up = dphi.raise(slot::_0, &MinkowskiMetric, &x).unwrap();
    let g = metric_tensor(&FieldMetric, &phi);
    let kin = (g.idx((I, J)) * dphi.idx((mu, I)) * dphi_up.idx((mu, J))).into_tensor();

    let gdiag = [1.0, (2.0f64 * 0.3).exp()];
    let eta = [1.0, -1.0, -1.0, -1.0];
    let mut want = 0.0;
    for (mu_, e) in eta.iter().enumerate() {
        for (i, gi) in gdiag.iter().enumerate() {
            want += gi * e * d[mu_][i] * d[mu_][i];
        }
    }
    // A handful of products and sums: 10 ε relative.
    assert!(rel(kin.value_in(one / (m * m)), want) <= 10.0 * f64::EPSILON);
}

// TC-I05 / CC-06: everything user-defined (unit system, space, chart, metric).
mod custom {
    use tensorcheck::utility::typelevel::{Cons, Nil, P1, Z0};
    use tensorcheck::{Dim, Point, Scalar, UnitSystem};
    pub struct Lengths;
    impl UnitSystem for Lengths {
        type Dimensionless = Cons<Z0, Nil>;
        const SYMBOLS: &'static [&'static str] = &["u"];
    }
    pub type L = Dim<Lengths, Cons<P1, Nil>>;
    pub type D0 = Dim<Lengths, Cons<Z0, Nil>>;
    tensorcheck::unit!(/** The base unit of length. */ pub u: L = 1.0);
    tensorcheck::unit!(/** One. */ pub unity: D0 = 1.0);
    tensorcheck::space!(pub Disk, 2, L);
    tensorcheck::chart!(pub DiskPolar on Disk { Rr: L, Tt: D0 });
    pub struct DiskMetric;
    impl tensorcheck::Metric for DiskMetric {
        type Chart = DiskPolar;
        type Dim = D0;
        fn components<N: Scalar>(&self, p: &Point<DiskPolar, N>) -> [N; 3] {
            let r = p.get(0);
            [N::one(), N::zero(), r * r]
        }
    }
}

#[test]
fn everything_user_defined() {
    use custom::*;
    let p = Point::<DiskPolar>::new((2.0 * u, 0.5 * unity));
    let v: Vector<L, DiskPolar> = Vector::from_components([1.5, 0.25], u);
    let v2 = (v.lower(slot::_0, &DiskMetric, &p).idx(a) * v.idx(a)).into_tensor();
    assert!(rel(v2.value_in(u * u), 2.5) <= 4.0 * f64::EPSILON);
    let _ = b;
}
