//! UC-03: a two-field nonlinear sigma model on Minkowski spacetime,
//! L_kin = G_IJ(φ) ∂_μφ^I ∂^μφ^J, mixing spacetime and field-space indices.
//!
//! Contracting a spacetime index with a field-space index would not compile.

use tensorcheck::labels::{I, J, mu};
use tensorcheck::spaces::minkowski::{Inertial, MinkowskiMetric};
use tensorcheck::units::si::{Dimensionless, Length, m, one};
use tensorcheck::{Down, Metric, Point, Quot, Rank2, Scalar, Up, metric_tensor, slot};

tensorcheck::space!(Fields, 2, Dimensionless);
tensorcheck::chart!(Phi on Fields { P1: Dimensionless, P2: Dimensionless });

/// Field-space metric G = diag(1, e^{2φ¹}).
struct FieldMetric;
impl Metric for FieldMetric {
    type Chart = Phi;
    type Dim = Dimensionless;
    fn components<N: Scalar>(&self, p: &Point<Phi, N>) -> [N; 3] {
        [N::one(), N::zero(), (p.get(0) * N::from_f64(2.0)).exp()]
    }
}

fn main() {
    let phi = Point::<Phi>::new((0.3 * one, -1.0 * one));
    let x = Point::<Inertial>::new((0.0 * m, 0.0 * m, 0.0 * m, 0.0 * m));
    let d = [[0.5, -0.2], [0.1, 0.3], [0.0, 0.4], [-0.6, 0.2]];
    let dphi: Rank2<Quot<Dimensionless, Length>, Down<Inertial>, Up<Phi>> =
        Rank2::from_components(d, one / m);
    let dphi_up = dphi.raise(slot::_0, &MinkowskiMetric, &x).unwrap();
    let g = metric_tensor(&FieldMetric, &phi);
    let kin = (g.idx((I, J)) * dphi.idx((mu, I)) * dphi_up.idx((mu, J))).into_tensor();
    let value = kin.value_in(one / (m * m));
    println!("G_IJ ∂_μφ^I ∂^μφ^J = {value} m⁻²");

    let gdiag = [1.0, (0.6f64).exp()];
    let eta = [1.0, -1.0, -1.0, -1.0];
    let want: f64 = (0..4)
        .map(|m_| {
            (0..2)
                .map(|i| gdiag[i] * eta[m_] * d[m_][i] * d[m_][i])
                .sum::<f64>()
        })
        .sum();
    assert!((value - want).abs() <= 1e-14 * want.abs());
}
