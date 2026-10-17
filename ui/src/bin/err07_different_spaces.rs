// TC-F07 / ERR-07: contracting a spacetime index with a field-space index.
// expect: error[E0277]: index `mu` belongs to different spaces
// expect: error[E0277]: index `mu` is in chart `Inertial` in one place and `Phi` in the other
use tensorcheck::labels::mu;
use tensorcheck::spaces::minkowski::Inertial;
use tensorcheck::units::si::{Dimensionless, Length, m, one};
use tensorcheck::{Covector, Vector};

tensorcheck::space!(Fields, 2, Dimensionless);
tensorcheck::chart!(Phi on Fields { P1: Dimensionless, P2: Dimensionless });

fn main() {
    let x: Vector<Length, Inertial> = Vector::from_components([1.0, 0.0, 0.0, 0.0], m);
    let f: Covector<Dimensionless, Phi> = Covector::from_components([1.0, 2.0], one);
    let _ = (x.idx(mu) * f.idx(mu)).into_tensor();
}
