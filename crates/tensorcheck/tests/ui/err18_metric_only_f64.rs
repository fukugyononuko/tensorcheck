// TC-F18 / ERR-18: metric components written for f64 only (the compiler's own error).
use tensorcheck::spaces::plane::Polar;
use tensorcheck::units::si::Dimensionless;
use tensorcheck::{Metric, Point};

struct OnlyF64;
impl Metric for OnlyF64 {
    type Chart = Polar;
    type Dim = Dimensionless;
    fn components(&self, p: &Point<Polar, f64>) -> [f64; 3] {
        let r = p.get(0);
        [1.0, 0.0, r * r]
    }
}

fn main() {}
