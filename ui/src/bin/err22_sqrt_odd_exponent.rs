// TC-F23 / ERR-22: the square root of a length (the exponent of m is odd).
// expect: error[E0277]: cannot take the square root: the dimension has an odd exponent
use tensorcheck::Quantity;
use tensorcheck::units::si::{Length, m};

fn main() {
    let a: Quantity<Length> = Quantity::new(4.0, m);
    let _ = a.sqrt();
}
