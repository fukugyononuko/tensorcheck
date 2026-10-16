// TC-F23 / ERR-22: the square root of a length (the exponent of m is odd).
use tensorcheck::Quantity;
use tensorcheck::units::si::{Length, m};

fn main() {
    let a: Quantity<Length> = Quantity::new(4.0, m);
    let _ = a.sqrt();
}
