// TC-F02 / ERR-02: a length created with a unit of time.
use tensorcheck::Quantity;
use tensorcheck::units::si::{Length, s};

fn main() {
    let _: Quantity<Length> = Quantity::new(1.0, s);
}
