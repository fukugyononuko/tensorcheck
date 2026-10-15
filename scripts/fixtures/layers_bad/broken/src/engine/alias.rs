// VIOLATION: an alias of a manager module, used later.
use crate::manager as m;

pub fn f() {
    m::tensor::g();
}
