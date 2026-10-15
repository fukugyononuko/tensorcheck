// VIOLATION: import the layer module, then call another engine relative to it.
use crate::engine;

pub fn f() {
    engine::index::g();
}
