// A manager may use engines, also through a module import.
use crate::engine;
use crate::engine::index::Slot;

pub fn f() {
    engine::metric::first("");
}
