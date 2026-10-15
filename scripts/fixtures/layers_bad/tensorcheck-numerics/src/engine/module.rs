// VIOLATION: import internal, then use its engine relative to it.
use tensorcheck::internal;

pub fn f() {
    internal::engine::central();
}
