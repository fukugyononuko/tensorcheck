// VIOLATION: a manager function through a facade module of the crate root.
pub fn f() {
    crate::geometry::christoffel();
}
