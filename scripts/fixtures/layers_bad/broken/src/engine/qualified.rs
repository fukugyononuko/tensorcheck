// VIOLATION: a fully qualified path in an expression.
pub fn f() {
    crate::manager::tensor::g();
}
