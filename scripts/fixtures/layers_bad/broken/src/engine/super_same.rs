// VIOLATION: another engine reached through super in an expression.
pub fn f() {
    super::index::g();
}
