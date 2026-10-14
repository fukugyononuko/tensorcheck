//! TC-F: mistakes must not compile, and the headline must be ours (DD-16).

#[test]
fn compile_fail() {
    let t = trybuild::TestCases::new();
    t.compile_fail("tests/ui/*.rs");
}
