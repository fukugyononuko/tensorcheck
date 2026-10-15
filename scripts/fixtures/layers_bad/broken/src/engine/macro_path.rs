// VIOLATION: a manager path inside a macro body.
macro_rules! m {
    () => {
        $crate::manager::tensor::g()
    };
}
