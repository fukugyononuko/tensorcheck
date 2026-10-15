// VIOLATION: a fully qualified path in a trait bound.
pub fn f<T: crate::engine::index::Foo>() {}
