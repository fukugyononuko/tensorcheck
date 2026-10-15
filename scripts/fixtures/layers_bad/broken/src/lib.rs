#![forbid(unsafe_code)]
pub mod contracts;
pub mod engine;
pub mod manager;
pub use manager::tensor::Tensor;

/// A public facade over a manager, as in the real crate.
pub mod geometry {
    pub use crate::manager::tensor::christoffel;
}
