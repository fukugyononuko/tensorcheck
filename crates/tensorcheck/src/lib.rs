//! Compile-time checks for physical dimensions and tensor structure.
//!
//! Work in progress. See `docs/` for the design.

#![forbid(unsafe_code)]
#![deny(missing_docs)]

pub mod access;
pub mod engine;
pub mod manager;
pub mod units;
pub mod utility;

pub use engine::dimension::{Dim, Pow, Prod, Quot, Unit, UnitSystem};
pub use engine::index::{ChartBound, Tensorial};
pub use engine::symmetry::{Antisymmetric, NoSym, Symmetric};
pub use manager::tensor::{Quantity, Tensor};
pub use utility::scalar::Scalar;

/// Dimension arithmetic on types (API-07).
pub mod dim {
    pub use crate::engine::dimension::{Pow, Prod, Quot};
}
