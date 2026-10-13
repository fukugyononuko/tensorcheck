//! Compile-time checks for physical dimensions and tensor structure.
//!
//! Work in progress. See `docs/` for the design.

#![forbid(unsafe_code)]
#![deny(missing_docs)]
// The checks are done by type-level computation, so long types are the norm here.
#![allow(clippy::type_complexity)]

pub mod access;
pub mod contracts;
pub mod engine;
pub mod labels;
mod macros;
pub mod manager;
pub mod spaces;
pub mod units;
pub mod utility;

pub use contracts::{
    Chart, ChartBound, CoordinateBasis, Down, Label, NonCoordinateBasis, Point, Space, Tensorial,
    Up, slot,
};
pub use engine::dimension::{Dim, Pow, Prod, Quot, Unit, UnitSystem};
pub use engine::index::ChartTransform;
pub use engine::metric::Metric;
pub use engine::symmetry::{Antisymmetric, NoSym, Symmetric};
pub use manager::tensor::{
    Covector, Labeled, Quantity, Rank2, Rank3, Rank4, Tensor, Vector, inverse_metric_tensor,
    metric_tensor,
};
pub use utility::linear::LinearError;
pub use utility::scalar::Scalar;

/// Dimension arithmetic on types (API-07).
pub mod dim {
    pub use crate::engine::dimension::{Pow, Prod, Quot};
}
