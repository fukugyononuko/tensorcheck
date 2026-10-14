//! Items for `tensorcheck-numerics` only (01 §1.7). Not part of the public API;
//! user code must not use this module (checked by `scripts/check_layers.py`).

/// Engines: only the numerics crate's manager may use these.
pub mod engine {
    pub use crate::engine::differentiation::{central, jacobian_fd, ridders};
}

/// Utilities: any layer of the numerics crate may use these.
pub mod utility {
    pub use crate::utility::linear::{LinearError, lu_decompose, lu_decompose_with, lu_solve};
}
