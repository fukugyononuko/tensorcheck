//! Dimension-aware numerical solvers for `tensorcheck`: definite integrals, ordinary
//! differential equations (later phases) and numerical differentiation (later
//! phases). Inputs and outputs keep their dimensions and tensor structure.

#![forbid(unsafe_code)]
#![deny(missing_docs)]

pub mod diff;
mod engine;
mod manager;
pub mod ode;
pub mod quad;
