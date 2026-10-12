//! Symmetry of tensors (CMP-E04, DD-07).

/// No symmetry.
#[derive(Clone, Copy, Debug)]
pub struct NoSym;
/// Symmetric in its two indices: `T_ab = T_ba`.
#[derive(Clone, Copy, Debug)]
pub struct Symmetric;
/// Antisymmetric in its two indices: `T_ab = -T_ba`.
#[derive(Clone, Copy, Debug)]
pub struct Antisymmetric;
