//! Natural units with `ħ = c = 1`: the only base dimension is mass (API-06).
//!
//! The base unit is the GeV. Energy, mass, momentum and inverse length all have
//! mass dimension.

use crate::engine::dimension::{Dim, UnitSystem};
use crate::unit;
use crate::utility::typelevel::{Cons, N1, N2, Nil, P1, P2, P3, P4, Z0};

/// The natural unit system (`ħ = c = 1`).
pub struct Natural;
impl UnitSystem for Natural {
    type Dimensionless = Cons<Z0, Nil>;
    const SYMBOLS: &'static [&'static str] = &["GeV"];
}

/// A natural-unit dimension: mass to the power `E`.
pub type MassDim<E> = Dim<Natural, Cons<E, Nil>>;
/// Dimensionless.
pub type Dimensionless = MassDim<Z0>;
/// Mass (and energy, momentum).
pub type Mass = MassDim<P1>;
/// Mass squared.
pub type Mass2 = MassDim<P2>;
/// Mass cubed.
pub type Mass3 = MassDim<P3>;
/// Mass to the fourth (energy density, Lagrangian density).
pub type Mass4 = MassDim<P4>;
/// Length and time (inverse mass).
pub type Length = MassDim<N1>;
/// Inverse mass squared (for example `G_N`).
pub type InverseMass2 = MassDim<N2>;

unit!(/** Gigaelectronvolt (the base unit). */ pub GeV: Mass = 1.0);
unit!(/** Electronvolt. */ pub eV: Mass = 1.0e-9);
unit!(/** Kiloelectronvolt. */ pub keV: Mass = 1.0e-6);
unit!(/** Megaelectronvolt. */ pub MeV: Mass = 1.0e-3);
unit!(/** Teraelectronvolt. */ pub TeV: Mass = 1.0e3);
unit!(
    /// Reduced Planck mass `M_P = (ħc / 8πG)^(1/2)`.
    ///
    /// Computed from `ħ`, `c` and `G_N` of PDG 2024 (CODATA 2018); the test
    /// `reduced_planck_mass_matches_constants` recomputes it.
    pub reduced_planck_mass: Mass = 2.435_323_459_338_203_6e18
);
