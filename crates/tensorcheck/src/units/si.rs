//! The International System of Units: seven base dimensions (API-05).
//!
//! Base dimensions, in order: length, mass, time, electric current, temperature,
//! amount of substance, luminous intensity.

use crate::engine::dimension::{Dim, UnitSystem};
use crate::units::unit;
use crate::utility::typelevel::{Cons, N1, N2, N3, Nil, P1, P2, P3, Z0};

/// The SI unit system.
pub struct Si;
impl UnitSystem for Si {
    type Dimensionless = Exps<Z0, Z0, Z0, Z0, Z0, Z0, Z0>;
    const SYMBOLS: &'static [&'static str] = &["m", "kg", "s", "A", "K", "mol", "cd"];
}

/// Exponents of the seven SI base dimensions.
pub type Exps<L, M, T, I, Th, Amt, J> =
    Cons<L, Cons<M, Cons<T, Cons<I, Cons<Th, Cons<Amt, Cons<J, Nil>>>>>>>;
/// An SI dimension from its seven exponents.
pub type SiDim<L, M, T, I = Z0, Th = Z0, Amt = Z0, J = Z0> = Dim<Si, Exps<L, M, T, I, Th, Amt, J>>;

/// Dimensionless.
pub type Dimensionless = SiDim<Z0, Z0, Z0>;
/// Length.
pub type Length = SiDim<P1, Z0, Z0>;
/// Mass.
pub type Mass = SiDim<Z0, P1, Z0>;
/// Time.
pub type Time = SiDim<Z0, Z0, P1>;
/// Electric current.
pub type Current = SiDim<Z0, Z0, Z0, P1>;
/// Thermodynamic temperature.
pub type Temperature = SiDim<Z0, Z0, Z0, Z0, P1>;
/// Amount of substance.
pub type Amount = SiDim<Z0, Z0, Z0, Z0, Z0, P1>;
/// Luminous intensity.
pub type Luminosity = SiDim<Z0, Z0, Z0, Z0, Z0, Z0, P1>;
/// Area.
pub type Area = SiDim<P2, Z0, Z0>;
/// Volume.
pub type Volume = SiDim<P3, Z0, Z0>;
/// Frequency.
pub type Frequency = SiDim<Z0, Z0, N1>;
/// Velocity.
pub type Velocity = SiDim<P1, Z0, N1>;
/// Acceleration.
pub type Acceleration = SiDim<P1, Z0, N2>;
/// Momentum.
pub type Momentum = SiDim<P1, P1, N1>;
/// Force.
pub type Force = SiDim<P1, P1, N2>;
/// Energy.
pub type Energy = SiDim<P2, P1, N2>;
/// Power.
pub type Power = SiDim<P2, P1, N3>;
/// Pressure.
pub type Pressure = SiDim<N1, P1, N2>;
/// Electric charge.
pub type Charge = SiDim<Z0, Z0, P1, P1>;
/// Electric potential.
pub type Voltage = SiDim<P2, P1, N3, N1>;

unit!(/** Metre. */ m: Length = 1.0);
unit!(/** Kilometre. */ km: Length = 1.0e3);
unit!(/** Centimetre. */ cm: Length = 1.0e-2);
unit!(/** Millimetre. */ mm: Length = 1.0e-3);
unit!(/** Astronomical unit (IAU 2012, exact). */ au: Length = 1.495_978_707e11);
unit!(/** Kilogram. */ kg: Mass = 1.0);
unit!(/** Gram. */ g: Mass = 1.0e-3);
unit!(/** Second. */ s: Time = 1.0);
unit!(/** Minute. */ min: Time = 60.0);
unit!(/** Hour. */ h: Time = 3600.0);
unit!(/** Ampere. */ A: Current = 1.0);
unit!(/** Kelvin. */ K: Temperature = 1.0);
unit!(/** Mole. */ mol: Amount = 1.0);
unit!(/** Candela. */ cd: Luminosity = 1.0);
unit!(/** Radian (dimensionless). */ rad: Dimensionless = 1.0);
unit!(/** The number one (dimensionless). */ one: Dimensionless = 1.0);
unit!(/** Hertz. */ Hz: Frequency = 1.0);
unit!(/** Newton. */ N: Force = 1.0);
unit!(/** Joule. */ J: Energy = 1.0);
unit!(/** Watt. */ W: Power = 1.0);
unit!(/** Pascal. */ Pa: Pressure = 1.0);
unit!(/** Coulomb. */ C: Charge = 1.0);
unit!(/** Volt. */ V: Voltage = 1.0);
unit!(/** Electronvolt (exact since the 2019 SI). */ eV: Energy = 1.602_176_634e-19);
