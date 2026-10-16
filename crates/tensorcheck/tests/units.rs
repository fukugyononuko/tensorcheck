//! TC-U02 (units), TC-U03 (size of rank 0), TC-I01 scalar part, TC-R01, TC-R02.

use tensorcheck::units::natural::{self, GeV, MeV, reduced_planck_mass};
use tensorcheck::units::si::{
    Acceleration, Area, Energy, Force, Frequency, J, Length, Mass, Time, Velocity, Volume, au, g,
    h, kg, km, m, min, one, s,
};
use tensorcheck::{Dim, Quantity, Unit, UnitSystem};

/// Relative difference.
fn rel(a: f64, b: f64) -> f64 {
    ((a - b) / b).abs()
}

// Round trips multiply and divide by the same factor: at most two roundings.
const ROUND_TRIP_TOL: f64 = 2.0 * f64::EPSILON;

#[test]
fn conversions_between_units() {
    let d: Quantity<Length> = Quantity::new(1.0, km);
    assert_eq!(d.value_in(m), 1000.0);
    let t: Quantity<Time> = Quantity::new(1.5, h);
    assert_eq!(t.value_in(min), 90.0);
    let mass: Quantity<Mass> = Quantity::new(250.0, g);
    assert_eq!(mass.value_in(kg), 0.25);
    let e: Quantity<natural::Mass> = Quantity::new(1.0, GeV);
    assert!(rel(e.value_in(MeV), 1000.0) <= ROUND_TRIP_TOL);
}

#[test]
fn round_trip_is_exact_to_rounding() {
    // Inputs stay well inside the f64 range after scaling by the factor (about 1.5e11).
    for &x in &[1e-290, 1.234_567_89, 6.02e23, 1e290] {
        let d: Quantity<Length> = Quantity::new(x, au);
        assert!(rel(d.value_in(au), x) <= ROUND_TRIP_TOL, "{x}");
    }
}

#[test]
fn composite_units() {
    let v: Quantity<Velocity> = Quantity::new(36.0, km / h);
    assert!(rel(v.value_in(m / s), 10.0) <= ROUND_TRIP_TOL);
    let a: Quantity<Acceleration> = 9.8 * (m / s.pow::<2>());
    assert_eq!(a.value_in(m * s.pow::<-2>()), 9.8);
}

#[test]
fn newtons_second_law_for_scalars() {
    let mass = 2.0 * kg;
    let acc = 9.8 * (m / s.pow::<2>());
    let force: Quantity<Force> = mass * acc;
    assert!(rel(force.value_in(tensorcheck::units::si::N), 19.6) <= ROUND_TRIP_TOL);
    let work: Quantity<Energy> = force * (3.0 * m);
    assert!(rel(work.value_in(J), 58.8) <= 4.0 * f64::EPSILON);
    let speed: Quantity<Velocity> = (10.0 * m) / (4.0 * s);
    assert_eq!(speed.value_in(m / s), 2.5);
}

#[test]
fn arithmetic_on_quantities() {
    let a = 3.0 * m;
    let b = 2.0 * m;
    assert_eq!((a + b).value_in(m), 5.0);
    assert_eq!((a - b).value_in(m), 1.0);
    assert_eq!((-a).value_in(m), -3.0);
    assert_eq!((2.0 * a).value_in(m), 6.0);
    assert_eq!((a / 2.0).value_in(m), 1.5);
    assert!(b < a);
}

#[test]
fn debug_prints_base_units() {
    let f = (2.0 * kg) * (3.0 * (m / s.pow::<2>()));
    assert_eq!(format!("{f:?}"), "6.0 m kg s^-2");
}

#[test]
fn reduced_planck_mass_matches_constants() {
    // PDG 2024 (CODATA 2018): ħ and c exact, G_N = 6.67430e-11 m^3 kg^-1 s^-2.
    let hbar = 1.054_571_817e-34_f64;
    let c = 299_792_458.0_f64;
    let g_n = 6.674_30e-11_f64;
    let ev = 1.602_176_634e-19_f64;
    let m_p_gev = (hbar * c / (8.0 * std::f64::consts::PI * g_n)).sqrt() * c * c / ev / 1e9;
    let m_p: Quantity<natural::Mass> = Quantity::new(1.0, reduced_planck_mass);
    assert!(rel(m_p.value_in(GeV), m_p_gev) <= 4.0 * f64::EPSILON);
}

#[test]
fn rank_zero_has_no_overhead() {
    // TC-U03 / NFR-03: a quantity is exactly one f64.
    assert_eq!(size_of::<Quantity<Length>>(), size_of::<f64>());
}

// TC-R01 (VOL-01): a user-defined geometric unit system (c = G = 1, length only).
#[allow(non_camel_case_types)]
mod geometric {
    use tensorcheck::utility::typelevel::{Cons, Nil, P1, Z0};
    use tensorcheck::{Dim, Unit, UnitSystem};

    pub struct Geometric;
    impl UnitSystem for Geometric {
        type Dimensionless = Cons<Z0, Nil>;
        const SYMBOLS: &'static [&'static str] = &["m"];
    }
    pub type GLength = Dim<Geometric, Cons<P1, Nil>>;

    #[derive(Clone, Copy)]
    pub struct solar_mass_length;
    impl Unit for solar_mass_length {
        type Dim = GLength;
        // G M_sun / c^2 in metres (IAU 2015 nominal G M_sun = 1.3271244e20 m^3 s^-2).
        const FACTOR: f64 = 1.327_124_4e20 / (299_792_458.0 * 299_792_458.0);
    }
    #[derive(Clone, Copy)]
    pub struct metre;
    impl Unit for metre {
        type Dim = GLength;
        const FACTOR: f64 = 1.0;
    }
}

#[test]
fn user_defined_unit_system() {
    use geometric::*;
    let r: Quantity<GLength> = Quantity::new(1.0, solar_mass_length);
    assert!(rel(r.value_in(metre), 1476.625) < 1e-6);
    let r2 = r + r;
    assert!(rel(r2.value_in(metre), 2953.25) < 1e-6);
}

// TC-R02 (VOL-02): a user-defined unit in an existing system.
#[derive(Clone, Copy)]
#[allow(non_camel_case_types)]
struct light_year;
impl Unit for light_year {
    type Dim = Length;
    const FACTOR: f64 = 9.460_730_472_580_8e15;
}

#[test]
fn user_defined_unit() {
    let d: Quantity<Length> = Quantity::new(1.0, light_year);
    assert!(rel(d.value_in(au), 63_241.077) < 1e-7);
}

#[allow(dead_code)]
fn unit_system_is_usable_in_bounds<S: UnitSystem>() {}
#[allow(dead_code)]
type CheckDim =
    Dim<tensorcheck::units::si::Si, <tensorcheck::units::si::Si as UnitSystem>::Dimensionless>;

// TC-U02 (issue #19): square roots and integer powers of quantities. The types in the
// annotations are the check on the dimensions; the values are one rounding each.
#[test]
fn sqrt_halves_every_exponent() {
    let side: Quantity<Length> = Quantity::of(9.0, m.pow::<2>()).sqrt();
    assert_eq!(side.value_in(m), 3.0);
    let v: Quantity<Velocity> = Quantity::of(16.0, m.pow::<2>() / s.pow::<2>()).sqrt();
    assert_eq!(v.value_in(m / s), 4.0);
    // Negative exponents: GM / a^3 has dimension s^-2, so its root is an angular rate.
    let omega: Quantity<Frequency> =
        (Quantity::of(4.0, m.pow::<3>() / s.pow::<2>()) / Quantity::of(1.0, m.pow::<3>())).sqrt();
    assert_eq!(omega.value_in(one / s), 2.0);
    // Natural units: GeV^2 -> GeV.
    let mass: Quantity<natural::Mass> = Quantity::of(2.25, GeV.pow::<2>()).sqrt();
    assert!(rel(mass.value_in(GeV), 1.5) <= ROUND_TRIP_TOL);
    // Dimensionless quantities, and NaN for negative values as for f64::sqrt.
    assert_eq!(Quantity::of(0.25, one).sqrt().value_in(one), 0.5);
    assert!(Quantity::of(-1.0, m.pow::<2>()).sqrt().value_in(m).is_nan());
}

#[test]
fn powi_multiplies_every_exponent() {
    let a = Quantity::of(2.0, m);
    let vol: Quantity<Volume> = a.powi::<3>();
    assert_eq!(vol.value_in(m.pow::<3>()), 8.0);
    let per_area = a.powi::<-2>();
    assert_eq!(per_area.value_in(m.pow::<-2>()), 0.25);
    let unity = a.powi::<0>();
    assert_eq!(unity.value_in(one), 1.0);
    let same: Quantity<Length> = a.powi::<1>();
    assert_eq!(same.value_in(m), 2.0);
    // sqrt undoes powi::<2> up to the sign.
    let back: Quantity<Length> = Quantity::of(-3.0, m).powi::<2>().sqrt();
    assert_eq!(back.value_in(m), 3.0);
}

#[test]
fn sqrt_and_powi_carry_derivatives() {
    use tensorcheck::spaces::plane::Plane;
    use tensorcheck::utility::autodiff::Dual;
    // d/dx sqrt(x) = 1 / (2 sqrt(x)), d/dx x^3 = 3 x^2 at x = 4.
    let x: Quantity<Area, Dual<Plane>> = Quantity::from_base_components(Dual::variable(4.0, 0));
    let r: Quantity<Length, Dual<Plane>> = x.sqrt();
    assert_eq!((r.base_value().re, r.base_value().d[0]), (2.0, 0.25));
    let c = x.powi::<3>().base_value();
    assert_eq!((c.re, c.d[0]), (64.0, 48.0));
}
