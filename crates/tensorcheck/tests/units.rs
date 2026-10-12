//! TC-U02 (units), TC-U03 (size of rank 0), TC-I01 scalar part, TC-R01, TC-R02.

use tensorcheck::units::natural::{self, GeV, MeV, reduced_planck_mass};
use tensorcheck::units::si::{
    Acceleration, Energy, Force, J, Length, Mass, Time, Velocity, au, g, h, kg, km, m, min, s,
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
