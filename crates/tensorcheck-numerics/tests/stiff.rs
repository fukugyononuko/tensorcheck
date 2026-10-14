//! TC-U25: stiff problems with the Rosenbrock method (RODAS).

use tensorcheck::units::si::{Dimensionless, Time, one, s};
use tensorcheck::{Quantity, Vector};
use tensorcheck_numerics::ode::{
    DormandPrince54, Jacobian, Rosenbrock, SolveError, Tolerance, solve_ivp,
};

tensorcheck::space!(Species, 3, Dimensionless);
tensorcheck::chart!(Conc on Species { Y1: Dimensionless, Y2: Dimensionless, Y3: Dimensionless });
type Y = Vector<Dimensionless, Conc>;
type Rate = Vector<tensorcheck::units::si::Frequency, Conc>;

/// ROBER (Robertson 1966; Hairer & Wanner), with rates per second.
fn robertson(_t: Quantity<Time>, y: &Y) -> Rate {
    let c = y.components_in(one);
    let (y1, y2, y3) = (c[0], c[1], c[2]);
    Vector::from_components(
        [
            -0.04 * y1 + 1e4 * y2 * y3,
            0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2 * y2,
            3e7 * y2 * y2,
        ],
        one / s,
    )
}

/// Reference at t = 1e11 s: Table II.10.3 of the ROBER report in the INdAM-Bari
/// test set for IVP solvers (computed by RADAU).
#[allow(clippy::excessive_precision)] // digits as published in the test-set report
const REF_1E11: [f64; 3] = [
    0.2083340149701255e-7,
    0.8333360770334713e-13,
    0.9999999791665050,
];

fn solve_robertson(
    rtol: f64,
    t_end: f64,
) -> Result<tensorcheck_numerics::ode::Solution<Time, Y>, SolveError<Time>> {
    let y0: Y = Vector::from_components([1.0, 0.0, 0.0], one);
    let atol: Y = Vector::from_components([1e-6 * rtol, 1e-12 * rtol, 1e-6 * rtol], one);
    solve_ivp(
        robertson,
        (0.0 * s)..(t_end * s),
        y0,
        Rosenbrock {
            jacobian: Jacobian::Numerical,
        },
        Tolerance { rtol, atol },
    )
}

#[test]
#[ignore]
fn measure_robertson() {
    for rtol in [1e-4, 1e-6, 1e-8] {
        let sol = solve_robertson(rtol, 1e11).unwrap();
        let c = sol.y.components_in(one);
        let rel: Vec<f64> = (0..3)
            .map(|i| ((c[i] - REF_1E11[i]) / REF_1E11[i]).abs())
            .collect();
        println!(
            "rtol {rtol:e}: steps {} rejected {} rel err {:?}  mass {:e}",
            sol.stats.steps,
            sol.stats.rejected,
            rel,
            c.iter().sum::<f64>() - 1.0
        );
    }
    let rk = solve_robertson(1e-6, 1.0).unwrap();
    let y0: Y = Vector::from_components([1.0, 0.0, 0.0], one);
    let atol: Y = Vector::from_components([1e-12, 1e-18, 1e-12], one);
    let dp = solve_ivp(
        robertson,
        (0.0 * s)..(1.0 * s),
        y0,
        DormandPrince54,
        Tolerance { rtol: 1e-6, atol },
    );
    match dp {
        Ok(d) => println!(
            "t ≤ 1 s: RODAS {} steps, Dormand–Prince {} steps",
            rk.stats.steps, d.stats.steps
        ),
        Err(e) => println!(
            "t ≤ 1 s: RODAS {} steps, Dormand–Prince failed: {e}",
            rk.stats.steps
        ),
    }
}

tensorcheck::space!(Phase, 2, Dimensionless);
tensorcheck::chart!(VdpChart on Phase { U: Dimensionless, V: Dimensionless });
type Z = Vector<Dimensionless, VdpChart>;

/// van der Pol with μ = 1000: u' = v, v' = μ (1 − u²) v − u.
fn vdp(_t: Quantity<Time>, z: &Z) -> Vector<tensorcheck::units::si::Frequency, VdpChart> {
    let c = z.components_in(one);
    let mu = 1000.0;
    Vector::from_components([c[1], mu * (1.0 - c[0] * c[0]) * c[1] - c[0]], one / s)
}

fn solve_vdp(
    rtol: f64,
    t_end: f64,
) -> Result<tensorcheck_numerics::ode::Solution<Time, Z>, SolveError<Time>> {
    let z0: Z = Vector::from_components([2.0, 0.0], one);
    let atol: Z = Vector::from_components([rtol, rtol], one);
    solve_ivp(
        vdp,
        (0.0 * s)..(t_end * s),
        z0,
        Rosenbrock {
            jacobian: Jacobian::Numerical,
        },
        Tolerance { rtol, atol },
    )
}

#[test]
#[ignore]
fn measure_steps_and_vdp() {
    let rk = solve_robertson(1e-6, 10.0).unwrap();
    let y0: Y = Vector::from_components([1.0, 0.0, 0.0], one);
    let atol: Y = Vector::from_components([1e-12, 1e-18, 1e-12], one);
    let dp = solve_ivp(
        robertson,
        (0.0 * s)..(10.0 * s),
        y0,
        DormandPrince54,
        Tolerance { rtol: 1e-6, atol },
    );
    println!(
        "t ≤ 10 s: RODAS {} steps; Dormand–Prince {:?}",
        rk.stats.steps,
        dp.map(|d| d.stats.steps)
    );
    let reference = solve_vdp(1e-11, 1000.0).unwrap();
    for rtol in [1e-4, 1e-6, 1e-8] {
        let sol = solve_vdp(rtol, 1000.0).unwrap();
        let (a, b) = (sol.y.components_in(one), reference.y.components_in(one));
        println!(
            "vdp rtol {rtol:e}: steps {} err {:e} {:e}",
            sol.stats.steps,
            (a[0] - b[0]).abs(),
            (a[1] - b[1]).abs()
        );
    }
}

#[test]
fn robertson_matches_the_published_reference() {
    // Measured (2026-10-14): relative error ≈ 0.14 × rtol for rtol = 1e-4, 1e-6,
    // 1e-8 (1533 steps at 1e-8). Allowance 10 × rtol (local error control, 04).
    for rtol in [1e-4, 1e-6, 1e-8] {
        let sol = solve_robertson(rtol, 1e11).unwrap();
        let c = sol.y.components_in(one);
        for i in 0..3 {
            let rel = ((c[i] - REF_1E11[i]) / REF_1E11[i]).abs();
            assert!(
                rel <= 10.0 * rtol,
                "rtol {rtol}: y{} relative error {rel}",
                i + 1
            );
        }
        // y1 + y2 + y3 = 1 is conserved; rounding accumulates at most about
        // steps × ε (measured ≤ 1.1e-15).
        assert!((c.iter().sum::<f64>() - 1.0).abs() <= 1e-12);
    }
}

#[test]
fn rosenbrock_needs_far_fewer_steps_than_an_explicit_method() {
    // Measured on [0, 10] s at rtol 1e-6: RODAS 93 steps, Dormand–Prince 7102.
    let rk = solve_robertson(1e-6, 10.0).unwrap();
    let y0: Y = Vector::from_components([1.0, 0.0, 0.0], one);
    let atol: Y = Vector::from_components([1e-12, 1e-18, 1e-12], one);
    let dp = solve_ivp(
        robertson,
        (0.0 * s)..(10.0 * s),
        y0,
        DormandPrince54,
        Tolerance { rtol: 1e-6, atol },
    )
    .unwrap();
    assert!(
        10 * rk.stats.steps <= dp.stats.steps,
        "{} vs {}",
        rk.stats.steps,
        dp.stats.steps
    );
}

#[test]
fn van_der_pol_converges_to_a_tight_solution() {
    // Self-convergence (no verified published reference was available): compare
    // with rtol = 1e-11. Measured error ≈ 1.3–2.2 × rtol; allowance 10 × rtol.
    let reference = solve_vdp(1e-11, 1000.0).unwrap().y.components_in(one);
    for rtol in [1e-4, 1e-6, 1e-8] {
        let got = solve_vdp(rtol, 1000.0).unwrap().y.components_in(one);
        assert!((got[0] - reference[0]).abs() <= 10.0 * rtol * reference[0].abs().max(1.0));
    }
}
