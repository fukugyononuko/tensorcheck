//! Event detection (issue #21). TC-U29.

use std::cell::Cell;
use std::f64::consts::PI;

use tensorcheck::Quantity;
use tensorcheck::units::si::{Acceleration, Length, Time, Velocity, m, s};
use tensorcheck_numerics::ode::{
    DormandPrince54, Event, Jacobian, Method, Rk4, Rosenbrock, SolveError, Tolerance, solve_ivp,
    solve_ivp_events, solve_ivp_events_with,
};

type Osc = (Quantity<Length>, Quantity<Velocity>);

const ROSENBROCK: Rosenbrock = Rosenbrock {
    jacobian: Jacobian::Numerical,
};

/// x'' = −x with ω = 1/s and x(0) = 1 m: x = cos t crosses zero at π/2 (falling),
/// 3π/2 (rising) and 5π/2 (falling) before 8 s.
fn oscillator(_t: Quantity<Time>, y: &Osc) -> (Quantity<Velocity>, Quantity<Acceleration>) {
    (y.1, -y.0 / (1.0 * (s * s)))
}

fn position(_t: Quantity<Time>, y: &Osc) -> Quantity<Length> {
    y.0
}

/// The events of `event` on the oscillator up to 8 s, as (t, x, v) in s, m, m/s.
fn oscillator_events<M: Method<Jac = Jacobian> + Copy, G>(
    method: M,
    rtol: f64,
    event: Event<G>,
) -> Vec<(f64, f64, f64)>
where
    G: FnMut(Quantity<Time>, &Osc) -> Quantity<Length>,
{
    let tol = Tolerance {
        rtol,
        atol: (rtol * m, rtol * (m / s)),
    };
    let y0 = (1.0 * m, 0.0 * (m / s));
    let sol = solve_ivp_events(oscillator, (0.0 * s)..(8.0 * s), y0, method, tol, event).unwrap();
    assert_eq!(sol.t.value_in(s), 8.0);
    sol.events
        .iter()
        .map(|e| (e.t.value_in(s), e.y.0.value_in(m), e.y.1.value_in(m / s)))
        .collect()
}

/// The largest error of the event times against `want`, after checking that the
/// state at each event is on the solution: x = 0 and v = −sin t.
fn time_error(events: &[(f64, f64, f64)], want: &[f64], v_tol: f64) -> f64 {
    assert_eq!(events.len(), want.len(), "{events:?}");
    let mut worst: f64 = 0.0;
    for (&(t, x, v), &w) in events.iter().zip(want) {
        assert!(x.abs() < 1e-15, "x = {x} at t = {t}");
        assert!((v + w.sin()).abs() <= v_tol, "v = {v} at t = {t}");
        worst = worst.max((t - w).abs());
    }
    worst
}

#[test]
fn adaptive_methods_find_the_zeros_within_the_tolerance() {
    // Global error ≤ 10 × tolerance, as in TC-U23.
    let want = [0.5 * PI, 1.5 * PI, 2.5 * PI];
    for rtol in [1e-6, 1e-9, 1e-12] {
        let e = time_error(
            &oscillator_events(DormandPrince54, rtol, Event::new(position)),
            &want,
            10.0 * rtol,
        );
        assert!(e <= 10.0 * rtol, "Dormand–Prince rtol {rtol}: {e:e}");
    }
    for rtol in [1e-6, 1e-9] {
        let e = time_error(
            &oscillator_events(ROSENBROCK, rtol, Event::new(position)),
            &want,
            10.0 * rtol,
        );
        assert!(e <= 10.0 * rtol, "Rosenbrock rtol {rtol}: {e:e}");
    }
}

#[test]
fn rk4_event_times_are_fourth_order() {
    // The interpolant must not lower the order of RK4: halving h divides the error by 16.
    let want = [0.5 * PI, 1.5 * PI, 2.5 * PI];
    let error = |h: f64| {
        time_error(
            &oscillator_events(Rk4 { h }, 1e-9, Event::new(position)),
            &want,
            1e-4,
        )
    };
    let order = (error(0.1) / error(0.05)).log2();
    assert!((order - 4.0).abs() <= 0.2, "order {order}");
}

#[test]
fn direction_selects_the_crossings() {
    fn check<M: Method<Jac = Jacobian> + Copy>(name: &str, method: M) {
        let times = |event| -> Vec<f64> {
            let ev = oscillator_events(method, 1e-9, event);
            ev.iter().map(|e| e.0).collect()
        };
        let close = |got: Vec<f64>, want: &[f64]| {
            got.len() == want.len() && got.iter().zip(want).all(|(g, w)| (g - w).abs() < 1e-6)
        };
        let g: fn(Quantity<Time>, &Osc) -> Quantity<Length> = position;
        assert!(
            close(times(Event::new(g).rising()), &[1.5 * PI]),
            "{name} rising"
        );
        assert!(
            close(times(Event::new(g).falling()), &[0.5 * PI, 2.5 * PI]),
            "{name} falling"
        );
    }
    check("Rk4", Rk4 { h: 0.01 });
    check("Dormand–Prince", DormandPrince54);
    check("Rosenbrock", ROSENBROCK);
}

#[test]
fn terminal_event_stops_the_solve() {
    fn check<M: Method<Jac = Jacobian> + Copy>(name: &str, method: M) {
        let tol = Tolerance {
            rtol: 1e-9,
            atol: (1e-9 * m, 1e-9 * (m / s)),
        };
        let mut seen = Vec::new();
        let sol = solve_ivp_events_with(
            oscillator,
            (0.0 * s)..(8.0 * s),
            (1.0 * m, 0.0 * (m / s)),
            method,
            tol,
            Event::new(position).terminal(),
            |t, y: &Osc| seen.push((t.value_in(s), y.0.value_in(m))),
        )
        .unwrap();
        assert_eq!(sol.events.len(), 1, "{name}");
        let e = &sol.events[0];
        assert!((e.t.value_in(s) - 0.5 * PI).abs() < 1e-6, "{name}");
        // The solution ends at the event, and on_step sees it last.
        assert_eq!(sol.t.value_in(s), e.t.value_in(s), "{name}");
        assert_eq!(sol.y.0.value_in(m), e.y.0.value_in(m), "{name}");
        assert_eq!(
            *seen.last().unwrap(),
            (e.t.value_in(s), e.y.0.value_in(m)),
            "{name}"
        );
        assert!(seen.iter().all(|&(t, _)| t <= e.t.value_in(s)), "{name}");
        assert_eq!(sol.stats.steps, seen.len(), "{name}");
    }
    check("Rk4", Rk4 { h: 0.01 });
    check("Dormand–Prince", DormandPrince54);
    check("Rosenbrock", ROSENBROCK);
}

// y = t (in metres per second): steps of 0.25 s are exact in binary.
fn ramp(_t: Quantity<Time>, _y: &Quantity<Length>) -> Quantity<Velocity> {
    1.0 * (m / s)
}

fn ramp_tol() -> Tolerance<Quantity<Length>> {
    Tolerance {
        rtol: 1e-9,
        atol: 1e-12 * m,
    }
}

#[test]
fn zero_at_a_step_end_is_found_once_and_zero_at_the_start_is_not_an_event() {
    let at_half = Event::new(|_t: Quantity<Time>, y: &Quantity<Length>| *y - 0.5 * m);
    let sol = solve_ivp_events(
        ramp,
        (0.0 * s)..(1.0 * s),
        0.0 * m,
        Rk4 { h: 0.25 },
        ramp_tol(),
        at_half,
    )
    .unwrap();
    let times: Vec<f64> = sol.events.iter().map(|e| e.t.value_in(s)).collect();
    assert_eq!(times, [0.5]);

    let at_start = Event::new(|_t: Quantity<Time>, y: &Quantity<Length>| *y);
    let sol = solve_ivp_events(
        ramp,
        (0.0 * s)..(1.0 * s),
        0.0 * m,
        Rk4 { h: 0.25 },
        ramp_tol(),
        at_start,
    )
    .unwrap();
    assert!(sol.events.is_empty());
}

#[test]
fn tiny_values_of_g_are_compared_by_sign_not_by_product() {
    // PR #41 review (Codex): g(0) g(0.75) = −0.5e-300 × 0.25e-300 underflows to −0.
    for scale in [1e-300, 1e300] {
        let g = move |_t: Quantity<Time>, y: &Quantity<Length>| (*y - 0.5 * m) * scale;
        let sol = solve_ivp_events(
            ramp,
            (0.0 * s)..(1.0 * s),
            0.0 * m,
            Rk4 { h: 0.75 },
            ramp_tol(),
            Event::new(g),
        )
        .unwrap();
        assert_eq!(sol.events.len(), 1, "scale {scale:e}");
        assert!(
            (sol.events[0].t.value_in(s) - 0.5).abs() < 1e-12,
            "scale {scale:e}"
        );
    }
}

#[test]
fn backward_in_time_the_signs_are_read_along_the_solve() {
    // From t = 1 s back to 0: y − 0.3 m goes from positive to negative.
    let g = |_t: Quantity<Time>, y: &Quantity<Length>| *y - 0.3 * m;
    for (name, falling, rising) in [
        (
            "Rk4",
            solve_ivp_events(
                ramp,
                (1.0 * s)..(0.0 * s),
                1.0 * m,
                Rk4 { h: 0.25 },
                ramp_tol(),
                Event::new(g).falling(),
            ),
            solve_ivp_events(
                ramp,
                (1.0 * s)..(0.0 * s),
                1.0 * m,
                Rk4 { h: 0.25 },
                ramp_tol(),
                Event::new(g).rising(),
            ),
        ),
        (
            "Dormand–Prince",
            solve_ivp_events(
                ramp,
                (1.0 * s)..(0.0 * s),
                1.0 * m,
                DormandPrince54,
                ramp_tol(),
                Event::new(g).falling(),
            ),
            solve_ivp_events(
                ramp,
                (1.0 * s)..(0.0 * s),
                1.0 * m,
                DormandPrince54,
                ramp_tol(),
                Event::new(g).rising(),
            ),
        ),
        (
            "Rosenbrock",
            solve_ivp_events(
                ramp,
                (1.0 * s)..(0.0 * s),
                1.0 * m,
                ROSENBROCK,
                ramp_tol(),
                Event::new(g).falling(),
            ),
            solve_ivp_events(
                ramp,
                (1.0 * s)..(0.0 * s),
                1.0 * m,
                ROSENBROCK,
                ramp_tol(),
                Event::new(g).rising(),
            ),
        ),
    ] {
        let falling = falling.unwrap();
        assert_eq!(falling.events.len(), 1, "{name}");
        assert!(
            (falling.events[0].t.value_in(s) - 0.3).abs() < 1e-12,
            "{name}"
        );
        assert!(rising.unwrap().events.is_empty(), "{name}");
    }
}

#[test]
fn rk4_evaluates_f_once_more_only_in_a_step_with_an_event() {
    let solve = |level: f64| {
        let g = move |_t: Quantity<Time>, y: &Quantity<Length>| *y - level * m;
        solve_ivp_events(
            ramp,
            (0.0 * s)..(1.0 * s),
            0.0 * m,
            Rk4 { h: 0.25 },
            ramp_tol(),
            Event::new(g),
        )
        .unwrap()
    };
    // Four steps of four evaluations each.
    assert_eq!(solve(2.0).stats.evaluations, 16);
    assert_eq!(solve(0.3).stats.evaluations, 17);
}

#[test]
fn adaptive_methods_evaluate_f_no_more_with_events() {
    // Dormand–Prince and RODAS already have f at the end of the step.
    fn check<M: Method<Jac = Jacobian> + Copy>(name: &str, method: M) {
        let tol = Tolerance {
            rtol: 1e-9,
            atol: (1e-9 * m, 1e-9 * (m / s)),
        };
        let y0 = (1.0 * m, 0.0 * (m / s));
        let plain = solve_ivp(oscillator, (0.0 * s)..(8.0 * s), y0, method, tol).unwrap();
        let with = solve_ivp_events(
            oscillator,
            (0.0 * s)..(8.0 * s),
            y0,
            method,
            tol,
            Event::new(position),
        )
        .unwrap();
        assert_eq!(with.events.len(), 3, "{name}");
        assert_eq!(with.stats, plain.stats, "{name}");
    }
    check("Dormand–Prince", DormandPrince54);
    check("Rosenbrock", ROSENBROCK);
}

#[test]
fn non_finite_g_is_a_failure_with_its_time() {
    let solve = |g: &dyn Fn(f64, f64) -> f64| {
        let event = Event::new(|t: Quantity<Time>, y: &Quantity<Length>| {
            g(t.value_in(s), y.value_in(m)) * m
        });
        solve_ivp_events(
            ramp,
            (0.0 * s)..(1.0 * s),
            0.0 * m,
            Rk4 { h: 0.25 },
            ramp_tol(),
            event,
        )
    };
    // At the start, at the end of the second step, and only inside the step with the
    // root (0.4 s < t < 0.45 s, root at 0.42 s), where the interpolant is searched.
    let at_start = solve(&|_, _| f64::NAN);
    assert!(matches!(at_start, Err(SolveError::NonFinite { t }) if t.value_in(s) == 0.0));
    let at_step = solve(&|t, _| if t > 0.3 { f64::INFINITY } else { 1.0 });
    assert!(matches!(at_step, Err(SolveError::NonFinite { t }) if t.value_in(s) == 0.5));
    let inside = solve(&|t, y| {
        if t > 0.4 && t < 0.45 {
            f64::NAN
        } else {
            y - 0.42
        }
    });
    let t = inside.unwrap_err().time().unwrap().value_in(s);
    assert!(t > 0.4 && t < 0.45, "t = {t}");
}

// PR #41 review (Codex): f evaluated only for the interpolation must be finite too.
type Pair = (Quantity<Length>, Quantity<Length>);

#[test]
fn non_finite_f_at_the_step_end_is_a_failure_for_rk4() {
    // x' = x, z' = NaN only for 2.70 m ≤ x < 2.72 m: the RK4 stages miss the band,
    // the extra f(t1, y1) for the event (x(1) ≈ 2.708 m) does not.
    let f = |_t: Quantity<Time>, y: &Pair| {
        let x = y.0.value_in(m);
        let dz = if (2.70..2.72).contains(&x) {
            f64::NAN
        } else {
            0.0
        };
        (y.0 / (1.0 * s), dz * (m / s))
    };
    let tol = Tolerance {
        rtol: 1e-9,
        atol: (1e-12 * m, 1e-12 * m),
    };
    let event = Event::new(|_t: Quantity<Time>, y: &Pair| y.0 - 2.0 * m).terminal();
    let r = solve_ivp_events(
        f,
        (0.0 * s)..(1.0 * s),
        (1.0 * m, 0.0 * m),
        Rk4 { h: 1.0 },
        tol,
        event,
    );
    assert!(
        matches!(r, Err(SolveError::NonFinite { t }) if t.value_in(s) == 1.0),
        "{r:?}"
    );
}

#[test]
fn non_finite_f_at_the_step_end_is_a_failure_for_rosenbrock() {
    // First find the first accepted point (t1, y1) without NaN; then make z' NaN
    // exactly there, so only the f(t1, y1) computed after acceptance is not finite.
    let tol = Tolerance {
        rtol: 1e-6,
        atol: (1e-9 * m, 1e-9 * m),
    };
    let clean = |_t: Quantity<Time>, y: &Pair| (y.0 / (1.0 * s), 0.0 * (m / s));
    let mut first = None;
    solve_ivp_events_with(
        clean,
        (0.0 * s)..(1.0 * s),
        (1.0 * m, 0.0 * m),
        ROSENBROCK,
        tol,
        Event::new(|_t: Quantity<Time>, y: &Pair| y.0 - 2.0 * m),
        |t, y: &Pair| {
            first.get_or_insert((t.value_in(s), y.0.value_in(m)));
        },
    )
    .unwrap();
    let (t1, x1) = first.unwrap();
    // The event must fall in the first step for the interpolation to use f(t1, y1).
    let level = 1.0 + 0.5 * (x1 - 1.0);
    let f = move |t: Quantity<Time>, y: &Pair| {
        let hit = t.value_in(s) == t1 && y.0.value_in(m) == x1;
        (y.0 / (1.0 * s), if hit { f64::NAN } else { 0.0 } * (m / s))
    };
    let event = Event::new(move |_t: Quantity<Time>, y: &Pair| y.0 - level * m).terminal();
    let r = solve_ivp_events(
        f,
        (0.0 * s)..(1.0 * s),
        (1.0 * m, 0.0 * m),
        ROSENBROCK,
        tol,
        event,
    );
    assert!(
        matches!(r, Err(SolveError::NonFinite { t }) if t.value_in(s) == t1),
        "{r:?}"
    );
}

#[test]
fn overflow_inside_the_interpolant_is_a_failure() {
    // PR #41 third review (Codex): the RK4 stages (1, 0, 0, 0) and f(t1, y1) = 1.5e308
    // are finite, but the Hermite slope term overflows at θ = 0.5.
    let f = |t: Quantity<Time>, y: &Quantity<Length>| {
        let (t, y) = (t.value_in(s), y.value_in(m));
        let v = if t == 0.0 {
            1.0
        } else if t == 10.0 && y != 0.0 {
            1.5e308
        } else {
            0.0
        };
        v * (m / s)
    };
    let event = Event::new(|t: Quantity<Time>, _y: &Quantity<Length>| t - 5.0 * s).terminal();
    let r = solve_ivp_events(
        f,
        (0.0 * s)..(10.0 * s),
        0.0 * m,
        Rk4 { h: 10.0 },
        ramp_tol(),
        event,
    );
    assert!(matches!(r, Err(SolveError::NonFinite { .. })), "{r:?}");
}

#[test]
fn a_steep_g_is_still_located() {
    // PR #41 fourth review (Codex): y = t, g = y¹⁰⁰ − 0.5¹⁰⁰ (in m), one RK4 step.
    let c = 0.5f64.powi(100);
    let g = move |_t: Quantity<Time>, y: &Quantity<Length>| (y.value_in(m).powi(100) - c) * m;
    let sol = solve_ivp_events(
        ramp,
        (0.0 * s)..(1.0 * s),
        0.0 * m,
        Rk4 { h: 1.0 },
        ramp_tol(),
        Event::new(g).terminal(),
    )
    .unwrap();
    assert!(
        (sol.t.value_in(s) - 0.5).abs() < 1e-12,
        "t = {}",
        sol.t.value_in(s)
    );
}

#[test]
fn a_discontinuous_g_gives_the_point_where_the_sign_jumps() {
    // Issue #42 (PR #41 sixth review, Codex): g has no zero, only a jump of sign. If g
    // is finite wherever it is evaluated, the event is either a point where g was
    // evaluated exactly 0, or the end of the closed bracket on the side of the sign
    // after the change, next to the jump. If the search evaluates g exactly at a pole,
    // the solve fails with NonFinite at that time (PR #47 reviews, Codex).
    // Returns whether the solve succeeded.
    fn check(name: &str, g: fn(f64) -> f64, jump: f64) -> bool {
        let r = solve_ivp_events(
            ramp,
            (0.0 * s)..(1.0 * s),
            0.0 * m,
            Rk4 { h: 1.0 },
            ramp_tol(),
            Event::new(move |t: Quantity<Time>, _y: &Quantity<Length>| g(t.value_in(s)) * m)
                .terminal(),
        );
        let sol = match r {
            Ok(sol) => sol,
            Err(SolveError::NonFinite { t }) => {
                let t = t.value_in(s);
                assert!(!g(t).is_finite(), "{name}: NonFinite at {t}, g = {}", g(t));
                return false;
            }
            Err(e) => panic!("{name}: {e}"),
        };
        assert_eq!(sol.events.len(), 1, "{name}");
        let b = sol.t.value_in(s);
        if g(b) == 0.0 {
            return true;
        }
        let after = g(1.0).signum();
        assert_eq!(g(b).signum(), after, "{name}: g({b}) = {}", g(b));
        assert_eq!(
            g(b - 4.0 * f64::EPSILON).signum(),
            -after,
            "{name}: bracket not closed at {b}"
        );
        assert!((b - jump).abs() <= 4.0 * f64::EPSILON, "{name}: t = {b}");
        true
    }
    assert!(check(
        "step at 0.5",
        |t| if t < 0.5 { -1.0 } else { 1.0 },
        0.5
    ));
    assert!(check(
        "falling step at 0.3",
        |t| if t < 0.3 { 2.0 } else { -3.0 },
        0.3,
    ));
    // An exact zero at 0.5 is the first secant point, so it is returned although the
    // sign jumps at 0.75 (PR #47 third review, Codex).
    assert!(check(
        "isolated zero before a jump",
        |t| {
            if t == 0.5 {
                0.0
            } else if t < 0.75 {
                -1.0
            } else {
                1.0
            }
        },
        0.75
    ));
    // Continuous, but g(1) = 0 ends the search at once (PR #47 second review, Codex).
    assert!(check(
        "zero plateau from 0.2",
        |t| if t < 0.2 { t - 0.2 } else { 0.0 },
        0.2
    ));
    // The search does not land on the pole at 0.7, but lands exactly on the one at
    // 0.5 (g = 2 and −2 at the ends, so the first secant point is 0.5).
    assert!(check("pole at 0.7", |t| 1.0 / (0.7 - t), 0.7));
    assert!(!check("pole at 0.5", |t| 1.0 / (0.5 - t), 0.5));

    // Backward in time, the bracket end is on the side of the step end, which is the
    // earlier time (PR #47 fifth review, Codex).
    let g = |t: f64| if t < 0.5 { -1.0 } else { 1.0 };
    let sol = solve_ivp_events(
        ramp,
        (1.0 * s)..(0.0 * s),
        1.0 * m,
        Rk4 { h: 1.0 },
        ramp_tol(),
        Event::new(move |t: Quantity<Time>, _y: &Quantity<Length>| g(t.value_in(s)) * m).terminal(),
    )
    .unwrap();
    let b = sol.t.value_in(s);
    assert_eq!(
        g(b),
        g(0.0),
        "backward: g({b}) is not on the side of the step end"
    );
    assert_eq!(
        g(b + 4.0 * f64::EPSILON),
        -g(0.0),
        "backward: bracket not closed at {b}"
    );
    assert!((b - 0.5).abs() <= 4.0 * f64::EPSILON, "backward: t = {b}");
}

#[test]
fn at_large_times_only_the_bracket_in_theta_is_bounded() {
    // Issue #48 (PR #47 sixth review, Codex): g is evaluated at t_n + θ h rounded to
    // f64. With t_n = 1e12 and h = ±1, θ h is exact and rounding is monotone, so the
    // two ends of a bracket in θ at most 4ε wide have times that are equal or
    // adjacent f64 values; they differ in sign, so they are adjacent. The time one
    // f64 before the event, in the direction of integration, has the other sign.
    let t0 = 1e12;
    let jump = t0 + 0.5;
    let g = move |t: f64| if t < jump { -1.0 } else { 1.0 };
    for (name, span, before) in [
        ("forward", t0..(t0 + 1.0), f64::next_down as fn(f64) -> f64),
        ("backward", (t0 + 1.0)..t0, f64::next_up as fn(f64) -> f64),
    ] {
        let sol = solve_ivp_events(
            ramp,
            (span.start * s)..(span.end * s),
            0.0 * m,
            Rk4 { h: 1.0 },
            ramp_tol(),
            Event::new(move |t: Quantity<Time>, _y: &Quantity<Length>| g(t.value_in(s)) * m)
                .terminal(),
        )
        .unwrap();
        let b = sol.t.value_in(s);
        let a = before(b);
        assert_eq!(
            g(b),
            g(span.end),
            "{name}: g({b}) is not on the side of the step end"
        );
        assert_eq!(
            g(a),
            -g(span.end),
            "{name}: g({a}) next to {b} has the same sign"
        );
        // The times are further apart than 4ε times the step, so the documentation
        // must not bound the width in time.
        assert!((b - a).abs() > 4.0 * f64::EPSILON, "{name}: {b} − {a}");
    }
}

#[test]
fn a_tiny_g_beyond_the_root_does_not_end_the_search() {
    // PR #41 fifth review (Codex): g(1) = 5e-324; halving that weight gives 0, which
    // was taken for a root at t = 1. The only root is t = 0.9.
    let g = |t: Quantity<Time>, _y: &Quantity<Length>| {
        let x = t.value_in(s) - 0.9;
        (x * (-7420.0 * x.max(0.0)).exp()) * m
    };
    let sol = solve_ivp_events(
        ramp,
        (0.0 * s)..(1.0 * s),
        0.0 * m,
        Rk4 { h: 1.0 },
        ramp_tol(),
        Event::new(g).terminal(),
    )
    .unwrap();
    assert!(
        (sol.t.value_in(s) - 0.9).abs() < 1e-12,
        "t = {}",
        sol.t.value_in(s)
    );
}

#[test]
fn rejected_input_calls_neither_f_nor_g() {
    let (f_calls, g_calls) = (Cell::new(0), Cell::new(0));
    let f = |t: Quantity<Time>, y: &Quantity<Length>| {
        f_calls.set(f_calls.get() + 1);
        ramp(t, y)
    };
    let event = || {
        Event::new(|_t: Quantity<Time>, y: &Quantity<Length>| {
            g_calls.set(g_calls.get() + 1);
            *y
        })
    };
    let bad_tol = Tolerance {
        rtol: -1.0,
        atol: 1e-12 * m,
    };
    let r = solve_ivp_events(
        f,
        (0.0 * s)..(1.0 * s),
        0.0 * m,
        DormandPrince54,
        bad_tol,
        event(),
    );
    assert!(matches!(r, Err(SolveError::InvalidTolerance)));
    let r = solve_ivp_events(
        f,
        (0.0 * s)..(1.0 * s),
        0.0 * m,
        DormandPrince54.max_steps(0),
        ramp_tol(),
        event(),
    );
    assert!(matches!(r, Err(SolveError::MaxSteps { t }) if t.value_in(s) == 0.0));
    assert_eq!((f_calls.get(), g_calls.get()), (0, 0));
}
