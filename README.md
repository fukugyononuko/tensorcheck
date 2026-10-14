# tensorcheck

A Rust library that checks physical dimensions and tensor structure (rank and shape, upper/lower indices, coordinate systems, and symmetry) at compile time. It works for tensors on coordinate spaces (space, spacetime) and on field spaces alike, and provides numerical solvers (ODEs, quadrature, numerical differentiation) that keep those checks intact.

Mistakes such as adding a force to an acceleration, contracting two lower indices, mixing polar and Cartesian components, or transforming Christoffel symbols as if they were a tensor are **compile errors** with a readable headline. At run time the checks cost nothing: contractions run as fast as hand-written `f64` loops (see [Performance](#performance)).

Design documents (in Japanese) are in [docs/](docs/).

- Author: [fukugyononuko](https://fukugyononuko.github.io/)

## A taste

```rust
use tensorcheck::labels::a;
use tensorcheck::spaces::plane::{Polar, PolarMetric, Theta};
use tensorcheck::units::si::{Frequency, Velocity, m, rad, s};
use tensorcheck::{Point, Quantity, Vector, slot};

// A velocity in polar coordinates: V^r in m/s, V^θ in 1/s.
let p = Point::<Polar>::new((2.0 * m, 0.5 * rad));
let v: Vector<Velocity, Polar> = Vector::from_parts((1.5 * (m / s), 0.25 * (rad / s)));
let v_theta: Quantity<Frequency> = v.at((Theta,));

// g_ab V^a V^b, written with the Einstein convention.
let v2 = (v.lower(slot::_0, &PolarMetric, &p).idx(a) * v.idx(a)).into_tensor();
assert!((v2.value_in(m * m / (s * s)) - 2.5).abs() < 1e-14);

// let wrong = (v.idx(a) * v.idx(a)).into_tensor();
//   error: cannot contract index `a`: it is upper in both places
```

## What it covers

| | |
| --- | --- |
| Dimensions and units | SI (7 base dimensions), natural units (`ħ = c = 1`, GeV, reduced Planck mass), user-defined systems and units (`unit!`) |
| Tensors | any rank; each index has a space, a chart and a variance; components may have different dimensions in curvilinear charts; symmetric and antisymmetric rank-2 tensors |
| Expressions | Einstein-convention contraction with labels (`t.idx((a, b)) * v.idx(b)`), traces, outer products, raising and lowering with a metric tensor, chart transforms |
| Geometry | Christoffel symbols, Riemann, Ricci, Einstein tensors, scalar curvature and covariant derivatives from the metric tensor alone, with exact derivatives by dual numbers |
| Field spaces | the same machinery for field-space tensors and mixed objects such as `∂_μ φ^I` |
| Numerics (`tensorcheck-numerics`) | ODEs (RK4, Dormand–Prince 5(4), Rosenbrock RODAS for stiff systems), quadrature (trapezoid, Simpson, Romberg, Gauss–Legendre, adaptive Gauss–Kronrod), numerical derivatives — all keeping dimensions |

Shipped spaces: the plane (Cartesian, polar), 3D Euclidean space (Cartesian, spherical, cylindrical), Minkowski spacetime `(+, −, −, −)`, the two-sphere (angular, stereographic).

## Examples

Each example checks its result against an analytic or published value and fails if it does not match.

| Example | Use case |
| --- | --- |
| `cargo run -p tensorcheck --example uc01_newton` | F = m a with checked dimensions |
| `cargo run -p tensorcheck --example uc02_polar` | speed² in polar coordinates |
| `cargo run -p tensorcheck --example uc03_sigma_model` | `G_IJ(φ) ∂_μφ^I ∂^μφ^J` on Minkowski spacetime |
| `cargo run -p tensorcheck --example uc04_transform` | polar to Cartesian |
| `cargo run -p tensorcheck --example uc09_sphere_curvature` | curvature of a sphere in two charts |
| `cargo run -p tensorcheck-numerics --example uc07_robertson` | a stiff ODE to t = 1e11 s |
| `cargo run -p tensorcheck-numerics --example uc08_work` | work of a spring, typed as energy |
| `cargo run -p tensorcheck-numerics --example uc09_geodesic` | geodesics on a sphere |

## Tips

- Unit names such as `m` and `s` are ordinary identifiers. If they clash with your variable names, import the module instead: `use tensorcheck::units::si as u;` and write `2.0 * u::m`.
- `Quantity::new(x, unit)` and `from_components` need the dimension from a type annotation. Without one, use `Quantity::of(x, unit)` or `x * unit`.
- Write metric tensors and tensor fields generically over `N: Scalar`; the library evaluates them with dual numbers to get derivatives.
- Higher-order ODEs are solved as first-order systems: for `x'' = a(x)` use the state `(x, v)`.

## Performance

`cargo bench -p tensorcheck --bench overhead` and `cargo bench -p tensorcheck-numerics --bench overhead` compare with the same computation written by hand on `f64` arrays (ratio = checked / hand-written, release build):

| Operation | Ratio |
| --- | --- |
| `g_ab V^b` (4D) | 1.00 |
| `R^a_bcd T^cd` (4D) | 0.98 |
| lowering an index (4D) | 0.39 |
| RK4 step, 2-component state, trivial right-hand side (solver overhead) | 2.19 |

Contraction plans are computed at compile time for each combination of types, so the index bookkeeping disappears at run time.

## Limitations

- Symmetry markers are for rank-2 tensors only; Riemann symmetries are not encoded in the type.
- Connections are Levi-Civita in a coordinate basis (no torsion, no orthonormal frames yet).
- `transform` needs every index in the source chart.
- The Rosenbrock solver uses a numerical Jacobian; user-supplied Jacobians are planned.
- ODE solves stop after 100,000 steps.
- Requires `std` (no `no_std` yet).

## How it differs from existing crates

As of 2026-10-12 (details: [docs/01_overview.md](docs/01_overview.md) §1.1, in Japanese):

| Area | Crates | Difference |
| --- | --- | --- |
| Dimensions | uom, dimensioned, diman | Scalar dimensions only; no tensor structure. diman requires unstable Rust |
| Tensors | differential-geometry | Closest match (dimension $`n`$ and coordinate systems as types), but no physical units, symmetry, or Einstein-convention contraction; not updated since 2019 |
| ODEs | diffsol, russell_ode, ode_solvers | State is a plain array of numbers, so dimension checks stop at the solver boundary |
| Quadrature | bilby, integrate, quad-rs | Integrands return plain numbers |

We found no crate that combines dimensions and tensor structure in one type, and no numerical solver that keeps dimensions through the computation.

## Development

Requires Rust (version pinned in `rust-toolchain.toml`) and Python 3 (for the check scripts).

```sh
cargo fmt --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo doc --workspace --no-deps
python3 scripts/check_layers.py
python3 scripts/check_deps.py
python3 scripts/check_english.py
python3 scripts/check_workflows.py
python3 scripts/check_privacy.py
python3 scripts/check_math.py
```

## License

Licensed under either of [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE), at your option.
