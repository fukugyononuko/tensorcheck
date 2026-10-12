# tensorcheck

A Rust library that checks physical dimensions and tensor structure (rank and shape, upper/lower indices, coordinate systems, and symmetry) at compile time. It works for tensors on coordinate spaces (space, spacetime) and on field spaces alike, and provides numerical solvers (ODEs, quadrature, numerical differentiation) that keep those checks intact.

**Work in progress.** There is no usable API yet. Design documents (in Japanese) are in [docs/](docs/).

- Author: [fukugyononuko](https://fukugyononuko.github.io/)

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
cargo clippy --all-targets -- -D warnings
cargo test
cargo doc --no-deps
python3 scripts/check_workflows.py
python3 scripts/check_privacy.py
python3 scripts/check_math.py
```

## License

Licensed under either of [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE), at your option.
