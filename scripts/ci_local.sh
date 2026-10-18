#!/bin/sh
# CI と同じ検査を手元で実行する。1 つでも失敗したら、その時点で止まる（コミットの前に実行する）。
set -eu
python3 scripts/check_workflows.py
python3 scripts/check_layers.py
python3 scripts/test_check_layers.py
python3 scripts/check_deps.py
python3 scripts/check_english.py
python3 scripts/check_math.py
python3 scripts/check_privacy.py
python3 scripts/test_check_privacy.py
cargo fmt --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
for e in uc01_newton uc02_polar uc03_sigma_model uc04_transform uc09_sphere_curvature; do
  cargo run --release -q -p tensorcheck --example "$e" >/dev/null
done
for e in uc07_robertson uc08_work uc09_geodesic grav01_schwarzschild_curvature grav02_kepler grav03_perihelion grav04_light_bending; do
  cargo run --release -q -p tensorcheck-numerics --example "$e" >/dev/null
done
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
echo "ci_local: すべて通りました"
