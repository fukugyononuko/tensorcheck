#!/bin/sh
# CI と同じ検査を手元で実行する。1 つでも失敗したら、その時点で止まる（コミットの前に実行する）。
set -eu
python3 scripts/check_workflows.py
python3 scripts/check_layers.py
if python3 scripts/check_layers.py scripts/fixtures/layers_bad >/dev/null 2>&1; then
  echo "check_layers が、わざと壊した例を見逃した" >&2
  exit 1
fi
python3 scripts/check_deps.py
python3 scripts/check_english.py
python3 scripts/check_math.py
python3 scripts/check_privacy.py
cargo fmt --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
echo "ci_local: すべて通りました"
