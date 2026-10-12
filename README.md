# tensorcheck

物理量の次元と、テンソルの構造（階数と形・添字の上下・座標系・対称性）を、Rust のコンパイル時に検査するライブラリです。座標（空間・時空）のテンソルにも、場空間のテンソルにも使えます。

**開発中です。** まだ使える API はありません。

- 制作：[副業のぬこ](https://fukugyononuko.github.io/)

## 開発

Rust（`rust-toolchain.toml` で版を固定）と Python 3（検査スクリプト）を使います。

```sh
cargo fmt --check
cargo clippy --all-targets -- -D warnings
cargo test
cargo doc --no-deps
python3 scripts/check_workflows.py
python3 scripts/check_privacy.py
```

## ライセンス

MIT または Apache-2.0 のどちらかを選べます（[LICENSE-MIT](LICENSE-MIT)、[LICENSE-APACHE](LICENSE-APACHE)）。
