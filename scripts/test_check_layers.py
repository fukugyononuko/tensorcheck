"""check_layers.py 自体のテスト（TC-A01）。

- fixtures/layers_bad：`VIOLATION` と書いたファイルは、1 つ残らず違反として見つかること。書いていないファイルは見つからないこと。
  すり抜けの書き方（まとめた use、モジュールを use してからの相対の呼び出し、self・super、別名、完全なパス、
  再公開、glob、マクロの中の $crate、先頭の ::）を 1 ファイルに 1 つずつ置いている（Codex のレビュー 2026-10-14、指摘 5）。
- fixtures/layers_good：正しい書き方（同じファイル、下の層、契約、コメントと文字列の中のパス）で、違反が 0 件であること。
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from check_layers import check  # noqa: E402

FIXTURES = Path(__file__).parent / "fixtures"


def main() -> int:
    failures = []
    bad_root = FIXTURES / "layers_bad"
    found = {line.split(": ", 1)[0] for line in check(bad_root)}
    for path in sorted(bad_root.rglob("*.rs")):
        marked = "VIOLATION" in path.read_text(encoding="utf-8")
        if marked and str(path) not in found:
            failures.append(f"見逃した: {path}")
        if not marked and str(path) in found:
            failures.append(f"違反でないのに見つけた: {path}")
    good = check(FIXTURES / "layers_good")
    failures += [f"正しい例を違反とした: {line}" for line in good]
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print(f"test_check_layers: 違反の例 {len(found)} ファイルをすべて見つけ、正しい例は 0 件でした")
    return 0


if __name__ == "__main__":
    sys.exit(main())
