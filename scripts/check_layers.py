"""呼び出しの規則（01 §1.7、DD-01）と unsafe の禁止を検査する（TC-A01、TC-A02）。

- manager/ は engine/・access/・utility/ だけを使う。ほかの manager/ のファイルを使わない。
- engine/ は access/・utility/ だけを使う。ほかの engine/ のファイルを使わない。
- access/ は utility/ だけを使う。utility/ は utility/ だけを使う。
- tensorcheck::internal を使ってよいのは tensorcheck-numerics の src/manager/ だけ。
- 両方のクレートの lib.rs に #![forbid(unsafe_code)] がある。

引数にディレクトリを渡すと、そこを crates/ の代わりに検査する（検査自体のテスト用）。
"""

import re
import sys
from pathlib import Path

LAYERS = ("manager", "engine", "access", "utility")
ALLOWED = {
    "manager": {"engine", "access", "utility"},
    "engine": {"access", "utility"},
    "access": {"utility"},
    "utility": {"utility"},
}
USE_CRATE = re.compile(r"crate::(manager|engine|access|utility)::(\w+)")


def check(root: Path) -> list[str]:
    problems = []
    for crate in sorted(p for p in root.iterdir() if p.is_dir()):
        src = crate / "src"
        lib = src / "lib.rs"
        if lib.exists() and "#![forbid(unsafe_code)]" not in lib.read_text(encoding="utf-8"):
            problems.append(f"{lib}: #![forbid(unsafe_code)] がない")
        for path in sorted(crate.rglob("*.rs")):
            text = path.read_text(encoding="utf-8")
            rel = path.relative_to(crate)
            parts = rel.parts
            layer = parts[1] if len(parts) > 2 and parts[0] == "src" and parts[1] in LAYERS else None
            module = Path(parts[-1]).stem if layer else None
            if "tensorcheck::internal" in text:
                if not (crate.name == "tensorcheck-numerics" and layer == "manager"):
                    problems.append(f"{path}: tensorcheck::internal は numerics の manager だけが使える")
            if layer is None:
                continue
            for target, target_module in USE_CRATE.findall(text):
                if target not in ALLOWED[layer]:
                    problems.append(f"{path}: {layer} から {target} を使っている（{target}::{target_module}）")
                elif target == layer and layer in ("manager", "engine") and target_module not in (module, "mod"):
                    problems.append(f"{path}: {layer} どうしの呼び出し（{target}::{target_module}）")
    return problems


def main() -> int:
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("crates")
    problems = check(root)
    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    print(f"check_layers: 問題はありません（{root}）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
