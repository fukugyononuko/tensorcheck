"""実行時の依存を検査する（TC-A03、TC-A04、NFR-01）。

- tensorcheck は実行時に何にも依存しない。
- tensorcheck-numerics は実行時に tensorcheck だけに依存する。
開発だけで使う依存（dev-dependencies）は対象外。
"""

import json
import subprocess
import sys

EXPECTED = {"tensorcheck": set(), "tensorcheck-numerics": {"tensorcheck"}}


def main() -> int:
    meta = json.loads(subprocess.run(
        ["cargo", "metadata", "--format-version", "1", "--no-deps"],
        check=True, capture_output=True, text=True).stdout)
    problems = []
    for pkg in meta["packages"]:
        runtime = {d["name"] for d in pkg["dependencies"] if d["kind"] is None}
        want = EXPECTED.get(pkg["name"])
        if want is None:
            problems.append(f"想定していないパッケージ: {pkg['name']}")
        elif runtime != want:
            problems.append(f"{pkg['name']} の実行時の依存が {sorted(runtime)}（期待: {sorted(want)}）")
    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    print("check_deps: 問題はありません")
    return 0


if __name__ == "__main__":
    sys.exit(main())
