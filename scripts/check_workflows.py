"""GitHub Actions のワークフローが決まりを守っているかを検査する（ホームページの TC-N10 と同じ決まり）。

- 外部の Action はコミット SHA（40 桁）で固定する
- ワークフロー全体の既定の権限は contents: read だけ
- pull_request_target を使わない
"""

import re
import sys
from pathlib import Path

USES = re.compile(r"^\s*(?:-\s*)?uses:\s*([^\s#]+)")


def main() -> int:
    problems = []
    for path in sorted(Path(".github/workflows").glob("*.y*ml")):
        lines = path.read_text(encoding="utf-8").split("\n")
        for i, line in enumerate(lines, 1):
            m = USES.match(line)
            if m and not m[1].startswith("./") and not re.search(r"@[0-9a-f]{40}$", m[1]):
                problems.append(f"{path}:{i} コミット SHA で固定されていない Action: {m[1]}")
            if re.match(r"^\s*pull_request_target\s*:?", line):
                problems.append(f"{path}:{i} pull_request_target は使わない")

        starts = [i for i, l in enumerate(lines) if re.match(r"^permissions:\s*$", l)]
        if not starts:
            problems.append(f"{path}: ワークフロー全体の permissions がない")
            continue
        block = []
        for l in lines[starts[0] + 1 :]:
            if not re.match(r"^\s+\S", l):
                break
            block.append(l.strip())
        if block != ["contents: read"]:
            problems.append(f"{path}: ワークフロー全体の permissions は contents: read だけにする（実際: {', '.join(block) or 'なし'}）")

    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    print("check_workflows: 問題はありません")
    return 0


if __name__ == "__main__":
    sys.exit(main())
