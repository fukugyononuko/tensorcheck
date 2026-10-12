"""文書の数式の書き方を検査する（docs/README.md の「表記」）。

GitHub は、`$...$` の形の数式を、隣に全角の記号があると数式として表示しない。
そのため、文中の数式は $`...`$ の形、独立した行の数式は ```math のブロックで書く。
コードブロックとコードの中を除いて、それ以外の `$` が残っていたら失敗にする。
"""

import re
import sys
from pathlib import Path

FILES = sorted(Path("docs").rglob("*.md")) + [Path("README.md")]
CODE_SPAN = re.compile(r"(`+)(.+?)\1")
INLINE_MATH = re.compile(r"\$`[^`]*`\$")


def main() -> int:
    problems = []
    for path in FILES:
        in_block = False
        for no, line in enumerate(path.read_text(encoding="utf-8").split("\n"), 1):
            if line.lstrip().startswith("```"):
                in_block = not in_block
                continue
            if in_block:
                continue
            rest = INLINE_MATH.sub("", line)
            rest = CODE_SPAN.sub("", rest)
            if "$" in rest:
                problems.append(f"{path}:{no}: 文中の数式は $`...`$ の形で書く")
    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    print(f"check_math: 問題はありません（{len(FILES)} ファイル）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
