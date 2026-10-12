"""crates/ の下と README.md が英語で書かれているかを検査する（TC-N05、NFR-07）。

ひらがな・カタカナ・漢字（CJK 統合漢字）が含まれていたら失敗にする。
"""

import re
import sys
from pathlib import Path

CJK = re.compile(r"[぀-ヿ㐀-䶿一-鿿＀-￯]")


def main() -> int:
    files = [p for p in Path("crates").rglob("*") if p.is_file() and "target" not in p.parts]
    files.append(Path("README.md"))
    problems = []
    for path in files:
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        for no, line in enumerate(text.split("\n"), 1):
            if CJK.search(line):
                problems.append(f"{path}:{no}: 日本語の文字がある")
    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    print(f"check_english: 問題はありません（{len(files)} ファイル）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
