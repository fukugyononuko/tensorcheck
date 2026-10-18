"""check_privacy.py 自体のテスト（issue #43）。

一時の git リポジトリ（と、その git worktree）を作り、check_privacy.py を実行して終了コードを確かめる。
守る条件：手元（CI でない）での実行は、禁止語を実際に照合したときだけ成功する。

- 元の作業ツリーにある .privacy-denylist の語を、元の作業ツリーでも worktree（その下の階層）でも見つけること
- 禁止語がどこにもないとき、手元では失敗し、CI（CI=true）では警告して成功すること
- 禁止語があり、文章に含まれないときは成功し、照合した語数を表示すること
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "check_privacy.py"
WORD = "zzdenylisttestwordzz"
GIT = ["git", "-c", "user.name=test", "-c", "user.email=noreply@github.com", "-c", "commit.gpgsign=false"]


def run(cwd: Path, ci: bool = False) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if k not in ("PRIVACY_DENYLIST", "PR_TEXT", "CI")}
    if ci:
        env["CI"] = "true"
    return subprocess.run([sys.executable, str(SCRIPT)], cwd=cwd, env=env, capture_output=True, text=True)


def make_repo(root: Path, text: str) -> tuple[Path, Path]:
    """元の作業ツリーと worktree を作り、その 2 つのパスを返す。"""
    main = root / "main"
    main.mkdir(parents=True)
    subprocess.run(["git", "init", "-q", "-b", "main"], cwd=main, check=True)
    (main / ".gitignore").write_text(".privacy-denylist\n", encoding="utf-8")
    (main / "README.md").write_text(text, encoding="utf-8")
    (main / "docs").mkdir()
    (main / "docs" / "a.md").write_text("nothing here\n", encoding="utf-8")
    subprocess.run(["git", "add", "."], cwd=main, check=True)
    subprocess.run([*GIT, "commit", "-q", "-m", "init"], cwd=main, check=True)
    tree = root / "tree"
    subprocess.run(["git", "worktree", "add", "-q", "-b", "other", str(tree)], cwd=main, check=True)
    return main, tree


def main() -> int:
    failures = []

    def expect(name: str, result: subprocess.CompletedProcess, code_ok: bool, needle: str = "") -> None:
        output = result.stdout + result.stderr
        if (result.returncode == 0) != code_ok or needle not in output:
            failures.append(f"{name}：終了コード {result.returncode}、出力 {output.strip()!r}")

    with tempfile.TemporaryDirectory() as tmp:
        main_tree, worktree = make_repo(Path(tmp) / "hit", f"contains {WORD}\n")
        (main_tree / ".privacy-denylist").write_text(f"{WORD}\n", encoding="utf-8")
        expect("元の作業ツリーで禁止語を見つける", run(main_tree), False, "禁止語（1 行目の語）")
        expect("worktree で元の作業ツリーの禁止語を見つける", run(worktree), False, "禁止語（1 行目の語）")
        # 禁止語は直下の README.md にある。下の階層から実行しても、作業ツリー全体を照合すること
        expect("worktree の下の階層から実行しても見つける", run(worktree / "docs"), False, "禁止語（1 行目の語）")

        main_tree, worktree = make_repo(Path(tmp) / "clean", "nothing here\n")
        (main_tree / ".privacy-denylist").write_text(f"{WORD}\nanother\n", encoding="utf-8")
        expect("禁止語を含まなければ成功する（worktree）", run(worktree), True, "禁止語 2 語")

        main_tree, worktree = make_repo(Path(tmp) / "none", "nothing here\n")
        expect("禁止語がなければ手元では失敗する（元の作業ツリー）", run(main_tree), False)
        expect("禁止語がなければ手元では失敗する（worktree）", run(worktree), False)
        expect("禁止語がなければ CI では警告して成功する", run(main_tree, ci=True), True, "飛ばしました")

    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print("test_check_privacy: 7 つの場合がすべて期待どおりでした")
    return 0


if __name__ == "__main__":
    sys.exit(main())
