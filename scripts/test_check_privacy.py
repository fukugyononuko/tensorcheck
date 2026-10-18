"""check_privacy.py 自体のテスト（issue #43）。

一時の git リポジトリ（と、その git worktree）を作り、check_privacy.py を実行して終了コードを確かめる。
守る条件：手元（CI でない）での実行は、禁止語を実際に照合したときだけ成功する。

- 元の作業ツリーにある .privacy-denylist の語を、元の作業ツリーでも worktree（その下の階層）でも見つけること。
  .git を作業ツリーの外に置いた構成（--separate-git-dir）では、core.worktree（絶対パスか、共通の .git を基準にした相対パス）があれば同じ。なければ（共通の .git の親に関係のない .privacy-denylist があっても）手元では成功しないこと
- bare リポジトリの worktree では、bare リポジトリの親の .privacy-denylist を使わないこと
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


def make_repo(root: Path, text: str, separate: bool = False) -> tuple[Path, Path]:
    """元の作業ツリーと worktree を作り、その 2 つのパスを返す。

    separate なら、.git を作業ツリーの外に置く（git init --separate-git-dir）。共通の .git の親が元の作業ツリーでない構成。
    """
    main = root / "main"
    main.mkdir(parents=True)
    (root / "meta").mkdir()
    init = ["git", "init", "-q", "-b", "main"] + (["--separate-git-dir", str(root / "meta" / "repo.git")] if separate else [])
    subprocess.run(init, cwd=main, check=True)
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

        main_tree, worktree = make_repo(Path(tmp) / "separate", f"contains {WORD}\n", separate=True)
        (main_tree / ".privacy-denylist").write_text(f"{WORD}\n", encoding="utf-8")
        # core.worktree がないと git に元の作業ツリーの記録がない。見つけられなくても、手元では成功しないこと
        expect(".git が外で core.worktree がなければ、手元では成功しない", run(worktree), False, "禁止語が見つかりません")
        # 共通の .git の親（meta）は元の作業ツリーではない。そこに関係のない禁止語があっても使わないこと
        decoy = main_tree.parent / "meta" / ".privacy-denylist"
        decoy.write_text("unrelatedword\n", encoding="utf-8")
        expect(".git が外で core.worktree がなければ、共通の .git の親の禁止語を使わない", run(worktree), False, "禁止語が見つかりません")
        # core.worktree は絶対パスでも、共通の .git（meta/repo.git）を基準にした相対パスでもよい
        for kind, value in (("絶対パス", str(main_tree)), ("相対パス", "../../main")):
            subprocess.run(["git", "config", "core.worktree", value], cwd=main_tree, check=True)
            expect(f".git が外でも、core.worktree（{kind}）があれば worktree で見つける", run(worktree), False, "禁止語（1 行目の語）")
            expect(f".git が外でも、core.worktree（{kind}）があれば CI で検査を飛ばさない", run(worktree, ci=True), False, "禁止語（1 行目の語）")

        # bare リポジトリ（meta/.git）の worktree。meta は作業ツリーではないので、そこの関係のない禁止語を使わないこと
        source, _ = make_repo(Path(tmp) / "bare" / "source", f"contains {WORD}\n")
        bare = Path(tmp) / "bare" / "meta" / ".git"
        subprocess.run(["git", "clone", "-q", "--bare", str(source), str(bare)], check=True)
        bare_tree = Path(tmp) / "bare" / "tree"
        subprocess.run(["git", "--git-dir", str(bare), "worktree", "add", "-q", "-b", "other2", str(bare_tree), "main"], check=True)
        (bare.parent / ".privacy-denylist").write_text("unrelatedword\n", encoding="utf-8")
        expect("bare リポジトリの親の禁止語を使わない", run(bare_tree), False, "禁止語が見つかりません")

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
    print("test_check_privacy: 14 の場合がすべて期待どおりでした")
    return 0


if __name__ == "__main__":
    sys.exit(main())
