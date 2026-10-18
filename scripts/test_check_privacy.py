"""check_privacy.py 自体のテスト（issue #43）。

一時の git リポジトリ（と、その git worktree）を作り、check_privacy.py を実行して終了コードを確かめる。
守る条件：手元（CI でない）での実行は、禁止語を実際に照合したときだけ成功する。

- 通常の構成で、元の作業ツリーにある .privacy-denylist の語を、元の作業ツリーでも worktree（その下の階層）でも見つけること
- 通常でない構成（.git を作業ツリーの外に置く --separate-git-dir、bare リポジトリ）では元の作業ツリーを探さず、
  共通の .git の親に関係のない .privacy-denylist があっても（その親が別のリポジトリの作業ツリーでも）、core.worktree があっても、手元では成功しないこと。
  PRIVACY_DENYLIST で与えれば照合すること
- core.worktree で、元の作業ツリーが .git の親の外側にある構成でも、その親の禁止語を使わないこと
- 空白だけの .privacy-denylist や PRIVACY_DENYLIST は、次の候補（元の作業ツリー）の禁止語を隠さないこと
- 優先順位の低い出どころは、上に有効な語があれば読まないこと
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


def run(cwd: Path, ci: bool = False, denylist: str = "") -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if k not in ("PRIVACY_DENYLIST", "PR_TEXT", "CI")}
    if ci:
        env["CI"] = "true"
    if denylist:
        env["PRIVACY_DENYLIST"] = denylist
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
        # 共通の .git の親が別のリポジトリの作業ツリーでも、元の作業ツリーとしては使わないこと
        subprocess.run(["git", "init", "-q"], cwd=decoy.parent, check=True)
        expect("共通の .git の親が別のリポジトリの作業ツリーでも、その禁止語を使わない", run(worktree), False, "禁止語が見つかりません")
        # core.worktree を設定しても、通常でない構成では元の作業ツリーを探さない（git の設定の解釈を真似ないため）
        subprocess.run(["git", "config", "core.worktree", str(main_tree)], cwd=main_tree, check=True)
        expect(".git が外なら、core.worktree があっても手元では成功しない", run(worktree), False, "禁止語が見つかりません")
        expect(".git が外でも、PRIVACY_DENYLIST で与えれば照合する", run(worktree, denylist=WORD), False, "禁止語（1 行目の語）")

        # bare リポジトリ（meta/.git）の worktree。meta は作業ツリーではないので、そこの関係のない禁止語を使わないこと
        source, _ = make_repo(Path(tmp) / "bare" / "source", f"contains {WORD}\n")
        bare = Path(tmp) / "bare" / "meta" / ".git"
        subprocess.run(["git", "clone", "-q", "--bare", str(source), str(bare)], check=True)
        bare_tree = Path(tmp) / "bare" / "tree"
        subprocess.run(["git", "--git-dir", str(bare), "worktree", "add", "-q", "-b", "other2", str(bare_tree), "main"], check=True)
        (bare.parent / ".privacy-denylist").write_text("unrelatedword\n", encoding="utf-8")
        expect("bare リポジトリの親の禁止語を使わない", run(bare_tree), False, "禁止語が見つかりません")

        # core.worktree で、通常の名前の .git を持つ main の外側（container）を作業ツリーにした構成。
        # main は作業ツリーの内側にあるだけで、直下ではない。main の関係のない禁止語を使わないこと
        main_tree, worktree = make_repo(Path(tmp) / "container" / "inner", f"contains {WORD}\n")
        subprocess.run(["git", "config", "core.worktree", str(main_tree.parent)], cwd=main_tree, check=True)
        (main_tree.parent / ".privacy-denylist").write_text(f"{WORD}\n", encoding="utf-8")
        (main_tree / ".privacy-denylist").write_text("unrelatedword\n", encoding="utf-8")
        expect("core.worktree の内側にある .git の親の禁止語を使わない", run(worktree), False, "禁止語が見つかりません")

        # 空白だけの .privacy-denylist や PRIVACY_DENYLIST は、元の作業ツリーの禁止語を隠さないこと
        main_tree, worktree = make_repo(Path(tmp) / "blank", f"contains {WORD}\n")
        (main_tree / ".privacy-denylist").write_text(f"{WORD}\n", encoding="utf-8")
        (worktree / ".privacy-denylist").write_text("  \n\n", encoding="utf-8")
        expect("空白だけの worktree のファイルがあっても、元の作業ツリーの禁止語を使う（CI）", run(worktree, ci=True), False, "禁止語（1 行目の語）")
        expect("空白だけの PRIVACY_DENYLIST があっても、元の作業ツリーの禁止語を使う", run(worktree, denylist=" \n"), False, "禁止語（1 行目の語）")

        main_tree, worktree = make_repo(Path(tmp) / "clean", "nothing here\n")
        (main_tree / ".privacy-denylist").write_text(f"{WORD}\nanother\n", encoding="utf-8")
        expect("禁止語を含まなければ成功する（worktree）", run(worktree), True, "禁止語 2 語")
        # 優先順位の低い出どころは、上に有効な語があれば読まない（読めないものが置いてあっても止まらない）
        (worktree / ".privacy-denylist").mkdir()
        expect("PRIVACY_DENYLIST があれば、読めない .privacy-denylist を読まない", run(worktree, denylist="unrelatedword"), True, "禁止語 1 語")

        main_tree, worktree = make_repo(Path(tmp) / "none", "nothing here\n")
        expect("禁止語がなければ手元では失敗する（元の作業ツリー）", run(main_tree), False)
        expect("禁止語がなければ手元では失敗する（worktree）", run(worktree), False)
        expect("禁止語がなければ CI では警告して成功する", run(main_tree, ci=True), True, "飛ばしました")

    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print("test_check_privacy: 17 の場合がすべて期待どおりでした")
    return 0


if __name__ == "__main__":
    sys.exit(main())
