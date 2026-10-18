"""個人情報が公開リポジトリに入っていないかを検査する（ホームページの TC-N06 と同じ決まり）。

1. コミットの作者・コミッターのアドレスが、許可したものだけであること
2. 公開される文章（追跡しているファイル、すべてのコミットメッセージ、PR のタイトルと本文）に禁止語がないこと
   禁止語は公開リポジトリに書かず、環境変数 PRIVACY_DENYLIST（GitHub の Secret）か、
   git 管理外のファイル .privacy-denylist（作業ツリーの直下、なければ git worktree の元の作業ツリーの直下。main_worktree を参照）から読む。
   どれもなければ、CI では警告して検査を飛ばし、手元では失敗する（禁止語を照合せずに成功としないため。issue #43）。
"""

import os
import re
import subprocess
import sys
from pathlib import Path

ALLOWED_EMAILS = {
    "340638804+fukugyononuko@users.noreply.github.com",
    "fukugyononuko@gmail.com",
    "noreply@github.com",
    "49699333+dependabot[bot]@users.noreply.github.com",
}
TEXT_FILE = re.compile(r"\.(md|py|rs|toml|ya?ml|txt|json|csv|cfg|ini|lock)$")


def git(*args: str) -> str:
    return subprocess.run(["git", *args], check=True, capture_output=True, text=True).stdout


def main_worktree() -> Path | None:
    """git worktree の元の作業ツリー。通常の構成でなければ None。

    通常の構成とは、共通の .git がその親の .git そのもので、その親を git が作業ツリーとして扱うこと。
    それ以外（.git を外に置く --separate-git-dir、bare リポジトリ、core.worktree を使う構成）では探さない。
    git の設定の解釈を真似ると、外れたときに関係のない .privacy-denylist を使うおそれがあるため（PR #45 のレビュー）。
    そのときは手元では失敗し、PRIVACY_DENYLIST で与えてもらう。
    """
    common = Path(git("rev-parse", "--path-format=absolute", "--git-common-dir").strip()).resolve()
    if (common.parent / ".git").resolve() != common:
        return None
    # その親が、作業ツリーの直下そのものであること（親で動く git は、上で確かめた共通の .git を使う）。
    # bare リポジトリ meta/.git の親 meta は作業ツリーでなく、core.worktree で外側を作業ツリーにした構成では、
    # 親は作業ツリーの内側にあるだけ
    top = subprocess.run(["git", "-C", str(common.parent), "rev-parse", "--show-toplevel"], capture_output=True, text=True)
    return common.parent if top.returncode == 0 and Path(top.stdout.strip()).resolve() == common.parent else None


def denylist_sources():
    """禁止語の出どころを優先順に 1 つずつ返す。後の出どころは、前に有効な語がなかったときだけ読む。"""
    yield os.environ.get("PRIVACY_DENYLIST", "")
    yield Path(".privacy-denylist").read_text(encoding="utf-8") if Path(".privacy-denylist").exists() else ""
    # git worktree の作業ツリーには管理外のファイルがないので、元の作業ツリーでも探す
    main_tree = main_worktree()
    if main_tree and (main_tree / ".privacy-denylist").exists():
        yield (main_tree / ".privacy-denylist").read_text(encoding="utf-8")


def read_denylist() -> list[str]:
    """禁止語の一覧。有効な語が 1 つもない出どころ（空白だけのファイルなど）は飛ばし、次の候補を見る。"""
    for raw in denylist_sources():
        if words := [w.strip().lower() for w in raw.splitlines() if w.strip()]:
            return words
    return []


def main() -> int:
    # 下の階層から実行しても、作業ツリー全体を照合する
    os.chdir(git("rev-parse", "--show-toplevel").strip())
    problems = []
    # On push, CI passes the blank between the empty title and body, so strip first.
    pr_text = os.environ.get("PR_TEXT", "").strip()
    identities = {e for e in git("log", "--format=%ae%n%ce").splitlines() if e}
    problems += [f"許可していないアドレスのコミットがあります：{e}" for e in sorted(identities - ALLOWED_EMAILS)]

    denylist = read_denylist()
    if not denylist and os.environ.get("CI") != "true":
        problems.append(
            f"禁止語が見つかりません。PRIVACY_DENYLIST か .privacy-denylist（作業ツリーの直下か、git worktree の元の作業ツリーの直下）で与えてください"
        )
    elif not denylist:
        print("check_privacy: 禁止語が設定されていないため、文章の検査を飛ばしました（PRIVACY_DENYLIST か .privacy-denylist）")
    else:
        sources = [(f, Path(f).read_text(encoding="utf-8", errors="ignore")) for f in git("ls-files").splitlines() if TEXT_FILE.search(f)]
        log = git("log", "--format=%H%x00%B%x00").split("\0")
        sources += [(f"コミット {log[i].strip()[:7]} のメッセージ", log[i + 1]) for i in range(0, len(log) - 1, 2)]
        if pr_text:
            sources.append(("PR のタイトルと本文", pr_text))
        for where, text in sources:
            lowered = text.lower()
            # 禁止語そのものはログに出さない（CI のログも公開されるため）
            problems += [f"{where}：禁止語（{i + 1} 行目の語）が含まれています" for i, w in enumerate(denylist) if w in lowered]

    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    pr = "あり" if pr_text else "なし"
    print(f"check_privacy: 問題はありません（コミットのアドレス {len(identities)} 種類、禁止語 {len(denylist)} 語、PR の本文 {pr}）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
