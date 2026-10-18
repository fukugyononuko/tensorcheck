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
    """git worktree の元の作業ツリー。わからなければ None。

    core.worktree があればそれ（相対パスなら git と同じく共通の .git を基準にする）。なければ、共通の .git が
    その親の .git そのものであるとき（通常の構成）だけ、その親を元の作業ツリーとする。.git を外に置き core.worktree も
    ない構成（git init --separate-git-dir の既定）では、git に元の作業ツリーの記録がないので None を返す。
    そこにある関係のない .privacy-denylist を使わず、手元では失敗させるため。
    """
    configured = subprocess.run(["git", "config", "--get", "core.worktree"], capture_output=True, text=True).stdout.strip()
    common = Path(git("rev-parse", "--path-format=absolute", "--git-common-dir").strip()).resolve()
    if configured:
        # 絶対パスなら、つないでもそのまま
        return (common / configured).resolve()
    if (common.parent / ".git").resolve() == common:
        return common.parent
    return None


def read_denylist() -> str:
    if os.environ.get("PRIVACY_DENYLIST"):
        return os.environ["PRIVACY_DENYLIST"]
    # git worktree の作業ツリーには管理外のファイルがないので、元の作業ツリーでも探す
    main_tree = main_worktree()
    for path in [Path(".privacy-denylist")] + ([main_tree / ".privacy-denylist"] if main_tree else []):
        if path.exists():
            return path.read_text(encoding="utf-8")
    return ""


def main() -> int:
    # 下の階層から実行しても、作業ツリー全体を照合する
    os.chdir(git("rev-parse", "--show-toplevel").strip())
    problems = []
    # On push, CI passes the blank between the empty title and body, so strip first.
    pr_text = os.environ.get("PR_TEXT", "").strip()
    identities = {e for e in git("log", "--format=%ae%n%ce").splitlines() if e}
    problems += [f"許可していないアドレスのコミットがあります：{e}" for e in sorted(identities - ALLOWED_EMAILS)]

    denylist = [w.strip().lower() for w in read_denylist().splitlines() if w.strip()]
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
