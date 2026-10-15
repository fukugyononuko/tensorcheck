"""呼び出しの規則（01 §1.7、DD-01）と unsafe の禁止を検査する（TC-A01、TC-A02）。

- manager/ は engine/・access/・utility/ だけを使う。ほかの manager/ のファイルを使わない。
- engine/ は access/・utility/ だけを使う。ほかの engine/ のファイルを使わない。
- access/ は utility/ だけを使う。utility/ は utility/ だけを使う。
- src/contracts.rs（層をまたいで共有する型）は、どの層からも使える。契約自身は access/・utility/ だけを使う。
- tensorcheck::internal::engine（E08）を使ってよいのは tensorcheck-numerics の src/manager/ だけ。tensorcheck::internal::utility（U04）は numerics のどの層からも使える。examples/ と tests/ は使わない。
  numerics から tensorcheck::engine を直接使うのも、internal::engine と同じ扱いにする。
- 両方のクレートの lib.rs に #![forbid(unsafe_code)] がある。

パスは、正規表現ではなく字句に分けてから名前を解決して調べる（Codex のレビュー 2026-10-14、指摘 5）。
次の書き方を、すべて同じパスとして扱う。

- まとめた use（`use crate::{manager::tensor};`）、別名（`as`）、`self`・`super`・`$crate`、先頭の `::`
- モジュールを use してから相対で呼ぶ形（`use crate::manager; manager::tensor::f()`）
- 式や型の中の完全なパス（`crate::manager::tensor::f()`）
- 子モジュール（`mod x;`）や glob（`use super::*;`）から入った名前
- 層の外のモジュールでの再公開（lib.rs の `pub use manager::tensor::Tensor;`、`crate::geometry::…` など）を通した参照

コメントと文字列の中は見ない。glob で入った先が層のモジュールで、名前を決められないときは、違反として扱う（安全側）。

引数にディレクトリを渡すと、そこを crates/ の代わりに検査する（検査自体のテスト用。scripts/test_check_layers.py）。
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
EXTERN = {"core", "std", "alloc", "tensorcheck", "tensorcheck_numerics"}

TOKEN = re.compile(
    r"""
    (?P<ws>\s+)
  | (?P<line>//[^\n]*)
  | (?P<block>/\*)
  | (?P<raw>(?:br|r)(?P<hashes>\#*)")
  | (?P<str>b?"(?:\\.|[^"\\])*")
  | (?P<char>b?'(?:\\(?:x[0-9a-fA-F]{2}|u\{[0-9a-fA-F]+\}|.)|[^'\\\n])')
  | (?P<life>'[A-Za-z_]\w*)
  | (?P<ident>\$?(?:r\#)?[A-Za-z_]\w*)
  | (?P<path>::)
  | (?P<num>\d\w*(?:\.\d\w*)?)
  | (?P<punct>.)
    """,
    re.VERBOSE | re.DOTALL,
)


def tokenize(text: str) -> list[str]:
    """識別子・`::`・記号の列にする。コメント・文字列・文字・数は捨てる。"""
    out = []
    i = 0
    while i < len(text):
        m = TOKEN.match(text, i)
        assert m is not None
        kind = m.lastgroup if m.lastgroup != "hashes" else "raw"
        if kind == "block":
            depth, j = 1, m.end()
            while depth and j < len(text):
                if text.startswith("/*", j):
                    depth, j = depth + 1, j + 2
                elif text.startswith("*/", j):
                    depth, j = depth - 1, j + 2
                else:
                    j += 1
            i = j
            continue
        if m.group("raw"):
            end = text.find('"' + m.group("hashes"), m.end())
            i = len(text) if end < 0 else end + 1 + len(m.group("hashes"))
            continue
        if kind in ("ident", "path", "punct"):
            tok = m.group(0)
            out.append(tok[2:] if tok.startswith("r#") else tok)
        i = m.end()
    return out


class Module:
    """モジュール 1 つ（ファイル、またはファイルの中の `mod x { … }`）。"""

    def __init__(self, path: tuple[str, ...], file: Path, layer_file: tuple[str, str] | None):
        self.path = path  # crate からの絶対パス（("crate", "engine", "metric") など）
        self.file = file
        self.layer_file = layer_file  # このモジュールを含むファイルの（層, ファイル名）
        self.children: set[str] = set()  # `mod x;` / `mod x { … }`
        self.aliases: dict[str, tuple[str, ...]] = {}  # use で入れた名前 → 解決前のパス
        self.globs: list[tuple[str, ...]] = []  # use x::*; の x（解決前）
        self.pub_aliases: dict[str, tuple[str, ...]] = {}  # pub use（再公開）
        self.refs: list[tuple[str, ...]] = []  # 使っているパス（解決前）


def parse_use_tree(toks: list[str], i: int, prefix: tuple[str, ...]):
    """use の木を読み、(名前, パス, glob か) の列と、次の位置を返す。"""
    items = []
    path = list(prefix)
    if i < len(toks) and toks[i] == "::":
        i += 1  # 先頭の `::`（外部クレート）
    while i < len(toks):
        t = toks[i]
        if t == "{":
            i += 1
            while toks[i] != "}":
                sub, i = parse_use_tree(toks, i, tuple(path))
                items += sub
                if toks[i] == ",":
                    i += 1
            return items, i + 1
        if t == "*":
            items.append((None, tuple(path), True))
            return items, i + 1
        path.append(t)
        i += 1
        if i < len(toks) and toks[i] == "::":
            i += 1
            continue
        name = path[-1]
        if i + 1 < len(toks) and toks[i] == "as":
            name = toks[i + 1]
            i += 2
        if name == "self":  # use a::b::{self} は b を入れる
            path.pop()
            name = path[-1]
        items.append((name, tuple(path), False))
        return items, i
    return items, i


def module_path_of(rel: Path) -> tuple[str, ...]:
    parts = list(rel.with_suffix("").parts[1:])  # src/ を除く
    if parts and parts[-1] in ("mod", "lib", "main"):
        parts.pop()
    return ("crate", *parts)


def scan_file(path: Path, rel: Path, modules: dict[tuple[str, ...], Module]) -> None:
    parts = rel.parts
    layer = parts[1] if len(parts) > 2 and parts[0] == "src" and parts[1] in LAYERS else None
    layer_file = (layer, Path(parts[-1]).stem) if layer else None
    toks = tokenize(path.read_text(encoding="utf-8"))
    top = Module(module_path_of(rel), path, layer_file)
    modules[top.path] = top
    stack = [(top, 0)]  # (モジュール, そのモジュールが始まった波括弧の深さ)
    depth = 0
    i = 0
    while i < len(toks):
        mod = stack[-1][0]
        t = toks[i]
        if t == "{":
            depth += 1
        elif t == "}":
            depth -= 1
            if len(stack) > 1 and depth < stack[-1][1]:
                stack.pop()
        elif t == "mod" and i + 1 < len(toks) and re.match(r"\w+$", toks[i + 1]):
            name = toks[i + 1]
            mod.children.add(name)
            if i + 2 < len(toks) and toks[i + 2] == "{":
                inner = Module((*mod.path, name), path, layer_file)
                modules[inner.path] = inner
                depth += 1
                stack.append((inner, depth))
                i += 3
                continue
        elif t == "use":
            is_pub = i > 0 and toks[i - 1] in ("pub", ")")
            items, i = parse_use_tree(toks, i + 1, ())
            for name, p, glob in items:
                mod.refs.append(p)
                if glob:
                    mod.globs.append(p)
                else:
                    mod.aliases[name] = p
                    if is_pub:
                        mod.pub_aliases[name] = p
            continue
        elif re.match(r"\$?\w+$", t) and (i == 0 or toks[i - 1] != "::") and i + 1 < len(toks) and toks[i + 1] == "::":
            # 式・型の中のパス a::b::c（`::<` の総称の引数で切る）
            p = [t]
            j = i + 1
            while j + 1 < len(toks) and toks[j] == "::" and re.match(r"\w+$", toks[j + 1]):
                p.append(toks[j + 1])
                j += 2
            mod.refs.append(tuple(p))
            i = j
            continue
        i += 1


class Resolver:
    def __init__(self, modules: dict[tuple[str, ...], Module]):
        self.modules = modules

    def absolute(
        self, mod: Module, p: tuple[str, ...], seen=frozenset(), globs=True
    ) -> tuple[str, ...] | None:
        """mod の中で書かれたパス p を、crate からの絶対パス（外部なら先頭がクレート名）にする。"""
        if not p:
            return None
        head, rest = p[0], p[1:]
        if head in ("crate", "$crate"):
            return ("crate", *rest)
        if head in ("self", "super"):
            base = mod.path
            while p and p[0] in ("self", "super"):
                if p[0] == "super":
                    base = base[:-1] if len(base) > 1 else base
                p = p[1:]
            return (*base, *p)
        if head in mod.children:
            return (*mod.path, *p)
        if head in mod.aliases and (mod.path, head) not in seen:
            target = self.absolute(mod, mod.aliases[head], seen | {(mod.path, head)})
            return None if target is None else (*target, *rest)
        for g in mod.globs if globs else []:
            src = self.absolute(mod, g, seen, globs=False)
            if src is None:
                continue
            m = self.modules.get(src)
            if m is not None and (head in m.children or head in m.aliases):
                return (*src, *p)
        if head in EXTERN:
            return p
        return None  # ローカルの変数・総称の型・プレリュードなど

    def canonical(self, p: tuple[str, ...]) -> tuple[str, ...]:
        """層の外のモジュールでの再公開をたどる（crate::Tensor → crate::manager::tensor::Tensor）。"""
        for _ in range(16):
            if len(p) < 2 or p[0] != "crate" or (len(p) > 1 and p[1] in LAYERS):
                return p
            for k in range(len(p) - 1, 0, -1):
                m = self.modules.get(p[:k])
                if m is not None and p[k] in m.pub_aliases:
                    target = self.absolute(m, m.pub_aliases[p[k]])
                    if target is not None:
                        p = (*target, *p[k + 1 :])
                        break
            else:
                return p
        return p


def check(root: Path) -> list[str]:
    problems = []
    for crate in sorted(p for p in root.iterdir() if p.is_dir()):
        src = crate / "src"
        lib = src / "lib.rs"
        if lib.exists() and "#![forbid(unsafe_code)]" not in lib.read_text(encoding="utf-8"):
            problems.append(f"{lib}: #![forbid(unsafe_code)] がない")
        modules: dict[tuple[str, ...], Module] = {}
        for path in sorted(src.rglob("*.rs")) if src.exists() else []:
            scan_file(path, path.relative_to(crate), modules)
        resolver = Resolver(modules)
        for mod in modules.values():
            for p in mod.refs:
                problems += check_ref(crate, mod, p, resolver)
        # src/ の外（tests/・examples/・benches/）は internal だけを見る
        for path in sorted(crate.rglob("*.rs")):
            rel = path.relative_to(crate)
            if rel.parts[0] in ("src", "target"):
                continue
            outside: dict[tuple[str, ...], Module] = {}
            scan_file(path, Path("src") / rel.name, outside)
            r = Resolver(outside)
            for m in outside.values():
                for p in m.refs:
                    a = r.absolute(m, p)
                    if a and a[0] == "tensorcheck" and len(a) > 1 and a[1] == "internal":
                        problems.append(f"{path}: tensorcheck::internal は numerics のソースだけが使える")
    return sorted(set(problems))


def check_ref(crate: Path, mod: Module, p: tuple[str, ...], resolver: Resolver) -> list[str]:
    a = resolver.absolute(mod, p)
    if a is None:
        return []
    where = f"{mod.file}"
    shown = "::".join(p)
    layer = mod.layer_file[0] if mod.layer_file else None
    # numerics から core クレートへ
    if a[0] == "tensorcheck" and len(a) > 1:
        if crate.name != "tensorcheck-numerics":
            return [f"{where}: tensorcheck::internal は numerics のソースだけが使える"] if a[1] == "internal" else []
        kind = None
        if a[1] == "internal":
            if len(a) == 2:
                return []  # モジュールを入れただけ。使う側のパスで調べる
            kind = a[2]
        elif a[1] == "engine":
            kind = "engine"
        if kind is None:
            return []
        if kind == "engine" and layer != "manager":
            return [f"{where}: tensorcheck::internal::engine は numerics の manager だけが使える（{shown}）"]
        if kind not in ("engine", "utility"):
            return [f"{where}: tensorcheck::internal::{kind} は存在しない区分（{shown}）"]
        return []
    if a[0] != "crate":
        return []
    a = resolver.canonical(a)
    if len(a) < 2:
        return []
    is_contract = mod.path[:2] == ("crate", "contracts")
    if a[1] == "internal" and len(a) > 2 and crate.name == "tensorcheck":
        return []  # core クレートの中で internal を組み立てる側
    if a[1] not in LAYERS:
        return []
    target = a[1]
    target_module = a[2] if len(a) > 2 else None
    if is_contract:
        if target not in ("access", "utility"):
            return [f"{where}: 契約が {target} を使っている（{shown}）"]
        return []
    if layer is None:
        return []
    own = mod.layer_file[1]
    if target == layer and layer in ("manager", "engine"):
        if target_module is None:
            # 層のモジュールそのもの。名前を入れただけなら、使う側のパスで調べる。
            # glob（use crate::engine::*）は中身を決められないので違反とする
            return [f"{where}: {layer} どうしの呼び出し（{shown}::*）"] if p in mod.globs else []
        if target_module != own:
            return [f"{where}: {layer} どうしの呼び出し（{shown}）"]
        return []
    if target not in ALLOWED[layer]:
        return [f"{where}: {layer} から {target} を使っている（{shown}）"]
    return []


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
