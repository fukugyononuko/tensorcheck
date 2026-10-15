"""重力の例題（grav01〜grav04）を実行し、その CSV から docs/figures/ の SVG を作る。

標準ライブラリだけで SVG を書く（依存を増やさないため）。図は例題の出力から
作るので、例題を直したらこのスクリプトを実行し直して、図もコミットする。

    python3 scripts/make_figures.py
"""

import csv
import math
import subprocess
import sys
import tempfile
from pathlib import Path

OUT = Path("docs/figures")
W, H = 640, 480
PAD_L, PAD_R, PAD_T, PAD_B = 72, 24, 68, 56

# 配色は dataviz のリファレンスの配色（カテゴリの 1〜3 番、青の順序の段階）。
# ライトとダークを、SVG の中の CSS で切り替える。CSS の変数は使わない
# （rsvg などの描画系が対応していないため）。色ごとに st-（線）と fi-（塗り）のクラスを作る。
LIGHT = {"surface": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e", "grid": "#e4e3df",
         "s1": "#2a78d6", "s2": "#eb6834", "s3": "#1baf7a",
         "q1": "#6da7ec", "q2": "#3987e5", "q3": "#2a78d6", "q4": "#256abf", "q5": "#184f95",
         "q6": "#0d366b"}
DARK = {"surface": "#1a1a19", "ink": "#f0efec", "ink2": "#c3c2b7", "grid": "#3a3a37",
        "s1": "#3987e5", "s2": "#d95926", "s3": "#199e70",
        "q1": "#184f95", "q2": "#256abf", "q3": "#2a78d6", "q4": "#3987e5", "q5": "#6da7ec",
        "q6": "#b7d3f6"}


def color_rules(pal: dict[str, str]) -> str:
    rules = [f".st-{k}{{stroke:{v}}} .fi-{k}{{fill:{v}}}" for k, v in pal.items()]
    rules.append(f"text{{fill:{pal['ink2']}}} .title{{fill:{pal['ink']}}}")
    rules.append(f".axis{{stroke:{pal['ink2']}}} .grid{{stroke:{pal['grid']}}}")
    return "\n".join(rules)


STYLE = f"""
<style>
  text {{ font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; font-size: 13px; }}
  .title {{ font-size: 16px; font-weight: 600; }}
  .axis {{ stroke-width: 1; fill: none; }}
  .grid {{ stroke-width: 1; }}
  .line {{ fill: none; stroke-width: 2; stroke-linejoin: round; stroke-linecap: round; }}
  {color_rules(LIGHT)}
  @media (prefers-color-scheme: dark) {{
  {color_rules(DARK)}
  }}
</style>
"""


def run_example(name: str, tmp: Path) -> list[dict]:
    path = tmp / f"{name}.csv"
    subprocess.run(
        ["cargo", "run", "--release", "-q", "-p", "tensorcheck-numerics",
         "--example", name, "--", "--csv", str(path)],
        check=True, stdout=subprocess.DEVNULL)
    with path.open() as f:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]


class Axes:
    """データ座標から SVG の座標への変換。log=True の軸は常用対数で写す。"""

    def __init__(self, xr, yr, xlog=False, ylog=False, equal=False, height=H):
        self.height = height
        self.xlog, self.ylog = xlog, ylog
        self.x0, self.x1 = (math.log10(v) if xlog else v for v in xr)
        self.y0, self.y1 = (math.log10(v) if ylog else v for v in yr)
        self.left, self.top = PAD_L, PAD_T
        self.w, self.h = W - PAD_L - PAD_R, height - PAD_T - PAD_B
        if equal:  # 1 単位の長さを x と y で同じにし、余った幅は中央に寄せる
            sx, sy = self.w / (self.x1 - self.x0), self.h / (self.y1 - self.y0)
            s = min(sx, sy)
            self.left += (self.w - s * (self.x1 - self.x0)) / 2
            self.top += (self.h - s * (self.y1 - self.y0)) / 2
            self.w, self.h = s * (self.x1 - self.x0), s * (self.y1 - self.y0)

    def px(self, x):
        x = math.log10(x) if self.xlog else x
        return self.left + (x - self.x0) / (self.x1 - self.x0) * self.w

    def py(self, y):
        y = math.log10(y) if self.ylog else y
        return self.top + (self.y1 - y) / (self.y1 - self.y0) * self.h

    def frame(self, xticks, yticks, xlabel, ylabel):
        out = [f'<clipPath id="plot"><rect x="{self.left:.1f}" y="{self.top:.1f}" '
               f'width="{self.w:.1f}" height="{self.h:.1f}"/></clipPath>']
        for v, label in xticks:
            x = self.px(v)
            out.append(f'<line class="grid" x1="{x:.1f}" x2="{x:.1f}" y1="{self.top:.1f}" '
                       f'y2="{self.top + self.h:.1f}"/>')
            out.append(f'<text x="{x:.1f}" y="{self.top + self.h + 18:.1f}" '
                       f'text-anchor="middle">{label}</text>')
        for v, label in yticks:
            y = self.py(v)
            out.append(f'<line class="grid" x1="{self.left:.1f}" x2="{self.left + self.w:.1f}" '
                       f'y1="{y:.1f}" y2="{y:.1f}"/>')
            out.append(f'<text x="{self.left - 8:.1f}" y="{y + 4:.1f}" '
                       f'text-anchor="end">{label}</text>')
        out.append(f'<rect class="axis" x="{self.left:.1f}" y="{self.top:.1f}" '
                   f'width="{self.w:.1f}" height="{self.h:.1f}"/>')
        out.append(f'<text x="{self.left + self.w / 2:.1f}" y="{self.height - 14}" '
                   f'text-anchor="middle">{xlabel}</text>')
        out.append(f'<text transform="translate(18 {self.top + self.h / 2:.1f}) rotate(-90)" '
                   f'text-anchor="middle">{ylabel}</text>')
        return out

    def polyline(self, xs, ys, color, extra=""):
        pts = " ".join(f"{self.px(x):.1f},{self.py(y):.1f}" for x, y in zip(xs, ys))
        return (f'<polyline class="line st-{color}" clip-path="url(#plot)" '
                f'{extra} points="{pts}"/>')

    def disk(self, x, y, r, cls, extra=""):
        rx = abs(self.px(x + r) - self.px(x))
        return (f'<circle cx="{self.px(x):.1f}" cy="{self.py(y):.1f}" r="{rx:.1f}" '
                f'class="{cls}" {extra}/>')

    def label(self, x, y, text, anchor="start", dx=6, dy=4):
        return (f'<text x="{self.px(x) + dx:.1f}" y="{self.py(y) + dy:.1f}" '
                f'text-anchor="{anchor}">{text}</text>')


def svg(title: str, subtitle: str, body: list[str], desc: str, h: int = H) -> str:
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {h}" '
            f'width="{W}" height="{h}" role="img" aria-label="{title}">\n'
            f"<title>{title}</title><desc>{desc}</desc>{STYLE}"
            f'<rect width="{W}" height="{h}" class="fi-surface"/>\n'
            f'<text class="title" x="{PAD_L}" y="28">{title}</text>\n'
            f'<text x="{PAD_L}" y="50">{subtitle}</text>\n'
            + "\n".join(body) + "\n</svg>\n")


def legend(items: list[tuple[str, str, str]]) -> list[str]:
    """右上の凡例。items は (色, 種類 "line" か "dot", 文)。"""
    out, x, y = [], W - PAD_R - 230, PAD_T + 18
    for color, kind, text in items:
        if kind == "line":
            out.append(f'<line x1="{x}" x2="{x + 22}" y1="{y - 4}" y2="{y - 4}" '
                       f'class="st-{color}" stroke-width="2"/>')
        else:
            out.append(f'<circle cx="{x + 11}" cy="{y - 4}" r="4.5" '
                       f'class="fi-{color} st-surface" stroke-width="2"/>')
        out.append(f'<text x="{x + 30}" y="{y}">{text}</text>')
        y += 20
    return out


def fig_curvature(rows):
    ax = Axes((1, 150), (1e-13, 30), xlog=True, ylog=True)
    body = ax.frame([(10 ** k, f"10{sup(k)}") for k in range(0, 3)],
                    [(10 ** k, f"10{sup(k)}") for k in range(-12, 2, 3)],
                    "r / r_s", "K r_s⁴")
    xs = [r["r_over_rs"] for r in rows]
    body.append(ax.polyline(xs, [r["k_exact_rs4"] for r in rows], "s1"))
    for r in rows[::2]:
        body.append(f'<circle cx="{ax.px(r["r_over_rs"]):.1f}" cy="{ax.py(r["k_rs4"]):.1f}" r="4.5" '
                    f'class="fi-s2 st-surface" stroke-width="2"/>')
    body += legend([("s1", "line", "12 r_s² / r⁶ (exact)"),
                    ("s2", "dot", "computed from g_ab")])
    return svg("Curvature of Schwarzschild spacetime",
               "R_ab = 0 everywhere (vacuum), yet K = R_abcd R^abcd ≠ 0", body,
               "Kretschmann scalar computed from the metric tensor by automatic "
               "differentiation, against the exact 12 r_s^2 / r^6.")


def fig_kepler(rows):
    ax = Axes((-37, 3), (-6, 6), equal=True, height=300)
    body = ax.frame([(v, str(v)) for v in range(-35, 1, 5)],
                    [(v, str(v)) for v in range(-5, 6, 5)], "x [AU]", "y [AU]")
    body.append(ax.polyline([r["x_au"] for r in rows], [r["y_au"] for r in rows], "s1"))
    n = 200
    body.append(ax.polyline([math.cos(2 * math.pi * k / n) for k in range(n + 1)],
                            [math.sin(2 * math.pi * k / n) for k in range(n + 1)], "s2"))
    body.append(ax.disk(0, 0, 0.25, "fi-ink"))
    body.append(ax.label(-17.4, -4.5, "Halley's comet", anchor="middle", dx=0, dy=-8))
    body.append(ax.label(-1.0, -1.0, "Earth's orbit", anchor="end", dx=-4, dy=14))
    body.append(ax.label(0, 0, "Sun", anchor="end", dx=-6, dy=-8))
    return svg("Halley's comet with Newton's gravity, in SI units",
               "one period (75.3 years) from Kepler's third law; returns to perihelion within 4e-8 of its distance", body,
               "Orbit of Halley's comet over one period, solved with Dormand-Prince "
               "from d2x/dt2 = -GM x / r^3, with Earth's orbit for scale.", h=300)


def fig_perihelion(rows):
    ax = Axes((-34, 34), (-34, 34), equal=True)
    body = ax.frame([(v, str(v)) for v in range(-30, 31, 10)],
                    [(v, str(v)) for v in range(-30, 31, 10)], "x / r_s", "y / r_s")
    body.append(ax.polyline([r["x_rs"] for r in rows], [r["y_rs"] for r in rows], "s1"))
    body.append(ax.disk(0, 0, 1, "fi-ink"))
    return svg("A precessing orbit near a compact object",
               "r_p = 10 r_s, r_a = 30 r_s: the perihelion advances 0.7439 rad per orbit", body,
               "Timelike geodesic of Schwarzschild spacetime over six orbits; the "
               "perihelion advances by the same angle each orbit.")


def fig_rays(rows):
    ax = Axes((-25, 25), (-18, 18), equal=True)
    body = ax.frame([(v, str(v)) for v in range(-20, 21, 10)],
                    [(v, str(v)) for v in range(-15, 16, 5)], "x / r_s", "y / r_s")
    rays = sorted({int(r["ray"]) for r in rows})
    for i in rays:
        pts = [r for r in rows if int(r["ray"]) == i]
        color = f"q{i + 1}"
        body.append(ax.polyline([p["x_rs"] for p in pts], [p["y_rs"] for p in pts], color))
        # 光線の出口側（枠の中の最後の点）にラベルを置く
        inside = [p for p in pts if -25 <= p["x_rs"] <= 25 and -18 <= p["y_rs"] <= 18]
        end = inside[-1]
        bottom = end["y_rs"] < -17
        body.append(ax.label(end["x_rs"], end["y_rs"], f"b = {pts[0]['b_rs']:g}",
                             anchor="end" if end["x_rs"] > 20 or bottom else "start",
                             dx=-6 if bottom else 6, dy=-6 if bottom else (-12 if end["x_rs"] > 20 else 4)))
    body.append(ax.disk(0, 0, 1, "fi-ink"))
    body.append(ax.disk(0, 0, 1.5, "st-ink2", 'fill="none" stroke-dasharray="3 3"'))
    return svg("Light rays bent by a black hole",
               "impact parameter b in r_s; dashed circle: photon sphere at 1.5 r_s", body,
               "Null geodesics of Schwarzschild spacetime entering from the left with "
               "impact parameters from 2.7 to 8 Schwarzschild radii.")


def sup(k: int) -> str:
    table = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return str(k).translate(table)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    figures = [
        ("grav01_schwarzschild_curvature", "schwarzschild_curvature.svg", fig_curvature),
        ("grav02_kepler", "halley.svg", fig_kepler),
        ("grav03_perihelion", "precession.svg", fig_perihelion),
        ("grav04_light_bending", "light_bending.svg", fig_rays),
    ]
    with tempfile.TemporaryDirectory() as tmp:
        for example, name, draw in figures:
            rows = run_example(example, Path(tmp))
            (OUT / name).write_text(draw(rows), encoding="utf-8")
            print(f"make_figures: {OUT / name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
