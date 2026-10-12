# 2a. 基本設計（外部設計）

外部設計は、インターフェースの設計である（[README](README.md)）。次の 2 つを一覧にする。

- **外部に提供するインターフェース**：公開する API（API）、利用者が実装する拡張点（EXT）、コンパイルエラーの文（ERR）
- **外部から利用するもの**（USE）

名前と書き方はここで決める。型の定義の細部と実装の方式は 03 で決める。コードの名前・ドキュメント・エラー文は英語で書く（NFR-07）。

## 2a.1 前提の確認（試作の結果）

決まり 4（フレームワークに依存する前提は、小さく試して確かめる）に従い、2a の前提を stable の Rust 1.99.0 で試した。

| 前提 | 結果 |
| --- | --- |
| ラベルが同じかを型の計算で判定し、縮約の相手を見つけて、結果の枠の並びを型として求められる | できた。$`g_{ab} V^b`$ は「$`a`$ が下付き」、$`T^a{}_c W^c`$ は「$`a`$ が上付き」の型になった。ラベルに型レベルの番号を持たせ、番号の一致を trait で判定する |
| `#[diagnostic::on_unimplemented]` で、誤りの見出しを決めた文にできる | できた。下付きどうしの縮約で、決めた見出し・ラベル・補足がエラーの先頭に出た |
| エラー文に、利用者のつけたラベルの名前が出る | 番号の型（`S<Z>` など）のままでは読めない。ラベルを名前つきの型（`b` など）にし、その型をエラー文に渡す（2a.4） |

制約：見出しは決めた文にできるが、その下に続く補足（どの trait が実装されていないか、など）には内部の型の名前が出る。NFR-05 は「見出しだけで誤りの種類がわかる」ことを求めるので、見出しに誤りの種類と、利用者の書いた名前（ラベル・座標系）を入れる。

## 2a.2 外部に提供するインターフェース：公開する API（API）

### クレートとモジュール

| ID | 場所 | 内容 | 構成要素 |
| --- | --- | --- | --- |
| API-01 | `tensorcheck::Tensor` | テンソルの型 `Tensor<D, S, Sym = NoSym, K = Tensorial, N = f64>`。`D` は次元、`S` は枠の並び、`Sym` は対称性、`K` はテンソル性（`Tensorial` か `ChartBound`）、`N` は数値の型 | M01 |
| API-02 | `tensorcheck::{Quantity, Vector, Covector}` | よく使う形の別名。`Quantity<D>` は階数 0、`Vector<D, C>` は `Up<C>` が 1 つ、`Covector<D, C>` は `Down<C>` が 1 つ | M01 |
| API-03 | `tensorcheck::{Up, Down}` | 枠の型。`Up<C>` は座標系 `C` の上付き、`Down<C>` は下付き。空間は `C` から決まる | E02 |
| API-04 | `tensorcheck::{NoSym, Symmetric, Antisymmetric}` | 対称性の型 | E04 |
| API-05 | `tensorcheck::units::si` | SI の単位系、次元の別名（`Length`、`Time`、`Velocity`、`Force` など）、単位（`m`、`km`、`s`、`kg`、`N`、`J` など） | E01 |
| API-06 | `tensorcheck::units::natural` | 自然単位系（$`\hbar = c = 1`$）。次元の別名（`Mass`、`Dimensionless` など）と単位（`eV`、`GeV`、`reduced_planck_mass` など） | E01 |
| API-07 | `tensorcheck::dim::{Prod, Quot, Pow}` | 次元から次元を作る型（例：`Quot<Length, Time>` は速度） | E01 |
| API-08 | `tensorcheck::labels` | 用意したラベル：`a`〜`h`、`i`〜`n`、`mu`、`nu`、`rho`、`sigma`、`alpha`、`beta`、`I`〜`N`。利用者がラベルを足すマクロ `label!` | E02 |
| API-09 | `tensorcheck::slot::{_0, _1, _2, _3}` | 上げ下げで、どの枠かを指定する印 | E02 |
| API-10 | `tensorcheck::Point<C, N = f64>` | 座標系 `C` の点。座標ごとに次元つきの値を持つ | E02 |
| API-11 | `tensorcheck::spaces` | 用意した空間・座標系・計量テンソル・座標変換：3 次元 Euclid 空間（デカルト、極（2 次元）、円柱、球）、Minkowski 時空（デカルト）、半径を指定する 2 次元球面、それぞれの計量テンソルと、デカルト座標との座標変換 | E02, E03 |
| API-12 | `tensorcheck::VectorSpace` | 数値解法の状態の契約（U05）。テンソル、点、それらの組（要素が 2〜6 個）が実装する | U05 |
| API-16 | `tensorcheck::geometry` | 接続・曲率・共変微分（API-50〜API-55） | M03 |
| API-17 | `tensorcheck::autodiff::{Dual, Dual2}` | 自動微分の数値の型（U02）。利用者が直接使うことは少なく、計量テンソルやテンソル場を数値の型によらない形で書けば、ライブラリが使う | U02 |
| API-13 | `tensorcheck_numerics::ode` | 常微分方程式：`solve_ivp`、`solve_second_order`、解法の型、許容値、解、失敗 | M02 |
| API-14 | `tensorcheck_numerics::quad` | 定積分：`integrate`、解法の型、結果、失敗 | M02 |
| API-15 | `tensorcheck_numerics::diff` | 数値微分：`derivative`、`jacobian`、方式の型 | M02 |

### テンソルの操作（M01）

| ID | 書き方 | 意味 | 対応する操作（2b） |
| --- | --- | --- | --- |
| API-20 | `Quantity::new(9.8, m / s.pow::<2>())`、`9.8 * (m / s.pow::<2>())` | 単位を指定して量を作る | `create` |
| API-21 | `Tensor::from_components([...], unit)` | 単位を指定してテンソルを作る。成分は行優先の並び | `create` |
| API-34 | `Vector::from_parts((1.0 * (m / s), 0.25 * (rad / s)))` | 成分ごとに次元が違うとき（座標基底の成分など）は、成分を量の組で渡す。各成分の次元は型で確かめる | `create` |
| API-22 | `Tensor::<_, _, Symmetric>::from_independent([...], unit)` | 対称・反対称のテンソルを、独立な成分だけから作る | `create` |
| API-23 | `q.value_in(km)`、`t.components_in(unit)` | 単位を指定して値を取り出す | `read` |
| API-24 | `t.at((Theta,))`、`t.at((R, Theta))` | 座標を型で指定して成分を 1 つ取り出す。戻り値は `Quantity` | `component` |
| API-25 | `a + b`、`a - b`、`k * t`、`t / k` | 足し算・引き算・スカラー倍 | `add` / `sub` / `mul` / `div` |
| API-26 | `a.outer(&b)` | 外積 | `outer` |
| API-27 | `t.idx((a, b))` | 枠にラベルをつける | `label` |
| API-28 | `g.idx((a, b)) * v.idx(b)` | ラベルつきの式の積。同じラベルが上下に 1 組あれば縮約する（Einstein の規約） | `contract` |
| API-29 | `t.idx((a, a))` | 1 つのテンソルの中の縮約（トレース） | `trace` |
| API-30 | `expr.into_tensor()` | ラベルつきの式をテンソルに戻す。残った枠は、式の中で現れた順に並ぶ | — |
| API-31 | `t.raise(slot::_0, &metric, &point)`、`t.lower(...)` | 計量テンソルで添字を上げ下げする | `raise` / `lower` |
| API-32 | `t.transform(&transform, &point)` | 座標系を変換する | `transform` |
| API-33 | `t.symmetrize()`、`t.antisymmetrize()` | 対称化・反対称化 | `symmetrize` / `antisymmetrize` |

### 幾何（M03）

| ID | 書き方 | 意味 | 対応する操作（2b） |
| --- | --- | --- | --- |
| API-50 | `Christoffel<C, N>` | Christoffel 記号の型。`Tensor<Dimensionless, (Up<C>, Down<C>, Down<C>), NoSym, ChartBound, N>` の別名 | — |
| API-51 | `geometry::christoffel(&metric, &p)` | 点 `p` での $`\Gamma^a{}_{bc}`$ | `christoffel` |
| API-52 | `geometry::{riemann, ricci, ricci_scalar, einstein}(&metric, &p)` | 点 `p` での $`R^a{}_{bcd}`$、$`R_{ab}`$、$`R`$、$`G_{ab}`$ | `riemann` ほか |
| API-53 | `geometry::curvature(&metric, &p)` | 同じ点の $`g`$、$`\partial g`$、$`\partial\partial g`$ を 1 回だけ求め、$`\Gamma`$ と曲率をまとめて返す | `riemann` ほか |
| API-54 | `geometry::covariant_derivative(&field, &metric, &p)` | 点 `p` での $`\nabla_c T^{a\cdots}{}_{b\cdots}`$。結果は下付きの枠が 1 つ増えたテンソル | `covariant_derivative` |
| API-55 | `field.at(&p)` | テンソル場を点で評価する（`f64` で） | — |

### 数値計算（M02）

| ID | 書き方 | 意味 |
| --- | --- | --- |
| API-40 | `ode::solve_ivp(f, t0..t1, y0, method, tol)` | 1 階の常微分方程式 $`dy/dt = f(t, y)`$ を解く。`f` の戻り値の型は `<Y as VectorSpace>::Per<T>`（$`y`$ ÷ $`[t]`$） |
| API-41 | `ode::solve_second_order(a, t0..t1, (x0, v0), method)` | 2 階の常微分方程式 $`d^2x/dt^2 = a(t, x)`$ を解く（leapfrog など） |
| API-42 | `ode::{Rk4 { h }, DormandPrince54, Rosenbrock { jacobian }, Leapfrog { h }}` | 解法。`Rosenbrock` の `jacobian` は `Jacobian::Numerical` か、利用者の関数 |
| API-43 | `ode::Tolerance { rtol, atol }` | 許容値。`rtol` は無次元の数、`atol` は状態と同じ型 |
| API-44 | `ode::Solution { t, y, stats }`、`ode::SolveError` | 解と、失敗（`MaxSteps`、`StepTooSmall`、`NonFinite`）。失敗には起きた時刻を含める |
| API-45 | `ode::solve_ivp_with(f, ..., on_step)` | 受け入れた各歩で呼ばれる関数を渡す（途中の値を記録するため） |
| API-46 | `quad::integrate(f, a..b, method, tol)` | 定積分 $`\int_a^b f(x)\,dx`$。結果の型は `<Y as VectorSpace>::Times<X>`（$`[f]\,[x]`$） |
| API-47 | `quad::{Trapezoid { n }, Simpson { n }, Romberg, GaussLegendre { n }, GaussKronrod}`、`quad::Tolerance { rtol, atol }` | 解法と許容値。`GaussKronrod` は適応型（7 点 Gauss、15 点 Kronrod）。`atol` は結果と同じ型 |
| API-48 | `quad::Integral { value, error_estimate, evaluations }`、`quad::QuadError` | 結果と、失敗（`ToleranceNotMet`、`NonFinite`） |
| API-49 | `diff::derivative(f, x, method)`、`diff::jacobian(f, y)` | 数値微分と Jacobi 行列。方式は `Central`、`Richardson` |

### 使い方の例

UC-01（F = m a）：

```rust
use tensorcheck::units::si::{kg, m, s, N};
use tensorcheck::spaces::euclid3::Cartesian;
use tensorcheck::{Quantity, Vector};

let mass = Quantity::new(2.0, kg);
let acc: Vector<_, Cartesian> = Vector::from_components([0.0, 0.0, -9.8], m / s.pow::<2>());
let force = mass * acc;                    // Vector<Force, Cartesian>
let fz = force.components_in(N)[2];        // -19.6
// let bad = force + acc;                  // compile error: dimension mismatch (ERR-01)
```

UC-02（極座標での速度の 2 乗）：

```rust
use tensorcheck::labels::{a, b};
use tensorcheck::spaces::polar::{Polar, PolarMetric, R, Theta};
use tensorcheck::units::si::{m, rad, s, Velocity};
use tensorcheck::{Point, Vector};

let p = Point::<Polar>::new((2.0 * m, 0.5 * rad));
let v: Vector<Velocity, Polar> = Vector::from_parts((1.0 * (m / s), 0.25 * (rad / s)));
let g = PolarMetric.at(&p);                               // g_ab, symmetric
let v2 = (g.idx((a, b)) * v.idx(a) * v.idx(b)).into_tensor(); // Quantity<Quot<Area, Pow<Time, 2>>>
let v_theta = v.at((Theta,));                             // Quantity with dimension 1/s
```

UC-03（場空間と時空）：

```rust
use tensorcheck::labels::{mu, nu, I, J};

// d_phi: Tensor<D, (Down<Minkowski>, Up<FieldChart>)>
let d_phi_up = d_phi.raise(slot::_0, &MinkowskiMetric, &x);
let kinetic = (field_metric.at(&phi).idx((I, J)) * d_phi.idx((mu, I)) * d_phi_up.idx((mu, J)))
    .into_tensor();
// d_phi.idx((mu, mu))  // compile error: index `mu` belongs to different spaces (ERR-07)
```

UC-07（硬い方程式）：

```rust
use tensorcheck_numerics::ode::{solve_ivp, Rosenbrock, Jacobian, Tolerance};

let sol = solve_ivp(rhs, t0..t1, y0, Rosenbrock { jacobian: Jacobian::Numerical },
                    Tolerance { rtol: 1e-6, atol: atol_like_y })?;
```

UC-08（定積分）：

```rust
use tensorcheck_numerics::quad::{integrate, GaussKronrod, Tolerance};

let work = integrate(|x| force_at(x), x0..x1, GaussKronrod, Tolerance { rtol: 1e-10, atol: 1e-12 * J })?;
// work.value has the dimension of energy
```

UC-09（球面の曲率と測地線）：

```rust
use tensorcheck::geometry::{christoffel, ricci_scalar};
use tensorcheck::labels::{a, b, c};
use tensorcheck::spaces::sphere2::{Sphere2, SphereMetric};

let g = SphereMetric::new(1.0 * m);                       // g = a^2 (dθ^2 + sin^2θ dφ^2)
let p = Point::<Sphere2>::new((0.7 * rad, 0.0 * rad));
let r = ricci_scalar(&g, &p);                             // 2 / a^2, dimension 1/m^2
let gamma = christoffel(&g, &p);                          // Christoffel<Sphere2>: not a tensor
// gamma.transform(&ToCartesian, &p);                     // compile error: not a tensor (ERR-16)

// geodesic: dx/dλ = u, du/dλ = -Γ^a_bc u^b u^c
let rhs = |_l, (x, u): &(Point<Sphere2>, Vector<_, Sphere2>)| {
    let gamma = christoffel(&g, x);
    (u.clone(), -(gamma.idx((a, b, c)) * u.idx(b) * u.idx(c)).into_tensor())
};
```

UC-09 の計量テンソルの定義（利用者が書く場合）：

```rust
impl Metric for SphereMetric {
    type Chart = Sphere2;
    type Dim = Area;
    // generic over the scalar type N, so the library can evaluate it with dual numbers
    fn components<N: Scalar>(&self, p: &Point<Sphere2, N>) -> [[N; 2]; 2] {
        let a2 = N::from(self.radius_squared());
        let s = p.coord::<Theta>().sin();
        [[a2, N::zero()], [N::zero(), a2 * s * s]]
    }
}
```

### 縮約の書き方（OPEN-03 の結論）

メソッド（`idx`）と演算子（`*`）で書く。マクロは初期には用意しない。

- 理由：IDE の補完と `rustfmt` がそのまま効く。エラーの位置がマクロの展開先ではなく、利用者の書いた式を指す。
- ラベルは、名前つきの型（`a`、`mu`、`I` など）として用意する。エラー文には、この名前が出る。
- 用意したラベルで足りないときは、利用者が `label!(kappa);` で足せる。ラベルの番号は自動で振らないので、`label!` は番号を明示する形にする（03 で決める）。

## 2a.3 外部に提供するインターフェース：拡張点（EXT）

利用者が実装する trait。定義は Engine の側に置く（01 §1.7）。

| ID | trait | 利用者が与えるもの | 定義する構成要素 | 変動性 |
| --- | --- | --- | --- | --- |
| EXT-01 | `UnitSystem` | 基本次元の数（型レベルの整数）と、次元の名前（ドキュメントとエラー文用） | E01 | VOL-01 |
| EXT-02 | `Unit` | 属する単位系、次元、基準の単位への換算係数 | E01 | VOL-02 |
| EXT-03 | `Space` | 次元数 $`n`$（型レベルの整数）、名前 | E02 | VOL-03 |
| EXT-04 | `Chart` | 属する空間、座標の型の並び、各座標の次元、基底の種類（`CoordinateBasis` か `NonCoordinateBasis`。`chart!` マクロの既定は座標基底）、名前 | E02 | VOL-04 |
| EXT-05 | `ChartTransform` | 変換元と変換先の座標系、点の写し方、Jacobi 行列 $`\partial x'^a / \partial x^b`$。逆は省略できる（省略すると U04 で求める） | E02 | VOL-04 |
| EXT-06 | `Metric` | 座標系、全体の次元、点での成分 $`g_{ab}(x)`$ を数値の型 `N` によらない形で。逆は省略できる | E03 | VOL-05 |
| EXT-08 | `TensorField` | 座標系、返すテンソルの型、点での成分を数値の型 `N` によらない形で | M03 | — |
| EXT-07 | 関数（`Fn(T, &Y) -> Y::Per<T>` など） | 右辺、加速度、被積分関数、微分する関数、Jacobi 行列 | M02 | — |

単位系・空間・座標系の定義は、構造体の定義と trait の実装を毎回書くと長くなるので、それを生成するマクロ（`unit_system!`、`space!`、`chart!`）も用意する。マクロが生成するのは EXT の実装だけで、ライブラリの内部には触れない。

## 2a.4 外部に提供するインターフェース：コンパイルエラーの文（ERR）

各誤りに専用の trait を用意し、その trait に `#[diagnostic::on_unimplemented]` で見出し（message）、位置の説明（label）、補足（note）をつける。`{…}` には利用者の書いた型の名前（ラベル、座標系、単位）が入る。

| ID | 誤り | 見出し（message） | 補足（note） | 判定する構成要素 |
| --- | --- | --- | --- | --- |
| ERR-01 | 次元の違う量の足し算・引き算 | `cannot add or subtract quantities with different dimensions` | `both operands must have the same physical dimension` | E01 |
| ERR-02 | 単位の次元が量の次元と違う | `` unit `{U}` does not have the dimension of this quantity `` | `use a unit of the same dimension, or convert the quantity first` | E01 |
| ERR-03 | 枠の並びの違うテンソルの足し算（階数・形・上下・座標系） | `cannot add tensors with different index structure` | `rank, index positions (upper/lower) and coordinate systems must all match` | E02 |
| ERR-04 | ラベルの数が階数と違う | `the number of index labels does not match the rank of the tensor` | — | E02 |
| ERR-05 | 同じラベルが上と上、または下と下 | `` cannot contract index `{L}`: it is upper in both places ``（下と下のときは `lower`） | `a repeated index must be upper in one place and lower in the other; raise or lower one with a metric first` | E02 |
| ERR-06 | 同じラベルが 3 回以上 | `` index `{L}` appears more than twice `` | `use a different label for each pair of contracted indices` | E02 |
| ERR-07 | 縮約するラベルの空間が違う | `` index `{L}` belongs to different spaces `` | `contracted indices must belong to the same space (e.g. both spacetime, or both field space)` | E02 |
| ERR-08 | 縮約するラベルの座標系が違う | `` index `{L}` is in chart `{C1}` in one place and `{C2}` in the other `` | `transform one tensor to the other chart first` | E02 |
| ERR-09 | 計量テンソルの座標系と枠の座標系が違う | `` metric `{G}` is defined on chart `{C1}`, but the index is in chart `{C2}` `` | — | E02 |
| ERR-10 | 座標変換の変換元と、テンソルの座標系が違う | `` transform `{X}` starts from chart `{C1}`, but the tensor is in chart `{C2}` `` | — | E02 |
| ERR-11 | 指定した座標が、その枠の座標系のものでない | `` coordinate `{K}` does not belong to chart `{C}` `` | — | E02 |
| ERR-12 | 枠が同じでないテンソルを対称・反対称にする | `only tensors whose indices are all of the same kind can be symmetric or antisymmetric` | `all indices must be in the same space and chart, and all upper or all lower` | E04 |
| ERR-13 | 対称性の違うテンソルの足し算 | `cannot add tensors with different symmetry` | `symmetrize or antisymmetrize one of them first` | E04 |
| ERR-14 | 右辺の型が「状態 ÷ 時間」でない | `the right-hand side must return the state type divided by the time dimension` | `for dy/dt = f(t, y), f must return a value whose dimension is [y]/[t]` | M02 |
| ERR-15 | 許容値の絶対誤差の型が状態と違う | `the absolute tolerance must have the same type as the state` | — | M02 |
| ERR-16 | テンソルでないもの（Christoffel 記号など）の座標変換 | `` `{T}` is not a tensor and cannot be transformed to another chart `` | `Christoffel symbols and partial derivatives depend on the chart; compute them in the target chart instead` | E02 |
| ERR-17 | テンソルとテンソルでないものの足し算 | `cannot add a tensor and a chart-bound (non-tensor) quantity` | `use geometry::covariant_derivative to combine partial derivatives and Christoffel symbols` | E02 |
| ERR-18 | 計量テンソル・テンソル場の成分が、数値の型によらない形で書かれていない | （差し替えられない。コンパイラ自身のエラー E0049「型引数の数が trait の宣言と合わない」が出る） | API ドキュメントと EXT-06・EXT-08 の説明に、数値の型によらない形で書く理由と書き方を載せる | コンパイラ |
| ERR-19 | 座標基底でない座標系で、接続・曲率・共変微分を求める | `` `{C}` is not a coordinate basis; Christoffel symbols and curvature need a coordinate basis `` | `connection coefficients in non-coordinate bases (e.g. orthonormal frames) are not supported yet` | M03 |

- 見出しの文言は実装時に、04 のコンパイル失敗のテスト（TC-F）で実際の出力を確かめて調整する。調整したらこの表も直す。
- 判定する構成要素が M02 の誤り（ERR-14、ERR-15）は、M02 の公開する関数の型の制約として表す。
- ERR-18 は、trait のメソッドの形の誤りなので、コンパイラ自身のエラーになり、見出しを差し替えられない。NFR-05 の例外として扱い、ドキュメントで補う。

## 2a.5 外部から利用するもの（USE）

| ID | 利用するもの | 用途 | 使う範囲 | 版の決め方 |
| --- | --- | --- | --- | --- |
| USE-01 | Rust の言語と標準ライブラリ（std） | すべて | 実行時 | `rust-toolchain.toml` で固定（1.99.0）。MSRV は OPEN-05 |
| USE-02 | コンパイラの診断機能 `#[diagnostic::on_unimplemented]` | エラー文（ERR） | コンパイル時 | USE-01 に含まれる（Rust 1.78 で安定化） |
| USE-03 | Cargo の workspace | 2 つのクレートの構成（01 §1.8） | 開発 | USE-01 に含まれる |
| USE-04 | コンパイル失敗のテスト用のクレート（trybuild） | ERR の文をテストする（NFR-08） | 開発だけ（`dev-dependencies`） | Dependabot で更新 |
| USE-05 | ベンチマーク用のクレート（criterion を候補とする） | NFR-02、NFR-04 の計測 | 開発だけ | 04 で決める |
| USE-06 | GitHub Actions | CI（fmt、clippy、test、doc、検査スクリプト） | 開発 | Action はコミット SHA で固定（NFR-14） |
| USE-07 | Python 3 | 検査スクリプト（個人情報、ワークフロー、呼び出しの規則） | 開発 | GitHub の実行環境のもの |
| USE-08 | Dependabot、secret scanning | 依存の更新、秘密情報の検査 | 開発 | — |

実行時に依存する外部のクレートはない（NFR-01）。`tensorcheck-numerics` は `tensorcheck` だけに依存する。
