# 3. 詳細設計

2b の契約と 2a のインターフェースを、Rust の型・trait・アルゴリズム・ファイル配置に落とす。
コードは設計を示すための抜粋で、細部（ライフタイム、`where` 節の全体など）は実装で決める。

## 3.0 詳細設計の項目一覧

| ID | 項目 | 構成要素 |
| --- | --- | --- |
| DD-01 | ファイル配置と呼び出しの規則 | すべて |
| DD-02 | 型レベルの整数・真偽・並び | U01 |
| DD-03 | 次元と単位 | E01 |
| DD-04 | 空間・座標系・点と、成分の持ち方 | E02, A01 |
| DD-05 | テンソルの型と基本の演算 | M01 |
| DD-06 | ラベルと縮約 | E02, M01 |
| DD-07 | 対称性 | E04, A01 |
| DD-08 | 計量テンソル・座標変換・逆行列 | E03, E02, U04 |
| DD-09 | 数値の型と自動微分 | U02, E08 |
| DD-10 | 接続と曲率、共変微分 | E09, M03 |
| DD-11 | 数値解法の状態 | U05 |
| DD-12 | 常微分方程式の解法 | E05, M02 |
| DD-13 | 誤差の測り方と刻み幅の調整 | E06 |
| DD-14 | 定積分 | E07, M02 |
| DD-15 | 差分による数値微分 | E08 |
| DD-16 | エラー文の実装 | U03 |
| DD-17 | 失敗の型 | M02, U04 |
| DD-18 | 利用者向けのマクロ | C01 向け |
| DD-19 | 試作で確かめたこと | — |
| DD-20 | ビルドと CI | — |

## 3.1 DD-01 ファイル配置と呼び出しの規則

```text
Cargo.toml                      workspace（members: crates/*）
crates/
  tensorcheck/
    src/
      lib.rs                    公開する名前の再輸出（API-01〜API-17）
      manager/   tensor.rs geometry.rs           M01, M03
      engine/    dimension.rs index.rs metric.rs symmetry.rs
                 differentiation.rs connection.rs E01〜E04, E08, E09
      access/    component.rs                    A01
      utility/   typelevel.rs scalar.rs autodiff.rs diagnostics.rs
                 linear.rs vector_space.rs       U01〜U05
      units/     si.rs natural.rs                API-05, API-06（EXT-01・EXT-02 の実装）
      spaces/    euclid3.rs polar.rs minkowski.rs sphere2.rs  API-11
      labels.rs  macros.rs
      internal.rs                                tensorcheck-numerics 向け（#[doc(hidden)]）
    tests/        ui/（コンパイル失敗のテスト）, *.rs
    examples/
  tensorcheck-numerics/
    src/
      lib.rs  ode.rs quad.rs diff.rs             API-13〜API-15
      manager/   solver.rs                       M02
      engine/    ode_step.rs error_control.rs quadrature.rs  E05〜E07
    tests/  examples/  benches/
scripts/  check_privacy.py check_workflows.py check_layers.py
```

- 層ごとにディレクトリを分ける。層の外に見せる項目は `pub(crate)`、クレートの外に見せる項目だけを `pub` にする。
- `scripts/check_layers.py`（NFR-11）は、各ファイルの `use crate::…` と `tensorcheck::internal::…` を読み、01 §1.7 の表に反しないかを検査する。
  - `manager/` は `engine/`・`access/`・`utility/` だけを使う。ほかの `manager/` を使わない。
  - `engine/` は `access/`・`utility/` だけを使う。ほかの `engine/` のファイルを使わない。拡張点の trait（EXT）は各 Engine のファイルに置く。
  - `access/` は `utility/` だけを使う。`utility/` は `utility/` だけを使う。
  - `internal` を使ってよいのは `tensorcheck-numerics/src/manager/` だけ。`examples/` と `tests/` は使わない。
- `#![forbid(unsafe_code)]` を両方のクレートの `lib.rs` に置く（NFR-12）。

## 3.2 DD-02 型レベルの整数・真偽・並び（U01）

| 種類 | 表現 | 演算 |
| --- | --- | --- |
| 自然数 | `Z`、`S<N>`（Peano 数） | 等しいか（`Eq`）、`Add` |
| 整数 | `Zero`、`Pos<S<N>>`、`Neg<S<N>>`（0 の表し方を 1 つにする） | `Add`、`Sub`、`Neg`、`MulNat`（自然数倍）、等しいか |
| 真偽 | `True`、`False` | `If<B, X, Y>` |
| 並び | `Nil`、`Cons<H, T>` | `Len`、`Concat`、`Zip`、`Map`（trait で表した関数を各要素に適用）、`RemoveMatch`（DD-06） |

- 型の計算は trait の関連型で表す（例：`<A as Add<B>>::Out`）。
- 次元の指数とラベルの番号は小さい数なので、Peano 数の再帰の深さは、コンパイラの既定の上限（128）に収まる。指数は $`\pm 20`$ まで、ラベルの番号は 63 までを保証の範囲とし、テストで確かめる。
- 空間の次元数 $`n`$ は型レベルの数にしない。空間ごとの関連定数 `N` と配列の型で持つ（DD-04）。

## 3.3 DD-03 次元と単位（E01）

```rust
pub trait UnitSystem: 'static {          // EXT-01
    type Count;                          // 基本次元の数（型レベルの自然数）
    const NAMES: &'static [&'static str];
}
pub struct Dim<Sys, Exps>(PhantomData<(Sys, Exps)>);   // Exps は整数の並び（長さ = Count）

pub trait Unit: Copy + 'static {         // EXT-02
    type Sys: UnitSystem;
    type Dim;                            // Dim<Self::Sys, …>
    const FACTOR: f64;                   // 1 この単位 = FACTOR × 基準の単位
    const SYMBOL: &'static str;
}
```

- 次元の積・商は、指数の並びを要素ごとに足す・引く（`Zip` と `Add`）。単位系が違えば、`SameSystem` を満たさずコンパイルエラーになる。
- 次元の一致は `SameDimension<D2>` で確かめる。`impl<D> SameDimension<D> for D` だけを実装し、この trait に ERR-01 の文をつける。
- 単位どうしの積・商・べき（`m / s`、`s.pow::<2>()`）は、`Prod<U1, U2>` などの型を作る。`FACTOR` は関連定数の浮動小数点の演算（`const` での演算は stable）で求める。
- 単位は比例するものだけを扱う（摂氏温度のような原点のずれる単位は扱わない）。
- 自然単位系（API-06）は、基本次元を「質量」の 1 つにする。$`\mathrm{eV}`$、$`\mathrm{GeV}`$、換算プランク質量 $`M_\mathrm{P} = 2.435 \times 10^{18}\,\mathrm{GeV}`$ を単位として用意する。基準の単位は $`\mathrm{GeV}`$。

## 3.4 DD-04 空間・座標系・点と、成分の持ち方（E02、A01）

空間の次元数を配列の長さに使うために、空間ごとに「長さ $`n`$ の配列の型」を関連型（GAT）として持たせる。stable の Rust では、型の計算で求めた数を配列の長さにできない（`generic_const_exprs` が要る）が、この方法なら長さの計算が要らない（DD-19）。

```rust
pub trait Space: 'static {               // EXT-03
    const N: usize;
    const NAME: &'static str;
    type Arr<T: Copy>: Copy + AsRef<[T]> + AsMut<[T]>;      // [T; N]
    type SymArr<T: Copy>: Copy + AsRef<[T]> + AsMut<[T]>;   // [T; N(N+1)/2]
    type AntiArr<T: Copy>: Copy + AsRef<[T]> + AsMut<[T]>;  // [T; N(N-1)/2]
    fn from_fn<T: Copy>(f: impl FnMut(usize) -> T) -> Self::Arr<T>;
}
pub trait Chart: 'static {               // EXT-04
    type Space: Space;
    type Coords;                         // 座標の印の型の並び（例：Cons<R, Cons<Theta, Nil>>）
    type CoordDims;                      // 各座標の次元の並び
    const NAME: &'static str;
}
pub struct Up<C>(PhantomData<C>);
pub struct Down<C>(PhantomData<C>);
pub struct Point<C: Chart, N: Scalar = f64> { x: <C::Space as Space>::Arr<N> }
```

- 具体的な配列の長さ（`[T; 4]` など）は、利用者が `space!` マクロ（DD-18）で空間を定義するときに書かれる。ライブラリの中では長さを計算しない。
- 枠の並び `S` から、成分の入れ物の型を入れ子の配列として組み立てる（A01）。

```rust
pub trait Slots {
    type Store<T: Copy>: Copy;           // Nil → T、Cons<枠, 残り> → <枠の空間>::Arr<残り::Store<T>>
    const LEN: usize;                    // 成分の数 = Π n_i（関連定数の演算。配列の長さには使わない）
    fn get<T: Copy>(s: &Self::Store<T>, flat: usize) -> T;
    fn set<T: Copy>(s: &mut Self::Store<T>, flat: usize, v: T);
    fn fill<T: Copy>(v: T) -> Self::Store<T>;
    fn extent(slot: usize) -> usize;     // slot 番目の枠の n
}
```

- 例：時空（4 次元）・場空間（2 次元）の枠を持つテンソルの入れ物は `[[f64; 2]; 4]`。Riemann テンソル（4 次元）は `[[[[f64; 4]; 4]; 4]; 4]`。
- 平らな番号 `flat` は行優先（最後の枠がいちばん速く回る）。`get` と `set` は `flat / 残りの LEN` と `flat % 残りの LEN` で入れ子をたどる。`unsafe` は使わない。
- 「階数と形の一致」は、枠の並びの型が等しいこと（空間・座標系・上下がすべて同じ）で確かめる。形（次元数）は空間で決まるので、別に比べない。

### 成分の次元（`component_dim`）

`t.at((Theta,))` は、座標の印 `Theta` が枠の座標系の `Coords` の何番目かを型の計算で求め（ないなら ERR-11）、その番号を平らな番号に、その座標の次元を戻り値の次元に使う。

```math
[T^{a_1 \cdots}{}_{b_1 \cdots}] = [T] \prod_i [x^{a_i}] \Big/ \prod_j [x^{b_j}]
```

## 3.5 DD-05 テンソルの型と基本の演算（M01）

```rust
pub struct Tensor<D, S: Slots, Sym = NoSym, K = Tensorial, N: Scalar = f64> {
    c: S::Store<N>,
    _m: PhantomData<(D, Sym, K)>,
}
pub type Quantity<D, N = f64> = Tensor<D, Nil, NoSym, Tensorial, N>;
pub type Vector<D, C, N = f64> = Tensor<D, Cons<Up<C>, Nil>, NoSym, Tensorial, N>;
pub struct Tensorial;
pub struct ChartBound;
```

| 演算 | 実装 | 型の制約（満たさないときのエラー） |
| --- | --- | --- |
| `a + b`、`a - b` | `impl Add for Tensor<…>` | `SameDimension`（ERR-01）、`SameSlots`（ERR-03）、`SameSymmetry`（ERR-13）、`SameKind`（ERR-17） |
| `k * t`、`t / k`（`k` は `Quantity` か `f64`） | `impl Mul`、`impl Div` | 次元は `Prod`・`Quot` |
| `q1 * q2`（スカラーどうし） | `impl Mul for Quantity` | — |
| `t.outer(&u)` | 入れ物を作り、$`a_i b_j`$ を置く | 枠の並びは `Concat`、テンソル性は両方が `Tensorial` なら `Tensorial` |
| `Quantity::new(x, unit)` | `x * U::FACTOR` を置く | `UnitMatches<D>`（ERR-02） |
| `q.value_in(unit)` | `c / U::FACTOR` | `UnitMatches<D>`（ERR-02） |
| `Tensor::from_components(arr, unit)` | 平らな配列（長さ `S::LEN`）を受け取る | 長さが違えば実行時ではなく型で止めたいが、配列の長さに `S::LEN` を使えないため、入れ子の配列 `S::Store<f64>` を受け取る形にする |
| `Vector::from_parts((q1, q2, …))` | 成分ごとに次元の違う量の組を受け取る | 各要素の次元が、DD-04 の式で求めた成分の次元と一致するか |

- 成分は基準の単位の数値で持つ。単位の換算は作るときと取り出すときだけに行う。
- `Copy` な型の入れ物なので、演算は値を返す（ヒープを使わない）。

## 3.6 DD-06 ラベルと縮約（E02、M01）

### ラベル

```rust
pub trait Label: 'static { type Id; const NAME: &'static str; }   // Id は型レベルの自然数
#[allow(non_camel_case_types)] pub struct mu;
impl Label for mu { type Id = /* 12 */; const NAME: &'static str = "mu"; }
```

- 用意するラベル（API-08）に、0〜39 の番号を振る。利用者は `label!(kappa = 40);` のように、40〜63 の番号を明示して足す。
- 2 つのラベルに同じ番号を振ると、別のラベルが同じものとして縮約されてしまう。これを型で止める方法は stable の Rust にないので、用意したラベルについては 04 のテストで重なりがないことを確かめ、利用者が足すラベルについては `label!` のドキュメントで注意する。
- エラー文に `{L}` として出るのは、ラベルの型の名前（`mu` など）。

### 縮約の型の計算

`t.idx((a, b))` は、枠の並びとラベルの並びを組にした「ラベルつきの枠の並び」 `Cons<Ix<a, Down<C>>, Cons<Ix<b, Up<C>>, Nil>>` を持つ `Labeled<'_, T, …>` を返す（ラベルの数が階数と違えば ERR-04）。

`Labeled<A> * Labeled<B>` は、次の型の計算で結果を求める（DD-19 の試作と同じ方式）。

1. A の各要素について、B の中から同じ番号のラベルを探す（`RemoveMatch`）。
2. 見つかったら、空間が同じか（`SameSpace`、ERR-07）、座標系が同じか（`SameChart`、ERR-08）、上下が逆か（`Opposite`、ERR-05）を確かめ、その組を取り除く。
3. 見つからなければ、その要素を結果に残す。
4. B の残りを結果の後ろにつなぐ。
5. 結果の中に同じラベルが残っていたら（3 回以上の出現、または同じ側での重複）、`Distinct`（ERR-06）で止める。

1 つの式の中の縮約（`t.idx((a, a))`、トレース）は、ラベルの並びの中で同じ番号の組を探し、同じ規則を当てはめる。

### 縮約の値の計算

型の計算は、A と B の各枠の「役割」（結果の何番目の枠か、何番目の縮約の組か）を返す関数も生成する。値の計算は、次の一般的な手順で行う。

1. 結果の平らな番号 $`k`$ ごとに、結果の多重添字を求める（各枠の $`n`$ で割った商と余り）。
2. 縮約する組の多重添字 $`j`$ を回し、A と B の平らな番号を、役割と各枠の歩幅から組み立てる。
3. $`\sum_j A_{(k, j)} B_{(k, j)}`$ を置く。

役割と歩幅は型から決まる関数なので、最適化でインライン展開され、定数になる。試作では、4 次元の $`g_{ab}V^b`$ が手書きの `f64` の計算と同じ速さだった（DD-19）。

結果のテンソル性：A と B の両方が `Tensorial` なら `Tensorial`、どちらかが `ChartBound` なら `ChartBound`。結果の対称性は `NoSym`（E04）。結果の次元は $`[A][B]`$。

`expr.into_tensor()` は、残った枠をラベルの並びの順のまま持つテンソルを返す。

## 3.7 DD-07 対称性（E04、A01）

- 初期版で対称・反対称を持てるのは、階数 2 で、2 つの枠が同じ（`SameSlot`、ERR-12）テンソルだけとする。計量テンソル、Ricci テンソル、Einstein テンソル、電磁場テンソルがこの形である。階数 3 以上の対称性と、一部の添字だけの対称性は将来（VOL-06）。
- 対称・反対称のテンソルも、成分は密な入れ子の配列で持つ（A01、初期）。`from_independent` は独立な成分（対称は上三角と対角の $`n(n+1)/2`$ 個、反対称は上三角の $`n(n-1)/2`$ 個、型は `Space::SymArr`・`Space::AntiArr`）を受け取り、$`T_{ji} = \pm T_{ij}`$ を埋めて置く。
- `set` を公開しないので、作ったあとで対称性が壊れることはない。演算の結果の対称性は 2b §CMP-E04 の表に従う。
- `symmetrize`：$`\tfrac{1}{2}(T_{ab} + T_{ba})`$、`antisymmetrize`：$`\tfrac{1}{2}(T_{ab} - T_{ba})`$。

## 3.8 DD-08 計量テンソル・座標変換・逆行列（E03、E02、U04）

```rust
pub trait Metric {                       // EXT-06
    type Chart: Chart;
    type Dim;                            // 全体の次元 [g]
    fn components<N: Scalar>(&self, p: &Point<Self::Chart, N>)
        -> <<Self::Chart as Chart>::Space as Space>::SymArr<N>;   // 独立な成分
    fn inverse<N: Scalar>(&self, _p: &Point<Self::Chart, N>)
        -> Option<<<Self::Chart as Chart>::Space as Space>::SymArr<N>> { None }
}
pub trait ChartTransform {               // EXT-05
    type From: Chart;
    type To: Chart;
    fn map<N: Scalar>(&self, p: &Point<Self::From, N>) -> Point<Self::To, N>;
    fn jacobian<N: Scalar>(&self, p: &Point<Self::From, N>) -> /* ∂x'^a/∂x^b の n'×n の入れ子の配列 */;
    fn inverse_jacobian<N: Scalar>(&self, _p: &Point<Self::From, N>) -> Option</* 同 */> { None }
}
```

- 計量テンソルは対称なので、独立な成分だけを返させる（DD-07）。
- 逆（$`g^{ab}`$、Jacobi 行列の逆）が与えられないときは、U04 で求める。

### 逆行列（U04）

| 大きさ $`n`$ | 方式 |
| --- | --- |
| 1 | $`1/a`$ |
| 2〜4 | 余因子行列 $`A^{-1} = \operatorname{adj}(A)/\det A`$。4 次は 2 次の小行列式を使う Laplace 展開で、$`\det`$ と余因子を共通の部分から求める |
| 5 以上 | 部分ピボット選択つき LU 分解。単位行列を右辺として解く |

- 方式は `Space::N`（定数）で選ぶ。`match` の分岐は定数なので、最適化で消える。
- 正則でないことの判定：Hadamard の不等式 $`\lvert \det A \rvert \le \prod_i \lVert a_i \rVert`$（$`a_i`$ は行）から、比 $`\rho = \lvert \det A \rvert / \prod_i \lVert a_i \rVert \in [0, 1]`$ は行の大きさによらない。$`\rho < n\,\varepsilon`$（$`\varepsilon`$ は `f64` の計算機イプシロン）のとき、数値的に正則でないとして `Singular` を返す。LU 分解ではピボットの絶対値を、その列の最大値と比べて同じ基準で判定する。しきい値は 04 のテストで、正則に近い行列（条件数 $`10^{8}`$〜$`10^{15}`$）を使って確かめ、必要なら直す（決まり 12）。
- 自動微分の数値の型でも同じ式を使う。判定は値の部分（`re`）だけで行う。

## 3.9 DD-09 数値の型と自動微分（U02、E08）

```rust
pub trait Scalar: Copy + Add<Output = Self> + Sub<…> + Mul<…> + Div<…> + Neg<…> + 'static {
    fn from_f64(x: f64) -> Self;
    fn re(self) -> f64;                  // 値の部分
    fn sqrt(self) -> Self; fn exp(self) -> Self; fn ln(self) -> Self;
    fn sin(self) -> Self; fn cos(self) -> Self; fn tan(self) -> Self;
    fn sinh(self) -> Self; fn cosh(self) -> Self; fn tanh(self) -> Self;
    fn asin(self) -> Self; fn acos(self) -> Self; fn atan(self) -> Self;
    fn powi(self, k: i32) -> Self; fn powf(self, p: f64) -> Self; fn abs(self) -> Self;
}
```

| 型 | 持つもの | 用途 |
| --- | --- | --- |
| `f64` | 値 | 通常の計算 |
| `Dual<Sp>` | 値 $`v`$、1 階微分 $`\partial_c v`$（`Sp::Arr<f64>`） | Christoffel 記号、共変微分、Jacobi 行列 |
| `Dual2<Sp>` | 値、1 階微分、2 階微分 $`\partial_c \partial_d v`$（`Sp::SymArr<f64>`） | Riemann テンソル |

- 演算は連鎖律で実装する。例：$`f(v)`$ の 2 階微分は $`\partial_c\partial_d f = f''(v)\,\partial_c v\,\partial_d v + f'(v)\,\partial_c\partial_d v`$。
- 数学関数は std の `f64` のものを使う（NFR-01）。
- E08 `gradient(場, 点)`：点の各座標 $`x^c`$ を、値 $`x^c`$・1 階微分 $`e_c`$（$`c`$ 番目だけ 1）の `Dual` にして場を 1 回評価する。`hessian` は `Dual2` で同様に 1 回評価する（2 階微分の初期値は 0）。
- 費用：$`n = 4`$ の `Dual2` は 1 つの数に 15 個の `f64` を持つ。計量テンソル（10 成分）を 1 回評価するだけなので、問題にならない大きさである。

## 3.10 DD-10 接続と曲率、共変微分（E09、M03）

M03 は、点 $`p`$ で計量テンソルを `Dual2` で 1 回評価し、$`g_{ab}`$、$`\partial_c g_{ab}`$、$`\partial_c\partial_d g_{ab}`$ を得る（API-53 `curvature` はこれを 1 回だけ行い、すべてを返す）。

| 量 | 式（E09 が計算する） |
| --- | --- |
| Christoffel 記号 | $`\Gamma^a{}_{bc} = \tfrac{1}{2} g^{ad}(\partial_b g_{dc} + \partial_c g_{db} - \partial_d g_{bc})`$ |
| 逆の微分 | $`\partial_e g^{ad} = -g^{af}(\partial_e g_{fh})\,g^{hd}`$ |
| Christoffel 記号の微分 | $`\partial_e \Gamma^a{}_{bc} = \tfrac{1}{2}(\partial_e g^{ad})(\cdots) + \tfrac{1}{2} g^{ad}(\partial_e\partial_b g_{dc} + \partial_e\partial_c g_{db} - \partial_e\partial_d g_{bc})`$ |
| Riemann テンソル | $`R^a{}_{bcd} = \partial_c \Gamma^a{}_{db} - \partial_d \Gamma^a{}_{cb} + \Gamma^a{}_{ce}\Gamma^e{}_{db} - \Gamma^a{}_{de}\Gamma^e{}_{cb}`$ |

| 量 | 求め方（M03 が E02 の縮約で） | 型 |
| --- | --- | --- |
| Ricci テンソル | $`R_{bd} = R^a{}_{bad}`$ | 階数 2、`Symmetric`（Levi-Civita 接続では対称。値の上でも対称になることを 04 で確かめる） |
| スカラー曲率 | $`R = g^{bd} R_{bd}`$ | 次元 $`1/[g]`$ |
| Einstein テンソル | $`G_{ab} = R_{ab} - \tfrac{1}{2} R\, g_{ab}`$ | 階数 2、`Symmetric` |

- 型：Christoffel 記号は `Tensor<Dimensionless, (Up<C>, Down<C>, Down<C>), NoSym, ChartBound>`（API-50）。Riemann テンソルは同じく 4 つの枠を持つ `Tensorial` なテンソル。全体の次元はどちらも無次元で、成分の次元は DD-04 の式で決まる。
- 共変微分 $`\nabla_c T`$：テンソル場を `Dual` で評価して $`\partial_c T`$（`ChartBound`）を得て、上付きの枠ごとに $`+\Gamma^{a}{}_{c e} T^{\cdots e \cdots}`$、下付きの枠ごとに $`-\Gamma^{e}{}_{c b} T_{\cdots e \cdots}`$ を足す。足し終えたものを M03 の中で `Tensorial` として返す。テンソル性を付け替えるのはこの関数だけで、`pub(crate)` の関数で行う。
- 共変微分の結果の枠の並びは、場の枠の並びの後ろに `Down<C>` を 1 つつないだもの。

## 3.11 DD-11 数値解法の状態（U05）

```rust
pub trait VectorSpace: Copy {
    type Raw: RawVec;                                    // 次元を外した数値の入れ物（同じ入れ子の配列）
    type Per<D>: VectorSpace;                            // 各要素の次元を D で割った型
    type Times<D>: VectorSpace;                          // 各要素の次元に D を掛けた型
    fn to_raw(&self) -> Self::Raw;
    fn from_raw(r: &Self::Raw) -> Self;
}
pub trait RawVec: Copy {                                 // E05〜E08 が使う、次元のない入れ物
    const LEN: usize;
    fn zero() -> Self;
    fn get(&self, i: usize) -> f64;
    fn set(&mut self, i: usize, v: f64);
    fn axpy(&mut self, a: f64, x: &Self);                // self += a x
    type Mat: RawMat;                                    // LEN × LEN の行列（Jacobi 行列用）
}
```

- `Tensor` は `Raw = S::Store<f64>` で実装する。`Point` は座標の値を並べたものとして実装する（測地線の状態に使う）。
- 要素が 2〜6 個の組（`(A, B)` など）は、要素ごとに実装する（`Raw = (A::Raw, B::Raw)`、`get` は前から順に番号を振る）。
- `RawVec::Mat` は `Raw` の各要素を `Raw` に置き換えた入れ物（入れ子の配列の入れ子）で、要素 $`(i, j)`$ の読み書きができる。Rosenbrock 法の線形方程式は、この行列を U04 に渡して解く（大きさ 5 以上は LU 分解）。
- `f` の戻り値が `Y::Per<T>` でなければ、`solve_ivp` の型の制約 `RhsReturns<Y, T>`（ERR-14）で止める。`atol` の型が `Y` でなければ ERR-15。

## 3.12 DD-12 常微分方程式の解法（E05、M02）

| 解法 | 種類 | 次数（誤差の見積もり） | 係数の出典 |
| --- | --- | --- | --- |
| `Rk4 { h }` | 陽的、固定刻み | 4（なし） | 古典的な Runge–Kutta 法 |
| `DormandPrince54` | 陽的、埋め込み、適応刻み、FSAL | 5（4） | Dormand & Prince (1980) |
| `Leapfrog { h }` | シンプレクティック、2 階の方程式、固定刻み | 2（なし） | kick–drift–kick 形式 |
| `Rosenbrock` | 線形陰的、L 安定、埋め込み、適応刻み | 4（3） | Hairer & Wanner, *Solving ODEs II*（1996）の RODAS4 |

- 係数は、出典の表をそのまま定数として書き、出典をコメントに書く。
- Rosenbrock 法の 1 歩：$`W = I - h\gamma J`$ を 1 回 LU 分解し（U04）、各段の線形方程式を同じ分解で解く。$`J`$ は利用者が与えるか、右辺が数値の型によらない形なら `Dual` で、そうでなければ差分（DD-15）で求める。
- 時間の向きは負でもよい（$`t_1 < t_0`$）。
- `solve_ivp_with` の `on_step` は、受け入れた各歩の $`(t, y)`$ を受け取る。密出力（任意の時刻の値の補間）は将来。

## 3.13 DD-13 誤差の測り方と刻み幅の調整（E06）

- 誤差のノルム（Hairer らの方式）：

```math
\mathrm{err} = \sqrt{\frac{1}{L}\sum_{i=1}^{L}\left(\frac{e_i}{\mathrm{atol}_i + \mathrm{rtol}\,\max(\lvert y_i \rvert, \lvert y'_i \rvert)}\right)^2}
```

  $`\mathrm{atol}_i`$ は、状態と同じ型の `atol` を `to_raw` した値なので、成分ごとに次元が合っている。
- 受け入れ：$`\mathrm{err} \le 1`$。
- 次の刻み：$`h_\mathrm{new} = h \cdot \min(f_\mathrm{max}, \max(f_\mathrm{min}, 0.9\,\mathrm{err}^{-1/(q+1)}))`$。$`q`$ は埋め込み公式の低いほうの次数。$`f_\mathrm{max} = 5`$、$`f_\mathrm{min} = 0.2`$。直前の歩を棄却したときは $`f_\mathrm{max} = 1`$。PI 制御は将来（VOL-08）。
- 最初の刻み：Hairer らの初期刻みの推定（右辺を 2 回評価する）。
- 失敗：歩数が上限（既定 100,000）を超えたら `MaxSteps`。$`\lvert h \rvert < 16\,\varepsilon\,\lvert t \rvert`$ になったら `StepTooSmall`。値が有限でなくなったら `NonFinite`。
- 適応型の積分：区間ごとの誤差の見積もりの和が $`\max(\mathrm{atol}, \mathrm{rtol}\,\lvert I \rvert)`$ 以下になれば受け入れ、そうでなければ誤差の最も大きい区間を選んで 2 分する（`pick_interval`）。区間の一覧は `Vec` に置く（RES-02）。区間の数の上限（既定 1,000）を超えたら `ToleranceNotMet`。

## 3.14 DD-14 定積分（E07、M02）

| 解法 | 方式 | 誤差の見積もり |
| --- | --- | --- |
| `Trapezoid { n }` | 複合台形公式（$`n`$ 区間） | $`n`$ と $`2n`$ の差 |
| `Simpson { n }` | 複合 Simpson 公式（$`n`$ は偶数。奇数なら型ではなく、作るときに切り上げる） | $`n`$ と $`2n`$ の差 |
| `Romberg` | 台形公式の区間を倍々にし、Richardson 補外の表を作る（最大 20 段） | 表の対角の隣どうしの差 |
| `GaussLegendre { n }` | $`n`$ 点の Gauss–Legendre 求積。節点は Legendre 多項式の根を Newton 法で求める（初期値 $`\cos(\pi(i - 1/4)/(n + 1/2))`$） | $`n`$ と $`n + 1`$ の差 |
| `GaussKronrod` | 7 点 Gauss・15 点 Kronrod（節点と重みは QUADPACK の値）の適応型 | Gauss と Kronrod の差（QUADPACK の補正式） |

- 被積分関数の値は `VectorSpace` の型でよい（ベクトル値の積分）。結果は `Y::Times<X>`。

## 3.15 DD-15 差分による数値微分（E08）

数値の型によらない形で書かれていない関数（`Fn(f64) -> f64` 相当）に使う。

- 中心差分：$`f'(x) \approx (f(x + h) - f(x - h)) / 2h`$、$`h = \varepsilon^{1/3}\max(\lvert x \rvert, x_\mathrm{scale})`$。$`h`$ は丸め誤差（$`\varepsilon/h`$）と打ち切り誤差（$`h^2`$）の釣り合いから決める。
- Richardson 補外：Ridders の方法（$`h`$ を一定の比で小さくしながら補外し、誤差の見積もりが増え始めたら止める）。誤差の見積もりも返す。
- Jacobi 行列：列ごとに中心差分。

## 3.16 DD-16 エラー文の実装（U03）

ERR ごとに専用の trait を 1 つ置き、`#[diagnostic::on_unimplemented]` をつける。

| ERR | trait | 実装されるとき |
| --- | --- | --- |
| ERR-01 | `SameDimension<D2>` | `D2 = Self` |
| ERR-02 | `UnitMatches<D>` | 単位の次元 = `D` |
| ERR-03 | `SameSlots<S2>` | `S2 = Self` |
| ERR-04 | `LabelCount<S>` | ラベルの並びの長さ = 枠の並びの長さ |
| ERR-05 | `Opposite<V2, L>` | 上と下、下と上 |
| ERR-06 | `Distinct<L>` | 結果の中にラベル `L` が残っていない |
| ERR-07 | `SameSpace<Sp2, L>` | 同じ空間 |
| ERR-08 | `SameChart<C2, L>` | 同じ座標系 |
| ERR-09 | `MetricOn<C>` | 計量テンソルの座標系 = `C` |
| ERR-10 | `TransformFrom<C>` | 座標変換の変換元 = `C` |
| ERR-11 | `CoordOf<C>` | 座標の印が `C::Coords` にある |
| ERR-12 | `SameSlot<S2>` | 2 つの枠が同じ |
| ERR-13 | `SameSymmetry<Sym2>` | 同じ対称性 |
| ERR-14 | `RhsReturns<Y, T>` | 戻り値の型 = `Y::Per<T>` |
| ERR-15 | `ToleranceFor<Y>` | `atol` の型 = `Y` |
| ERR-16 | `IsTensor` | `Tensorial` |
| ERR-17 | `SameKind<K2>` | 同じテンソル性 |

- 文は 2a §2a.4 のとおり。`{L}` などには、ラベル・座標系・単位の型の名前が入るように、その型を trait の型引数にする。
- 各 ERR について、コンパイル失敗のテスト（`tests/ui/*.rs` と、期待する出力 `*.stderr`）を置く（04 の TC-F）。期待する出力は見出しの行だけを比べる（下に続く補足の内部の型の名前は、コンパイラの版で変わるため）。

## 3.17 DD-17 失敗の型（M02、U04）

```rust
pub enum SolveError<T> { MaxSteps { t: T }, StepTooSmall { t: T, h: T }, NonFinite { t: T }, Singular { t: T } }
pub enum QuadError { ToleranceNotMet { estimate_error: f64 }, NonFinite }
pub enum LinearError { Singular }
```

- `T` は時間の型（`Quantity`）。失敗した時刻を、利用者の単位で取り出せる。
- `core::error::Error` と `Display` を実装する。メッセージは英語（NFR-07）。
- 幾何の操作（M03）で計量テンソルが正則でないときは、`Result<_, LinearError>` を返す。

## 3.18 DD-18 利用者向けのマクロ

| マクロ | 生成するもの |
| --- | --- |
| `space!(Spacetime, 4)` | 空の構造体と `Space` の実装（`Arr<T> = [T; 4]`、`SymArr<T> = [T; 10]`、`AntiArr<T> = [T; 6]` など、長さは数字から展開時に計算する） |
| `chart!(Polar in Plane: r: Length, theta: Dimensionless)` | 座標系の構造体、座標の印の構造体（`R`、`Theta`）、`Chart` の実装 |
| `unit_system!(…)` | 単位系の構造体と次元の別名 |
| `label!(kappa = 40)` | ラベルの構造体と `Label` の実装 |

- マクロは `macro_rules!` で書く（手続き型マクロのための外部クレートに依存しない。NFR-01）。
- `space!` の配列の長さ（$`n(n+1)/2`$ など）は、マクロの中で `{ $n * ($n + 1) / 2 }` という定数式として書く。これは具体的な数の式なので、stable で使える。

## 3.19 DD-19 試作で確かめたこと

決まり 4 に従い、詳細設計の前提を stable の Rust 1.99.0 で試した。

| 前提 | 結果 |
| --- | --- |
| 型レベルのラベルの照合で縮約の相手を見つけ、結果の枠の並びを求められる（2a §2a.1） | できた |
| 空間ごとの配列の型（GAT）を組み合わせて、入れ子の配列の入れ物を型の計算で作れる | できた。時空 × 場空間で `[[f64; 2]; 4]`（64 バイト）、4 次元の Riemann テンソルで 2048 バイト（成分の数 × 8 バイト。NFR-03） |
| `unsafe` なしで、平らな番号で成分を読み書きできる | できた |
| 一般的な縮約の手順が、手書きの `f64` の計算と同じ速さで動く | 4 次元の $`g_{ab}V^b`$ を 5,000 万回計算して、時間の比は 0.97〜1.00（release、`codegen-units = 1`）。NFR-02 の基準は 04 で、より多くの演算で計測して決める |
| `no_std` の `core` に数学関数がある | なかった（`core::f64::math` は unstable）。std を前提にした（NFR-01、OPEN-09） |

## 3.20 DD-20 ビルドと CI

- workspace の `Cargo.toml` に `[profile.release] codegen-units = 1` と `lto = "thin"` を置く（ベンチマークの値を安定させるため）。
- CI（`.github/workflows/ci.yml`）の手順：
  1. `check_workflows.py`、`check_privacy.py`、`check_layers.py`（NFR-11）
  2. `cargo fmt --check`、`cargo clippy --workspace --all-targets -- -D warnings`
  3. `cargo test --workspace`（コンパイル失敗のテストを含む）
  4. `cargo run --example …`（すべての例題。結果を解析解と比べて、外れたら失敗で終わる）
  5. `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps`
- ベンチマーク（`benches/`）は CI では実行しない（実行環境で時間がばらつくため）。手元で計測し、結果を 04 に記録する。
