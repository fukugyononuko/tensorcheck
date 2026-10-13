//! Connection, curvature and covariant derivatives (CMP-M03, DD-10, CC-09).
//!
//! Every operation evaluates the metric tensor once with second-order dual numbers
//! at the point, which gives `g`, `∂g` and `∂∂g` exactly (up to rounding); the
//! connection engine turns them into Christoffel symbols and the Riemann tensor, and
//! the index engine contracts the result.

use crate::access::component::Slots;
use crate::contracts::{
    Chart, ChartBound, Down, Ix, NoSym, Point, Rank2, Rank3, Rank4, Space, Symmetric, Tensor,
    Tensorial, Up,
};
use crate::engine::connection::{
    christoffel as gamma_of, christoffel_derivative, riemann as riemann_of,
};
use crate::engine::differentiation::{SpaceOf, seed1, seed2, unpack1, unpack2};
use crate::engine::dimension::{DimDiv, DimMul, DimPow, Pow, Prod, Quot};
use crate::engine::index::contract;
use crate::engine::metric::Metric;
use crate::labels::{a as La, b as Lb, c as Lc};
use crate::utility::diagnostics::{IsCoordinateBasis, MetricOn};
use crate::utility::extent::{Extent, sym_index};
use crate::utility::linear::{LinearError, inverse};
use crate::utility::scalar::Scalar;
use crate::utility::typelevel::{Concat, Cons, N1, N2, Nil};

/// The reference dimension `[X]` of the space of chart `C`.
pub type RefDimOf<C> = <<C as Chart>::Space as Space>::RefDim;

/// Christoffel symbols `Γ^a_bc` in chart `C`: chart-bound (not a tensor), overall
/// dimension `1/[X]` (API-50).
pub type Christoffel<C> = Rank3<Pow<RefDimOf<C>, N1>, Up<C>, Down<C>, Down<C>, NoSym, ChartBound>;
/// The Riemann tensor `R^a_bcd`, overall dimension `1/[X]²`.
pub type Riemann<C> = Rank4<Pow<RefDimOf<C>, N2>, Up<C>, Down<C>, Down<C>, Down<C>>;
/// The Ricci tensor `R_ab`, overall dimension `1/[X]²`.
pub type Ricci<C> = Rank2<Pow<RefDimOf<C>, N2>, Down<C>, Down<C>, Symmetric>;
/// The Einstein tensor `G_ab = R_ab − ½ R g_ab`, overall dimension `1/[X]²`.
pub type Einstein<C> = Ricci<C>;
/// The scalar curvature `R = g^ab R_ab` of a metric tensor `G`: dimension `1/([g][X]²)`.
pub type RicciScalar<G> = crate::contracts::Quantity<
    Prod<Pow<<G as Metric>::Dim, N1>, Pow<RefDimOf<<G as Metric>::Chart>, N2>>,
>;

/// Everything computed from one evaluation of the metric tensor at a point (API-53).
pub struct Curvature<G: Metric>
where
    RefDimOf<G::Chart>: DimPow<N1> + DimPow<N2>,
    G::Dim: DimPow<N1>,
    Pow<G::Dim, N1>: DimMul<Pow<RefDimOf<G::Chart>, N2>>,
{
    /// `Γ^a_bc`.
    pub christoffel: Christoffel<G::Chart>,
    /// `R^a_bcd`.
    pub riemann: Riemann<G::Chart>,
    /// `R_ab`.
    pub ricci: Ricci<G::Chart>,
    /// `R`.
    pub ricci_scalar: RicciScalar<G>,
    /// `G_ab`.
    pub einstein: Einstein<G::Chart>,
}

/// `g`, `∂g`, `∂∂g` and `g^{-1}` at a point, as flat arrays.
struct Jets {
    n: usize,
    g: Vec<f64>,
    ginv: Vec<f64>,
    dg: Vec<f64>,
    ddg: Vec<f64>,
}

fn jets<G: Metric>(g: &G, p: &Point<G::Chart>) -> Result<Jets, LinearError> {
    let n = <SpaceOf<G::Chart> as Extent>::N;
    let sym = g.components(&seed2(p));
    // Expand the independent components to full n × n before unpacking.
    let full: Vec<_> = (0..n * n)
        .map(|k| sym.as_ref()[sym_index(k / n, k % n, n)])
        .collect();
    let (gv, dg, ddg) = unpack2(&full);
    let mut ginv = vec![0.0; n * n];
    inverse(&gv, n, &mut ginv)?;
    Ok(Jets {
        n,
        g: gv,
        ginv,
        dg,
        ddg,
    })
}

fn ricci_values<C: Chart>(riem: &[f64]) -> Vec<f64> {
    // R_bd = R^a_bad: contract slots 0 and 2 (label a), keep b and c.
    type L<C> = Cons<
        Ix<La, Up<C>>,
        Cons<Ix<Lb, Down<C>>, Cons<Ix<La, Down<C>>, Cons<Ix<Lc, Down<C>>, Nil>>>,
    >;
    type S4<C> = Cons<Up<C>, Cons<Down<C>, Cons<Down<C>, Cons<Down<C>, Nil>>>>;
    type S2<C> = Cons<Down<C>, Cons<Down<C>, Nil>>;
    let store = S4::<C>::from_fn(|k| riem[k]);
    let raw = contract::<S4<C>, Nil, S2<C>, L<C>, Nil, f64>(&store, &1.0);
    let n = <SpaceOf<C> as Extent>::N;
    let r: Vec<f64> = (0..n * n).map(|k| S2::<C>::get(&raw, k)).collect();
    // Symmetric in theory (Levi-Civita); symmetrize so the marker holds exactly.
    (0..n * n)
        .map(|k| 0.5 * (r[k] + r[(k % n) * n + k / n]))
        .collect()
}

fn scalar_value<C: Chart>(ginv: &[f64], ricci: &[f64]) -> f64 {
    // R = g^ab R_ab, contracted by the index engine.
    type SU<C> = Cons<Up<C>, Cons<Up<C>, Nil>>;
    type SD<C> = Cons<Down<C>, Cons<Down<C>, Nil>>;
    type LU<C> = Cons<Ix<La, Up<C>>, Cons<Ix<Lb, Up<C>>, Nil>>;
    type LD<C> = Cons<Ix<La, Down<C>>, Cons<Ix<Lb, Down<C>>, Nil>>;
    let gs = SU::<C>::from_fn(|k| ginv[k]);
    let rs = SD::<C>::from_fn(|k| ricci[k]);
    contract::<SU<C>, SD<C>, Nil, LU<C>, LD<C>, f64>(&gs, &rs)
}

/// Christoffel symbols `Γ^a_bc` of `g` at `p` (API-51). Fails if the metric tensor
/// is singular at `p`. The chart must be a coordinate basis (ERR-19).
pub fn christoffel<G: Metric>(
    g: &G,
    p: &Point<G::Chart>,
) -> Result<Christoffel<G::Chart>, LinearError>
where
    <G::Chart as Chart>::Basis: IsCoordinateBasis<G::Chart>,
    RefDimOf<G::Chart>: DimPow<N1>,
{
    let n = <SpaceOf<G::Chart> as Extent>::N;
    let sym = g.components(&seed1(p));
    let full: Vec<_> = (0..n * n)
        .map(|k| sym.as_ref()[sym_index(k / n, k % n, n)])
        .collect();
    let (gv, dg) = unpack1(&full);
    let mut ginv = vec![0.0; n * n];
    inverse(&gv, n, &mut ginv)?;
    let gamma = gamma_of(n, &ginv, &dg);
    Ok(Tensor::from_flat_fn(|i| gamma[i]))
}

/// The Riemann, Ricci and Einstein tensors and the scalar curvature of `g` at `p`,
/// with the Christoffel symbols, from one evaluation of the metric tensor (API-53).
pub fn curvature<G: Metric>(g: &G, p: &Point<G::Chart>) -> Result<Curvature<G>, LinearError>
where
    <G::Chart as Chart>::Basis: IsCoordinateBasis<G::Chart>,
    RefDimOf<G::Chart>: DimPow<N1> + DimPow<N2>,
    G::Dim: DimPow<N1>,
    Pow<G::Dim, N1>: DimMul<Pow<RefDimOf<G::Chart>, N2>>,
{
    let j = jets(g, p)?;
    let n = j.n;
    let gamma = gamma_of(n, &j.ginv, &j.dg);
    let dgamma = christoffel_derivative(n, &j.ginv, &j.dg, &j.ddg);
    let r = riemann_of(n, &gamma, &dgamma);
    let riemann: Riemann<G::Chart> = Tensor::from_flat_fn(|i| r[i]);
    let ric = ricci_values::<G::Chart>(&r);
    let scalar = scalar_value::<G::Chart>(&j.ginv, &ric);
    let ein: Vec<f64> = (0..n * n).map(|k| ric[k] - 0.5 * scalar * j.g[k]).collect();
    Ok(Curvature {
        christoffel: Tensor::from_flat_fn(|i| gamma[i]),
        riemann,
        ricci: Tensor::from_flat_fn(|i| ric[i]),
        ricci_scalar: Tensor::from_flat_fn(|_| scalar),
        einstein: Tensor::from_flat_fn(|i| ein[i]),
    })
}

/// The Riemann tensor `R^a_bcd` of `g` at `p` (API-52).
pub fn riemann<G: Metric>(g: &G, p: &Point<G::Chart>) -> Result<Riemann<G::Chart>, LinearError>
where
    <G::Chart as Chart>::Basis: IsCoordinateBasis<G::Chart>,
    RefDimOf<G::Chart>: DimPow<N1> + DimPow<N2>,
    G::Dim: DimPow<N1>,
    Pow<G::Dim, N1>: DimMul<Pow<RefDimOf<G::Chart>, N2>>,
{
    Ok(curvature(g, p)?.riemann)
}

/// The Ricci tensor `R_ab` of `g` at `p` (API-52).
pub fn ricci<G: Metric>(g: &G, p: &Point<G::Chart>) -> Result<Ricci<G::Chart>, LinearError>
where
    <G::Chart as Chart>::Basis: IsCoordinateBasis<G::Chart>,
    RefDimOf<G::Chart>: DimPow<N1> + DimPow<N2>,
    G::Dim: DimPow<N1>,
    Pow<G::Dim, N1>: DimMul<Pow<RefDimOf<G::Chart>, N2>>,
{
    Ok(curvature(g, p)?.ricci)
}

/// The scalar curvature `R` of `g` at `p` (API-52).
pub fn ricci_scalar<G: Metric>(g: &G, p: &Point<G::Chart>) -> Result<RicciScalar<G>, LinearError>
where
    <G::Chart as Chart>::Basis: IsCoordinateBasis<G::Chart>,
    RefDimOf<G::Chart>: DimPow<N1> + DimPow<N2>,
    G::Dim: DimPow<N1>,
    Pow<G::Dim, N1>: DimMul<Pow<RefDimOf<G::Chart>, N2>>,
{
    Ok(curvature(g, p)?.ricci_scalar)
}

/// The Einstein tensor `G_ab` of `g` at `p` (API-52).
pub fn einstein<G: Metric>(g: &G, p: &Point<G::Chart>) -> Result<Einstein<G::Chart>, LinearError>
where
    <G::Chart as Chart>::Basis: IsCoordinateBasis<G::Chart>,
    RefDimOf<G::Chart>: DimPow<N1> + DimPow<N2>,
    G::Dim: DimPow<N1>,
    Pow<G::Dim, N1>: DimMul<Pow<RefDimOf<G::Chart>, N2>>,
{
    Ok(curvature(g, p)?.einstein)
}

// ---------------------------------------------------------------------------
// Tensor fields and covariant derivatives

/// A tensor field: a tensor at every point of a chart (EXT-08).
///
/// Write [`TensorField::at`] generically over the scalar type `N` (build the result
/// with `Tensor::from_base_components`), so the library can differentiate it.
pub trait TensorField {
    /// The chart.
    type Chart: Chart;
    /// The dimension of the field.
    type Dim;
    /// The index slots of the field.
    type Slots: Slots;
    /// The field at `p`.
    fn at<N: Scalar>(
        &self,
        p: &Point<Self::Chart, N>,
    ) -> Tensor<Self::Dim, Self::Slots, NoSym, Tensorial, N>;
}

/// Every slot is in chart `C` (the chart of the metric tensor `G`).
pub trait AllOnChart<C, G> {}
impl<C, G> AllOnChart<C, G> for Nil {}
impl<C, G, Cs, T: AllOnChart<C, G>> AllOnChart<C, G> for Cons<Up<Cs>, T> where C: MetricOn<Cs, G> {}
impl<C, G, Cs, T: AllOnChart<C, G>> AllOnChart<C, G> for Cons<Down<Cs>, T> where C: MetricOn<Cs, G> {}

/// The covariant derivative `∇_c T^{a…}_{b…}` of a tensor field at `p`, with the
/// Levi-Civita connection of `g` (API-54). The derivative index is appended as the
/// last (lower) index; the dimension is `[T]/[X]`.
pub fn covariant_derivative<F, G>(
    field: &F,
    g: &G,
    p: &Point<G::Chart>,
) -> Result<
    Tensor<
        Quot<F::Dim, RefDimOf<G::Chart>>,
        <F::Slots as Concat<Cons<Down<G::Chart>, Nil>>>::Out,
        NoSym,
        Tensorial,
    >,
    LinearError,
>
where
    F: TensorField<Chart = G::Chart>,
    G: Metric,
    <G::Chart as Chart>::Basis: IsCoordinateBasis<G::Chart>,
    F::Slots: AllOnChart<G::Chart, G>
        + crate::engine::index::SlotVariances
        + Concat<Cons<Down<G::Chart>, Nil>>,
    <F::Slots as Concat<Cons<Down<G::Chart>, Nil>>>::Out: Slots,
    F::Dim: DimDiv<RefDimOf<G::Chart>>,
{
    let n = <SpaceOf<G::Chart> as Extent>::N;
    // ∂_c T from one evaluation with first-order duals.
    let t = field.at(&seed1(p));
    let len = <F::Slots as Slots>::LEN;
    let rank = <F::Slots as Slots>::RANK;
    let vals: Vec<_> = (0..len).map(|k| t.flat(k)).collect();
    let (tv, dt) = unpack1(&vals);
    // Γ at p.
    let sym = g.components(&seed1(p));
    let full: Vec<_> = (0..n * n)
        .map(|k| sym.as_ref()[sym_index(k / n, k % n, n)])
        .collect();
    let (gv, dg) = unpack1(&full);
    let mut ginv = vec![0.0; n * n];
    inverse(&gv, n, &mut ginv)?;
    let gamma = gamma_of(n, &ginv, &dg);
    let gm = |a: usize, b: usize, c: usize| gamma[(a * n + b) * n + c];

    let mut ups = [false; crate::engine::index::MAX_RANK];
    <F::Slots as crate::engine::index::SlotVariances>::write(&mut ups[..rank]);
    let mut stride = vec![1usize; rank];
    for s in (0..rank.saturating_sub(1)).rev() {
        stride[s] = stride[s + 1] * n;
    }
    let mut out = vec![0.0; len * n];
    for k in 0..len {
        for cc in 0..n {
            let mut v = dt[cc * len + k];
            for s in 0..rank {
                let idx = (k / stride[s]) % n;
                for e in 0..n {
                    let ke = k - idx * stride[s] + e * stride[s];
                    if ups[s] {
                        v += gm(idx, cc, e) * tv[ke];
                    } else {
                        v -= gm(e, cc, idx) * tv[ke];
                    }
                }
            }
            out[k * n + cc] = v;
        }
    }
    Ok(Tensor::from_flat_fn(|i| out[i]))
}
