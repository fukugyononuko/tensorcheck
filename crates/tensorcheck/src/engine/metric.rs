//! Metric tensors and their inverses (CMP-E03, DD-08).

use crate::access::component::{Extent, sym_index};
use crate::contracts::{Chart, Point};
use crate::utility::linear::{LinearError, inverse};
use crate::utility::scalar::Scalar;

/// The array of independent components of a symmetric rank-2 tensor on chart `C`'s space.
pub type SymComponents<C, N> = <<C as Chart>::Space as Extent>::SymArr<N>;

/// A metric tensor `g_ab(x)` on a chart (EXT-06).
///
/// Write the components generically over the scalar type `N` so that the library
/// can evaluate them with dual numbers to obtain derivatives (FR-16).
pub trait Metric: 'static {
    /// The chart the components refer to.
    type Chart: Chart;
    /// The dimension `[g]` of the metric tensor (dimensionless for a line element
    /// with `ds²` of dimension `[X]²`).
    type Dim;
    /// Independent components `g_ab` with `a <= b`, row by row
    /// (`g_00, g_01, …, g_0(n-1), g_11, …`), in base units.
    fn components<N: Scalar>(&self, p: &Point<Self::Chart, N>) -> SymComponents<Self::Chart, N>;
    /// Independent components of the inverse `g^ab`, if known in closed form.
    fn inverse<N: Scalar>(
        &self,
        _p: &Point<Self::Chart, N>,
    ) -> Option<SymComponents<Self::Chart, N>> {
        None
    }
}

/// Expands independent components into a full row-major `n × n` matrix.
pub fn sym_to_full<N: Copy>(sym: &[N], n: usize, out: &mut [N]) {
    for i in 0..n {
        for j in 0..n {
            out[i * n + j] = sym[sym_index(i, j, n)];
        }
    }
}

/// The full matrix `g_ab` at `p`, row-major.
pub fn metric_at<G: Metric, N: Scalar>(g: &G, p: &Point<G::Chart, N>) -> Vec<N> {
    let n = <<G::Chart as Chart>::Space as Extent>::N;
    let sym = g.components(p);
    let mut full = vec![N::zero(); n * n];
    sym_to_full(sym.as_ref(), n, &mut full);
    full
}

/// The full matrix `g^ab` at `p`, row-major: the closed form if the metric gives
/// one, otherwise the numerical inverse (U04).
pub fn inverse_at<G: Metric, N: Scalar>(
    g: &G,
    p: &Point<G::Chart, N>,
) -> Result<Vec<N>, LinearError> {
    let n = <<G::Chart as Chart>::Space as Extent>::N;
    let mut full = vec![N::zero(); n * n];
    if let Some(sym) = g.inverse(p) {
        sym_to_full(sym.as_ref(), n, &mut full);
        return Ok(full);
    }
    let gm = metric_at(g, p);
    inverse(&gm, n, &mut full)?;
    Ok(full)
}
