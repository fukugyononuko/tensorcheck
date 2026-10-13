//! Inverses and linear solves of small square matrices (CMP-U04, DD-08).
//!
//! Matrices are flat, row-major slices of length `n * n`. For `n <= 4` the inverse
//! uses the closed-form adjugate formula (no branches, so it differentiates cleanly
//! with dual numbers); for `n >= 5` it uses LU decomposition with partial pivoting.

use core::fmt;

use crate::utility::scalar::Scalar;

/// A linear-algebra failure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LinearError {
    /// The matrix is singular to working precision.
    Singular,
}

impl fmt::Display for LinearError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LinearError::Singular => f.write_str("matrix is singular to working precision"),
        }
    }
}

impl std::error::Error for LinearError {}

/// Whether `|det A| / ∏ ‖row_i‖` is below `n ε` (DD-08).
///
/// By Hadamard's inequality the ratio lies in `[0, 1]` whatever the scale of the
/// rows, so a fixed threshold works for every matrix.
fn nearly_singular<N: Scalar>(a: &[N], n: usize, det: f64) -> bool {
    let mut prod = 1.0;
    for i in 0..n {
        let norm: f64 = (0..n)
            .map(|j| a[i * n + j].re().powi(2))
            .sum::<f64>()
            .sqrt();
        if norm == 0.0 {
            return true;
        }
        prod *= norm;
    }
    let ratio = det.abs() / prod;
    !ratio.is_finite() || ratio < n as f64 * f64::EPSILON
}

#[inline(always)]
fn det2<N: Scalar>(a: N, b: N, c: N, d: N) -> N {
    a * d - b * c
}

/// Writes the inverse of the `n × n` matrix `a` into `out`.
pub fn inverse<N: Scalar>(a: &[N], n: usize, out: &mut [N]) -> Result<(), LinearError> {
    debug_assert_eq!(a.len(), n * n);
    debug_assert_eq!(out.len(), n * n);
    match n {
        0 => Ok(()),
        1 => {
            if a[0].re() == 0.0 || !a[0].re().is_finite() {
                return Err(LinearError::Singular);
            }
            out[0] = N::one() / a[0];
            Ok(())
        }
        2 => {
            let det = det2(a[0], a[1], a[2], a[3]);
            if nearly_singular(a, 2, det.re()) {
                return Err(LinearError::Singular);
            }
            out[0] = a[3] / det;
            out[1] = -a[1] / det;
            out[2] = -a[2] / det;
            out[3] = a[0] / det;
            Ok(())
        }
        3 => {
            let m = |i: usize, j: usize| a[i * 3 + j];
            let c00 = det2(m(1, 1), m(1, 2), m(2, 1), m(2, 2));
            let c01 = -det2(m(1, 0), m(1, 2), m(2, 0), m(2, 2));
            let c02 = det2(m(1, 0), m(1, 1), m(2, 0), m(2, 1));
            let det = m(0, 0) * c00 + m(0, 1) * c01 + m(0, 2) * c02;
            if nearly_singular(a, 3, det.re()) {
                return Err(LinearError::Singular);
            }
            let c10 = -det2(m(0, 1), m(0, 2), m(2, 1), m(2, 2));
            let c11 = det2(m(0, 0), m(0, 2), m(2, 0), m(2, 2));
            let c12 = -det2(m(0, 0), m(0, 1), m(2, 0), m(2, 1));
            let c20 = det2(m(0, 1), m(0, 2), m(1, 1), m(1, 2));
            let c21 = -det2(m(0, 0), m(0, 2), m(1, 0), m(1, 2));
            let c22 = det2(m(0, 0), m(0, 1), m(1, 0), m(1, 1));
            // inverse = adj / det, adj = cofactor matrix transposed
            let adj = [c00, c10, c20, c01, c11, c21, c02, c12, c22];
            for (o, x) in out.iter_mut().zip(adj) {
                *o = x / det;
            }
            Ok(())
        }
        4 => inverse4(a, out),
        _ => inverse_lu(a, n, out),
    }
}

/// 4×4 inverse by Laplace expansion in 2×2 minors of the first two and last two rows.
fn inverse4<N: Scalar>(a: &[N], out: &mut [N]) -> Result<(), LinearError> {
    let m = |i: usize, j: usize| a[i * 4 + j];
    // 2×2 minors of rows 0,1 (s) and rows 2,3 (c), columns (j,k).
    let s0 = det2(m(0, 0), m(0, 1), m(1, 0), m(1, 1));
    let s1 = det2(m(0, 0), m(0, 2), m(1, 0), m(1, 2));
    let s2 = det2(m(0, 0), m(0, 3), m(1, 0), m(1, 3));
    let s3 = det2(m(0, 1), m(0, 2), m(1, 1), m(1, 2));
    let s4 = det2(m(0, 1), m(0, 3), m(1, 1), m(1, 3));
    let s5 = det2(m(0, 2), m(0, 3), m(1, 2), m(1, 3));
    let c5 = det2(m(2, 2), m(2, 3), m(3, 2), m(3, 3));
    let c4 = det2(m(2, 1), m(2, 3), m(3, 1), m(3, 3));
    let c3 = det2(m(2, 1), m(2, 2), m(3, 1), m(3, 2));
    let c2 = det2(m(2, 0), m(2, 3), m(3, 0), m(3, 3));
    let c1 = det2(m(2, 0), m(2, 2), m(3, 0), m(3, 2));
    let c0 = det2(m(2, 0), m(2, 1), m(3, 0), m(3, 1));
    let det = s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0;
    if nearly_singular(a, 4, det.re()) {
        return Err(LinearError::Singular);
    }
    let inv = [
        m(1, 1) * c5 - m(1, 2) * c4 + m(1, 3) * c3,
        -m(0, 1) * c5 + m(0, 2) * c4 - m(0, 3) * c3,
        m(3, 1) * s5 - m(3, 2) * s4 + m(3, 3) * s3,
        -m(2, 1) * s5 + m(2, 2) * s4 - m(2, 3) * s3,
        -m(1, 0) * c5 + m(1, 2) * c2 - m(1, 3) * c1,
        m(0, 0) * c5 - m(0, 2) * c2 + m(0, 3) * c1,
        -m(3, 0) * s5 + m(3, 2) * s2 - m(3, 3) * s1,
        m(2, 0) * s5 - m(2, 2) * s2 + m(2, 3) * s1,
        m(1, 0) * c4 - m(1, 1) * c2 + m(1, 3) * c0,
        -m(0, 0) * c4 + m(0, 1) * c2 - m(0, 3) * c0,
        m(3, 0) * s4 - m(3, 1) * s2 + m(3, 3) * s0,
        -m(2, 0) * s4 + m(2, 1) * s2 - m(2, 3) * s0,
        -m(1, 0) * c3 + m(1, 1) * c1 - m(1, 2) * c0,
        m(0, 0) * c3 - m(0, 1) * c1 + m(0, 2) * c0,
        -m(3, 0) * s3 + m(3, 1) * s1 - m(3, 2) * s0,
        m(2, 0) * s3 - m(2, 1) * s1 + m(2, 2) * s0,
    ];
    for (o, x) in out.iter_mut().zip(inv) {
        *o = x / det;
    }
    Ok(())
}

/// LU decomposition with partial pivoting, in place. Returns the row permutation
/// and the sign of the permutation.
pub fn lu_decompose<N: Scalar>(
    lu: &mut [N],
    n: usize,
    perm: &mut [usize],
) -> Result<f64, LinearError> {
    let original: Vec<N> = lu.to_vec();
    for (i, p) in perm.iter_mut().enumerate().take(n) {
        *p = i;
    }
    let mut sign = 1.0;
    for k in 0..n {
        let (piv, max) = (k..n)
            .map(|i| (i, lu[i * n + k].re().abs()))
            .fold((k, -1.0), |acc, x| if x.1 > acc.1 { x } else { acc });
        if max == 0.0 || !max.is_finite() {
            return Err(LinearError::Singular);
        }
        if piv != k {
            for j in 0..n {
                lu.swap(k * n + j, piv * n + j);
            }
            perm.swap(k, piv);
            sign = -sign;
        }
        let pivot = lu[k * n + k];
        for i in k + 1..n {
            let f = lu[i * n + k] / pivot;
            lu[i * n + k] = f;
            for j in k + 1..n {
                let v = lu[k * n + j];
                lu[i * n + j] = lu[i * n + j] - f * v;
            }
        }
    }
    let det: f64 = sign * (0..n).map(|i| lu[i * n + i].re()).product::<f64>();
    if nearly_singular(&original, n, det) {
        return Err(LinearError::Singular);
    }
    Ok(sign)
}

/// Solves `A x = b` given the LU factors from [`lu_decompose`]; `b` is overwritten by `x`.
pub fn lu_solve<N: Scalar>(lu: &[N], n: usize, perm: &[usize], b: &mut [N]) {
    let pb: Vec<N> = perm.iter().map(|&p| b[p]).collect();
    b.copy_from_slice(&pb);
    for i in 0..n {
        let mut s = b[i];
        for j in 0..i {
            s = s - lu[i * n + j] * b[j];
        }
        b[i] = s;
    }
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in i + 1..n {
            s = s - lu[i * n + j] * b[j];
        }
        b[i] = s / lu[i * n + i];
    }
}

fn inverse_lu<N: Scalar>(a: &[N], n: usize, out: &mut [N]) -> Result<(), LinearError> {
    let mut lu = a.to_vec();
    let mut perm = vec![0; n];
    lu_decompose(&mut lu, n, &mut perm)?;
    let mut col = vec![N::zero(); n];
    for j in 0..n {
        for (i, c) in col.iter_mut().enumerate() {
            *c = if i == j { N::one() } else { N::zero() };
        }
        lu_solve(&lu, n, &perm, &mut col);
        for i in 0..n {
            out[i * n + j] = col[i];
        }
    }
    Ok(())
}
