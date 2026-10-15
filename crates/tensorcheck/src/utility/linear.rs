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
/// rows, so a fixed threshold works for every matrix. The ratio is formed from
/// logarithms, so neither the determinant nor the product of the row norms has to
/// be representable: `log_det` is `ln |det A|`.
fn nearly_singular<N: Scalar>(a: &[N], n: usize, log_det: f64) -> bool {
    let mut log_prod = 0.0;
    for i in 0..n {
        let norm = row_norm(a, n, i);
        if norm == 0.0 || !norm.is_finite() {
            return true;
        }
        log_prod += norm.ln();
    }
    let log_ratio = log_det - log_prod;
    // A NaN ratio (from a non-finite entry) is singular too.
    log_ratio.is_nan() || log_ratio < (n as f64 * f64::EPSILON).ln()
}

/// The Euclidean norm of row `i`, scaled by its largest entry so that the squares
/// neither overflow nor underflow.
fn row_norm<N: Scalar>(a: &[N], n: usize, i: usize) -> f64 {
    let big = row_max(a, n, i);
    if big == 0.0 || !big.is_finite() {
        return big;
    }
    big * (0..n)
        .map(|j| (a[i * n + j].re() / big).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// The largest `|a_ij|` of row `i` (NaN if the row has one).
fn row_max<N: Scalar>(a: &[N], n: usize, i: usize) -> f64 {
    (0..n).fold(0.0, |m: f64, j| {
        let x = a[i * n + j].re().abs();
        if x.is_nan() || m.is_nan() {
            f64::NAN
        } else {
            m.max(x)
        }
    })
}

/// Fills `d` with `2^-e`, `2^e` the power of two nearest the largest entry of each
/// row, and `b` with `D A`, whose entries are then of order one. Multiplying by a
/// power of two is exact, so `A^{-1} = (D A)^{-1} D` loses nothing to the scaling.
/// Fails if a row is zero or not finite (DD-08).
fn scale_rows<N: Scalar>(a: &[N], n: usize, d: &mut [f64], b: &mut [N]) -> Result<(), LinearError> {
    for i in 0..n {
        let big = row_max(a, n, i);
        if big == 0.0 || !big.is_finite() {
            return Err(LinearError::Singular);
        }
        // A subnormal row cannot reach order one: 2^1024 is not representable.
        d[i] = pow2((-(big.log2().round() as i32)).min(f64::MAX_EXP - 1));
        for j in 0..n {
            b[i * n + j] = a[i * n + j].scale(d[i]);
        }
    }
    Ok(())
}

/// `2^k` for `-1074 <= k <= 1023`, exactly (subnormal powers included).
fn pow2(k: i32) -> f64 {
    // Split the exponent so that each half is a normal number.
    let half = k / 2;
    2f64.powi(half) * 2f64.powi(k - half)
}

#[inline(always)]
fn det2<N: Scalar>(a: N, b: N, c: N, d: N) -> N {
    a * d - b * c
}

/// Writes the inverse of the `n × n` matrix `a` into `out`.
///
/// The rows are first scaled by powers of two (exactly) to entries of order one, so
/// that a well-conditioned matrix is inverted whatever its scale (DD-08).
pub fn inverse<N: Scalar>(a: &[N], n: usize, out: &mut [N]) -> Result<(), LinearError> {
    debug_assert_eq!(a.len(), n * n);
    debug_assert_eq!(out.len(), n * n);
    if n == 0 {
        return Ok(());
    }
    // Small sizes (every space of the library) stay on the stack.
    let (mut d_small, mut b_small) = ([0.0; 4], [N::zero(); 16]);
    let (mut d_big, mut b_big) = (Vec::new(), Vec::new());
    let (d, b) = if n <= 4 {
        (&mut d_small[..n], &mut b_small[..n * n])
    } else {
        d_big.resize(n, 0.0);
        b_big.resize(n * n, N::zero());
        (&mut d_big[..], &mut b_big[..])
    };
    scale_rows(a, n, d, b)?;
    inverse_unscaled(b, n, out)?;
    // A^{-1} = B^{-1} D: column j of B^{-1} times d_j.
    for (k, o) in out.iter_mut().enumerate() {
        *o = o.scale(d[k % n]);
    }
    Ok(())
}

fn inverse_unscaled<N: Scalar>(a: &[N], n: usize, out: &mut [N]) -> Result<(), LinearError> {
    match n {
        1 => {
            if a[0].re() == 0.0 || !a[0].re().is_finite() {
                return Err(LinearError::Singular);
            }
            out[0] = N::one() / a[0];
            Ok(())
        }
        2 => {
            let det = det2(a[0], a[1], a[2], a[3]);
            if nearly_singular(a, 2, det.re().abs().ln()) {
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
            if nearly_singular(a, 3, det.re().abs().ln()) {
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
    if nearly_singular(a, 4, det.re().abs().ln()) {
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

/// LU decomposition with partial pivoting, in place, rejecting matrices that are
/// singular to working precision (Hadamard ratio below `n ε`, DD-08). Returns the
/// sign of the permutation.
pub fn lu_decompose<N: Scalar>(
    lu: &mut [N],
    n: usize,
    perm: &mut [usize],
) -> Result<f64, LinearError> {
    lu_decompose_with(lu, n, perm, true)
}

/// LU decomposition with partial pivoting. With `check_conditioning = false` only an
/// exactly zero or non-finite pivot is an error: the stiff solver's matrix
/// `I/(hγ) − J` becomes ill-conditioned for large steps when `J` has a zero
/// eigenvalue (a conserved quantity), yet the solves remain accurate enough there.
pub fn lu_decompose_with<N: Scalar>(
    lu: &mut [N],
    n: usize,
    perm: &mut [usize],
    check_conditioning: bool,
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
    let log_det: f64 = (0..n).map(|i| lu[i * n + i].re().abs().ln()).sum();
    if check_conditioning && nearly_singular(&original, n, log_det) {
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
