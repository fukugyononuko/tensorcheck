//! Fixed-size array types for a space's dimension (used by A01 and the dual numbers).

/// The memory side of a space: its number of dimensions and its array type.
///
/// `Space` extends this trait. It is a utility (a family of fixed-size array
/// types) so that the component layout and the dual numbers can both use it.
pub trait Extent: 'static {
    /// Number of dimensions `n`.
    const N: usize;
    /// `[T; n]`.
    type Arr<T: Copy>: Copy + AsRef<[T]> + AsMut<[T]>;
    /// `[T; n(n+1)/2]`, the independent components of a symmetric rank-2 tensor.
    type SymArr<T: Copy>: Copy + AsRef<[T]> + AsMut<[T]>;
    /// `[T; n(n-1)/2]`, the independent components of an antisymmetric rank-2 tensor.
    type AntiArr<T: Copy>: Copy + AsRef<[T]> + AsMut<[T]>;
    /// Builds an `[T; n]` from a function of the index.
    fn from_fn<T: Copy>(f: impl FnMut(usize) -> T) -> Self::Arr<T>;
    /// Builds an `[T; n(n+1)/2]` from a function of the index.
    fn sym_from_fn<T: Copy>(f: impl FnMut(usize) -> T) -> Self::SymArr<T>;
    /// Builds an `[T; n(n-1)/2]` from a function of the index.
    fn anti_from_fn<T: Copy>(f: impl FnMut(usize) -> T) -> Self::AntiArr<T>;
    /// Views a slice of `[T; n]` as one flat slice (no copy, no `unsafe`).
    fn flatten<T: Copy>(xs: &[Self::Arr<T>]) -> &[T];
}

/// Largest rank of a tensor (and of the slots of one contraction).
pub const MAX_RANK: usize = 12;

/// `[x, a[0], a[1], …]` (drops the last element), in a constant context.
pub const fn prepend(x: usize, a: [usize; MAX_RANK]) -> [usize; MAX_RANK] {
    let mut out = [0; MAX_RANK];
    out[0] = x;
    let mut i = 1;
    while i < MAX_RANK {
        out[i] = a[i - 1];
        i += 1;
    }
    out
}

/// Position of `(i, j)` (any order) in the independent components of a symmetric
/// `n × n` matrix stored row by row with `i <= j`.
#[inline(always)]
pub fn sym_index(i: usize, j: usize, n: usize) -> usize {
    let (i, j) = if i <= j { (i, j) } else { (j, i) };
    i * n - i * (i + 1) / 2 + j
}

/// Position of `(i, j)` with `i < j` in the independent components of an
/// antisymmetric `n × n` matrix stored row by row (strictly upper triangle).
#[inline(always)]
pub fn anti_index(i: usize, j: usize, n: usize) -> usize {
    debug_assert!(i < j);
    i * n - i * (i + 1) / 2 + (j - i - 1)
}
