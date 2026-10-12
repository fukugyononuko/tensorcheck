//! How tensor components are laid out in memory (CMP-A01, DD-04).
//!
//! Each index slot contributes one level of nesting, using the array type of the
//! slot's space (for example `[T; 4]` for spacetime). Stable Rust cannot use a
//! type-level product as an array length, but composing the per-space array types
//! needs no arithmetic at all.

use crate::utility::typelevel::{Cons, Nil};

/// The memory side of a space: its number of dimensions and its array type.
///
/// `Space` (in the index engine) extends this trait. It lives here so that the
/// component layout does not depend on the engine layer.
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
}

/// An index slot knows the extent of the space it belongs to.
pub trait SlotExtent {
    /// The extent.
    type Ext: Extent;
}

/// A list of index slots and the nested-array storage it implies.
pub trait Slots: 'static {
    /// The nested array holding components of type `T`.
    type Store<T: Copy>: Copy;
    /// Rank (number of slots).
    const RANK: usize;
    /// Number of components, `∏ n_i`.
    const LEN: usize;
    /// Reads the component with row-major flat index `flat`.
    fn get<T: Copy>(s: &Self::Store<T>, flat: usize) -> T;
    /// Writes the component with row-major flat index `flat`.
    fn set<T: Copy>(s: &mut Self::Store<T>, flat: usize, v: T);
    /// Builds a store from a function of the flat index.
    fn from_fn<T: Copy>(f: impl FnMut(usize) -> T) -> Self::Store<T>;
    /// The extent `n` of slot number `slot`.
    fn extent(slot: usize) -> usize;
}

impl Slots for Nil {
    type Store<T: Copy> = T;
    const RANK: usize = 0;
    const LEN: usize = 1;
    #[inline(always)]
    fn get<T: Copy>(s: &T, _: usize) -> T {
        *s
    }
    #[inline(always)]
    fn set<T: Copy>(s: &mut T, _: usize, v: T) {
        *s = v
    }
    #[inline(always)]
    fn from_fn<T: Copy>(mut f: impl FnMut(usize) -> T) -> T {
        f(0)
    }
    fn extent(_: usize) -> usize {
        panic!("a rank-0 tensor has no slots")
    }
}

impl<H: SlotExtent + 'static, Tl: Slots> Slots for Cons<H, Tl> {
    type Store<T: Copy> = <H::Ext as Extent>::Arr<Tl::Store<T>>;
    const RANK: usize = Tl::RANK + 1;
    const LEN: usize = <H::Ext as Extent>::N * Tl::LEN;
    #[inline(always)]
    fn get<T: Copy>(s: &Self::Store<T>, flat: usize) -> T {
        Tl::get(&s.as_ref()[flat / Tl::LEN], flat % Tl::LEN)
    }
    #[inline(always)]
    fn set<T: Copy>(s: &mut Self::Store<T>, flat: usize, v: T) {
        Tl::set(&mut s.as_mut()[flat / Tl::LEN], flat % Tl::LEN, v)
    }
    #[inline(always)]
    fn from_fn<T: Copy>(mut f: impl FnMut(usize) -> T) -> Self::Store<T> {
        <H::Ext as Extent>::from_fn(|i| Tl::from_fn(|j| f(i * Tl::LEN + j)))
    }
    fn extent(slot: usize) -> usize {
        if slot == 0 {
            <H::Ext as Extent>::N
        } else {
            Tl::extent(slot - 1)
        }
    }
}

/// Builds a store with every component equal to `v`.
#[inline(always)]
pub fn fill<S: Slots, T: Copy>(v: T) -> S::Store<T> {
    S::from_fn(|_| v)
}
