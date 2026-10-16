//! Type-level naturals, integers, booleans and lists (CMP-U01, DD-02).
//!
//! Everything here is computed by the trait solver at compile time. The types are
//! zero-sized markers and never exist at run time.

use core::marker::PhantomData;

// ---------------------------------------------------------------------------
// Booleans

/// Type-level `true`.
pub struct True;
/// Type-level `false`.
pub struct False;

/// A type-level boolean.
pub trait Bool {
    /// The value as a `bool`.
    const VALUE: bool;
}
impl Bool for True {
    const VALUE: bool = true;
}
impl Bool for False {
    const VALUE: bool = false;
}

// ---------------------------------------------------------------------------
// Naturals (Peano numbers)

/// Type-level zero.
pub struct Z;
/// Type-level successor of `N`.
pub struct S<N>(PhantomData<N>);

/// A type-level natural number.
pub trait Nat: 'static {
    /// The value as a `usize`.
    const VALUE: usize;
}
impl Nat for Z {
    const VALUE: usize = 0;
}
impl<N: Nat> Nat for S<N> {
    const VALUE: usize = N::VALUE + 1;
}

/// Equality of naturals: `Out` is [`True`] or [`False`].
pub trait NatEq<B> {
    /// The result.
    type Out: Bool;
}
impl NatEq<Z> for Z {
    type Out = True;
}
impl<N> NatEq<S<N>> for Z {
    type Out = False;
}
impl<N> NatEq<Z> for S<N> {
    type Out = False;
}
impl<N: NatEq<M>, M> NatEq<S<M>> for S<N> {
    type Out = <N as NatEq<M>>::Out;
}

/// Addition of naturals.
pub trait NatAdd<B> {
    /// The sum.
    type Out;
}
impl<B> NatAdd<B> for Z {
    type Out = B;
}
impl<N: NatAdd<B>, B> NatAdd<B> for S<N> {
    type Out = S<<N as NatAdd<B>>::Out>;
}

/// `Self - B` as a signed integer.
pub trait NatSub<B> {
    /// The difference.
    type Out;
}
impl NatSub<Z> for Z {
    type Out = Zero;
}
impl<N> NatSub<Z> for S<N> {
    type Out = Pos<S<N>>;
}
impl<N> NatSub<S<N>> for Z {
    type Out = Neg<S<N>>;
}
impl<N: NatSub<M>, M> NatSub<S<M>> for S<N> {
    type Out = <N as NatSub<M>>::Out;
}

// ---------------------------------------------------------------------------
// Integers. Zero has exactly one representation, so equal integers are equal types.

/// Type-level integer zero.
pub struct Zero;
/// Type-level positive integer `N` (`N >= 1`).
pub struct Pos<N>(PhantomData<N>);
/// Type-level negative integer `-N` (`N >= 1`).
pub struct Neg<N>(PhantomData<N>);

/// A type-level integer.
pub trait Int: 'static {
    /// The value as an `i32`.
    const VALUE: i32;
}
impl Int for Zero {
    const VALUE: i32 = 0;
}
impl<N: Nat> Int for Pos<S<N>> {
    const VALUE: i32 = <S<N> as Nat>::VALUE as i32;
}
impl<N: Nat> Int for Neg<S<N>> {
    const VALUE: i32 = -(<S<N> as Nat>::VALUE as i32);
}

/// Integer addition.
pub trait IntAdd<B> {
    /// The sum.
    type Out;
}
impl<B> IntAdd<B> for Zero {
    type Out = B;
}
impl<A> IntAdd<Zero> for Pos<A> {
    type Out = Pos<A>;
}
impl<A> IntAdd<Zero> for Neg<A> {
    type Out = Neg<A>;
}
impl<A: NatAdd<B>, B> IntAdd<Pos<B>> for Pos<A> {
    type Out = Pos<<A as NatAdd<B>>::Out>;
}
impl<A: NatAdd<B>, B> IntAdd<Neg<B>> for Neg<A> {
    type Out = Neg<<A as NatAdd<B>>::Out>;
}
impl<A: NatSub<B>, B> IntAdd<Neg<B>> for Pos<A> {
    type Out = <A as NatSub<B>>::Out;
}
impl<A, B: NatSub<A>> IntAdd<Pos<B>> for Neg<A> {
    type Out = <B as NatSub<A>>::Out;
}

/// Integer equality: `Out` is [`True`] or [`False`].
pub trait IntEq<B> {
    /// The result.
    type Out: Bool;
}
impl IntEq<Zero> for Zero {
    type Out = True;
}
impl<B> IntEq<Pos<B>> for Zero {
    type Out = False;
}
impl<B> IntEq<Neg<B>> for Zero {
    type Out = False;
}
impl<A> IntEq<Zero> for Pos<A> {
    type Out = False;
}
impl<A> IntEq<Zero> for Neg<A> {
    type Out = False;
}
impl<A, B> IntEq<Neg<B>> for Pos<A> {
    type Out = False;
}
impl<A, B> IntEq<Pos<B>> for Neg<A> {
    type Out = False;
}
impl<A: NatEq<B>, B> IntEq<Pos<B>> for Pos<A> {
    type Out = <A as NatEq<B>>::Out;
}
impl<A: NatEq<B>, B> IntEq<Neg<B>> for Neg<A> {
    type Out = <A as NatEq<B>>::Out;
}

/// Logical and.
pub trait And<B> {
    /// The result.
    type Out: Bool;
}
impl<B: Bool> And<B> for True {
    type Out = B;
}
impl<B> And<B> for False {
    type Out = False;
}

/// Logical or.
pub trait Or<B> {
    /// The result.
    type Out: Bool;
}
impl<B> Or<B> for True {
    type Out = True;
}
impl<B: Bool> Or<B> for False {
    type Out = B;
}

/// Element-wise equality of two lists of integers of equal length.
pub trait IntListEq<B> {
    /// The result.
    type Out: Bool;
}
impl IntListEq<Nil> for Nil {
    type Out = True;
}
impl<A, B, TA, TB> IntListEq<Cons<B, TB>> for Cons<A, TA>
where
    A: IntEq<B>,
    TA: IntListEq<TB>,
    <A as IntEq<B>>::Out: And<<TA as IntListEq<TB>>::Out>,
{
    type Out = <<A as IntEq<B>>::Out as And<<TA as IntListEq<TB>>::Out>>::Out;
}

/// Integer negation.
pub trait IntNeg {
    /// `-Self`.
    type Out;
}
impl IntNeg for Zero {
    type Out = Zero;
}
impl<A> IntNeg for Pos<A> {
    type Out = Neg<A>;
}
impl<A> IntNeg for Neg<A> {
    type Out = Pos<A>;
}

/// Integer subtraction.
pub trait IntSub<B> {
    /// `Self - B`.
    type Out;
}
impl<A, B: IntNeg> IntSub<B> for A
where
    A: IntAdd<<B as IntNeg>::Out>,
{
    type Out = <A as IntAdd<<B as IntNeg>::Out>>::Out;
}

/// Multiplication of an integer by a type-level integer `K` (used for powers of dimensions).
pub trait IntMul<K> {
    /// `Self * K`.
    type Out;
}
impl<A, K: ScaleInt<A>> IntMul<K> for A {
    type Out = <K as ScaleInt<A>>::Out;
}

/// `A * Self`: the recursion of [`IntMul`] on the factor. The factor is the `Self` type
/// so that, while it is still unknown (a method's `const` argument before it is
/// substituted), the compiler stops with "ambiguous" instead of recursing without end.
pub trait ScaleInt<A> {
    /// `A * Self`.
    type Out;
}
impl<A> ScaleInt<A> for Zero {
    type Out = Zero;
}
impl<A> ScaleInt<A> for Pos<S<Z>> {
    type Out = A;
}
impl<A, N> ScaleInt<A> for Pos<S<S<N>>>
where
    Pos<S<N>>: ScaleInt<A>,
    <Pos<S<N>> as ScaleInt<A>>::Out: IntAdd<A>,
{
    type Out = <<Pos<S<N>> as ScaleInt<A>>::Out as IntAdd<A>>::Out;
}
impl<A, N> ScaleInt<A> for Neg<N>
where
    Pos<N>: ScaleInt<A>,
    <Pos<N> as ScaleInt<A>>::Out: IntNeg,
{
    type Out = <<Pos<N> as ScaleInt<A>>::Out as IntNeg>::Out;
}

/// Maps a `const` integer to its type-level integer, for `-20..=20`.
pub struct ConstInt<const K: i32>;

/// Implemented by [`ConstInt<K>`] for supported `K`.
pub trait ToInt {
    /// The type-level integer.
    type Out: Int;
}

macro_rules! int_aliases {
    ($($p:ident $n:ident $v:literal = $nat:ty;)*) => {
        $(
            #[doc = concat!("Type-level integer `", stringify!($v), "`.")]
            pub type $p = Pos<$nat>;
            #[doc = concat!("Type-level integer `-", stringify!($v), "`.")]
            pub type $n = Neg<$nat>;
            impl ToInt for ConstInt<$v> { type Out = $p; }
            impl ToInt for ConstInt<{ -$v }> { type Out = $n; }
        )*
    };
}
/// Type-level integer `0`.
pub type Z0 = Zero;
impl ToInt for ConstInt<0> {
    type Out = Zero;
}
int_aliases! {
    P1 N1 1 = S<Z>;
    P2 N2 2 = S<S<Z>>;
    P3 N3 3 = S<S<S<Z>>>;
    P4 N4 4 = S<S<S<S<Z>>>>;
    P5 N5 5 = S<S<S<S<S<Z>>>>>;
    P6 N6 6 = S<S<S<S<S<S<Z>>>>>>;
    P7 N7 7 = S<S<S<S<S<S<S<Z>>>>>>>;
    P8 N8 8 = S<S<S<S<S<S<S<S<Z>>>>>>>>;
    P9 N9 9 = S<S<S<S<S<S<S<S<S<Z>>>>>>>>>;
    P10 N10 10 = S<S<S<S<S<S<S<S<S<S<Z>>>>>>>>>>;
    P11 N11 11 = S<S<S<S<S<S<S<S<S<S<S<Z>>>>>>>>>>>;
    P12 N12 12 = S<S<S<S<S<S<S<S<S<S<S<S<Z>>>>>>>>>>>>;
    P13 N13 13 = S<S<S<S<S<S<S<S<S<S<S<S<S<Z>>>>>>>>>>>>>;
    P14 N14 14 = S<S<S<S<S<S<S<S<S<S<S<S<S<S<Z>>>>>>>>>>>>>>;
    P15 N15 15 = S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<Z>>>>>>>>>>>>>>>;
    P16 N16 16 = S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<Z>>>>>>>>>>>>>>>>;
    P17 N17 17 = S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<Z>>>>>>>>>>>>>>>>>;
    P18 N18 18 = S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<Z>>>>>>>>>>>>>>>>>>;
    P19 N19 19 = S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<Z>>>>>>>>>>>>>>>>>>>;
    P20 N20 20 = S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<S<Z>>>>>>>>>>>>>>>>>>>>;
}

// ---------------------------------------------------------------------------
// Lists

/// The empty type-level list.
pub struct Nil;
/// A type-level list with head `H` and tail `T`.
pub struct Cons<H, T>(PhantomData<(H, T)>);

/// The length of a type-level list.
pub trait Len {
    /// The length.
    const LEN: usize;
}
impl Len for Nil {
    const LEN: usize = 0;
}
impl<H, T: Len> Len for Cons<H, T> {
    const LEN: usize = T::LEN + 1;
}

/// Concatenation of type-level lists.
pub trait Concat<B> {
    /// `Self ++ B`.
    type Out;
}
impl<B> Concat<B> for Nil {
    type Out = B;
}
impl<H, T: Concat<B>, B> Concat<B> for Cons<H, T> {
    type Out = Cons<H, <T as Concat<B>>::Out>;
}

/// Element-wise integer addition of two lists of equal length.
pub trait ZipAdd<B> {
    /// The element-wise sum.
    type Out;
}
impl ZipAdd<Nil> for Nil {
    type Out = Nil;
}
impl<A: IntAdd<B>, B, TA: ZipAdd<TB>, TB> ZipAdd<Cons<B, TB>> for Cons<A, TA> {
    type Out = Cons<<A as IntAdd<B>>::Out, <TA as ZipAdd<TB>>::Out>;
}

/// Element-wise integer subtraction of two lists of equal length.
pub trait ZipSub<B> {
    /// The element-wise difference.
    type Out;
}
impl ZipSub<Nil> for Nil {
    type Out = Nil;
}
impl<A: IntSub<B>, B, TA: ZipSub<TB>, TB> ZipSub<Cons<B, TB>> for Cons<A, TA> {
    type Out = Cons<<A as IntSub<B>>::Out, <TA as ZipSub<TB>>::Out>;
}

/// Multiplies every integer in a list by `K`.
pub trait MapMul<K> {
    /// The scaled list.
    type Out;
}
impl<K> MapMul<K> for Nil {
    type Out = Nil;
}
impl<A: IntMul<K>, T: MapMul<K>, K> MapMul<K> for Cons<A, T> {
    type Out = Cons<<A as IntMul<K>>::Out, <T as MapMul<K>>::Out>;
}

/// Half of a natural, rounded down. Total, so that a square root can name its result
/// type before the evenness check ([`IsEven`]) reports an odd exponent.
pub trait NatHalf {
    /// `Self / 2`, rounded down.
    type Out;
}
impl NatHalf for Z {
    type Out = Z;
}
impl NatHalf for S<Z> {
    type Out = Z;
}
impl<N: NatHalf> NatHalf for S<S<N>> {
    type Out = S<<N as NatHalf>::Out>;
}

/// Whether a natural or an integer is even: `Out` is [`True`] or [`False`].
pub trait IsEven {
    /// The result.
    type Out: Bool;
}
impl IsEven for Z {
    type Out = True;
}
impl IsEven for S<Z> {
    type Out = False;
}
impl<N: IsEven> IsEven for S<S<N>> {
    type Out = <N as IsEven>::Out;
}
impl IsEven for Zero {
    type Out = True;
}
impl<N: IsEven> IsEven for Pos<N> {
    type Out = <N as IsEven>::Out;
}
impl<N: IsEven> IsEven for Neg<N> {
    type Out = <N as IsEven>::Out;
}

/// Half of an integer, rounded toward zero (used for square roots of dimensions).
pub trait IntHalf {
    /// `Self / 2`, rounded toward zero.
    type Out;
}
impl IntHalf for Zero {
    type Out = Zero;
}
impl IntHalf for Pos<S<Z>> {
    type Out = Zero;
}
impl<N: NatHalf> IntHalf for Pos<S<S<N>>> {
    type Out = Pos<S<<N as NatHalf>::Out>>;
}
impl IntHalf for Neg<S<Z>> {
    type Out = Zero;
}
impl<N: NatHalf> IntHalf for Neg<S<S<N>>> {
    type Out = Neg<S<<N as NatHalf>::Out>>;
}

/// Halves every integer in a list.
pub trait MapHalf {
    /// The halved list.
    type Out;
}
impl MapHalf for Nil {
    type Out = Nil;
}
impl<A: IntHalf, T: MapHalf> MapHalf for Cons<A, T> {
    type Out = Cons<<A as IntHalf>::Out, <T as MapHalf>::Out>;
}

/// Whether every integer in a list is even.
pub trait AllEven {
    /// The result.
    type Out: Bool;
}
impl AllEven for Nil {
    type Out = True;
}
impl<A: IsEven, T: AllEven> AllEven for Cons<A, T>
where
    <A as IsEven>::Out: And<<T as AllEven>::Out>,
{
    type Out = <<A as IsEven>::Out as And<<T as AllEven>::Out>>::Out;
}

/// Reads a list of type-level integers into a slice.
pub trait IntList {
    /// Writes the values into `out` (which must have the list's length).
    fn write(out: &mut [i32]);
}
impl IntList for Nil {
    fn write(_: &mut [i32]) {}
}
impl<A: Int, T: IntList> IntList for Cons<A, T> {
    fn write(out: &mut [i32]) {
        out[0] = A::VALUE;
        T::write(&mut out[1..]);
    }
}

/// Maps a `const` natural number to its type-level natural, for `0..64`.
pub struct ConstNat<const K: usize>;

/// Implemented by [`ConstNat<K>`] for supported `K`.
pub trait ToNat {
    /// The type-level natural.
    type Out: Nat + NatEq<Self::Out>;
}

macro_rules! const_nats {
    ($prev:ty; $k:literal $(, $rest:literal)*) => {
        impl ToNat for ConstNat<$k> { type Out = S<$prev>; }
        const_nats!(S<$prev>; $($rest),*);
    };
    ($prev:ty;) => {};
}
impl ToNat for ConstNat<0> {
    type Out = Z;
}
const_nats!(Z; 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
    51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63);
