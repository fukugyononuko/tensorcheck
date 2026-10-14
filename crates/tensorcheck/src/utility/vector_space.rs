//! States of numerical solvers (CMP-U05, DD-11).
//!
//! A state is anything whose components can be read and written as numbers in base
//! units: tensors, points, and tuples of them. Solvers keep their work arrays as
//! plain `f64` slices and convert at the boundaries, so engines never see
//! dimensions while the user-facing types keep them.

/// A solver state (API-12).
pub trait VectorSpace: Copy {
    /// Number of components.
    const LEN: usize;
    /// Component number `i`, in base units.
    fn component(&self, i: usize) -> f64;
    /// A state from its components in base units.
    fn from_components(f: impl FnMut(usize) -> f64) -> Self;
}

/// The type of `self / [D]`: the time derivative of a state has type
/// `<Y as PerDim<Time>>::Out`.
pub trait PerDim<D> {
    /// The resulting state type.
    type Out: VectorSpace;
}

/// The type of `self · [D]`: the integral of a state over `x` has type
/// `<Y as TimesDim<X>>::Out`.
pub trait TimesDim<D> {
    /// The resulting state type.
    type Out: VectorSpace;
}

macro_rules! tuple_state {
    ($($t:ident $i:tt),+) => {
        impl<$($t: VectorSpace),+> VectorSpace for ($($t,)+) {
            const LEN: usize = 0 $(+ $t::LEN)+;
            fn component(&self, mut i: usize) -> f64 {
                $(
                    if i < $t::LEN {
                        return self.$i.component(i);
                    }
                    i -= $t::LEN;
                )+
                let _ = i;
                panic!("component index out of range")
            }
            #[allow(unused_assignments)]
            fn from_components(mut f: impl FnMut(usize) -> f64) -> Self {
                let mut offset = 0;
                ($(
                    {
                        let o = offset;
                        offset += $t::LEN;
                        $t::from_components(|k| f(o + k))
                    },
                )+)
            }
        }
        impl<D, $($t: PerDim<D>),+> PerDim<D> for ($($t,)+) {
            type Out = ($(<$t as PerDim<D>>::Out,)+);
        }
        impl<D, $($t: TimesDim<D>),+> TimesDim<D> for ($($t,)+) {
            type Out = ($(<$t as TimesDim<D>>::Out,)+);
        }
    };
}
tuple_state!(A 0, B 1);
tuple_state!(A 0, B 1, C 2);
tuple_state!(A 0, B 1, C 2, E 3);
tuple_state!(A 0, B 1, C 2, E 3, F 4);
tuple_state!(A 0, B 1, C 2, E 3, F 4, G 5);
