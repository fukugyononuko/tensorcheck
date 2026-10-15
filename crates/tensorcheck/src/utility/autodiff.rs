//! Dual numbers for automatic differentiation (CMP-U02, DD-09).
//!
//! `Dual<E>` carries a value and its first derivatives with respect to the `n`
//! coordinates of a space with extent `E`; `Dual2<E>` also carries the second
//! derivatives. Every operation applies the chain rule, so a function written
//! generically over `N: Scalar` yields exact derivatives (up to rounding).

// The product and quotient rules add terms inside `Mul` and `Div` impls.
#![allow(clippy::suspicious_arithmetic_impl)]

use core::fmt;
use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::utility::extent::{Extent, sym_index};
use crate::utility::scalar::Scalar;

/// A value with its first derivatives `∂_c v` (API-17).
pub struct Dual<E: Extent> {
    /// The value.
    pub re: f64,
    /// The first derivatives.
    pub d: E::Arr<f64>,
}

/// A value with its first and second derivatives (API-17). The second derivatives
/// `∂_c ∂_d v` are stored for `c <= d`, row by row.
pub struct Dual2<E: Extent> {
    /// The value.
    pub re: f64,
    /// The first derivatives.
    pub d: E::Arr<f64>,
    /// The second derivatives, independent components.
    pub h: E::SymArr<f64>,
}

impl<E: Extent> Clone for Dual<E> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<E: Extent> Copy for Dual<E> {}
impl<E: Extent> Clone for Dual2<E> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<E: Extent> Copy for Dual2<E> {}

impl<E: Extent> fmt::Debug for Dual<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Dual({}, {:?})", self.re, self.d.as_ref())
    }
}
impl<E: Extent> fmt::Debug for Dual2<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Dual2({}, {:?}, {:?})",
            self.re,
            self.d.as_ref(),
            self.h.as_ref()
        )
    }
}
impl<E: Extent> PartialEq for Dual<E> {
    fn eq(&self, o: &Self) -> bool {
        self.re == o.re && self.d.as_ref() == o.d.as_ref()
    }
}
impl<E: Extent> PartialEq for Dual2<E> {
    fn eq(&self, o: &Self) -> bool {
        self.re == o.re && self.d.as_ref() == o.d.as_ref() && self.h.as_ref() == o.h.as_ref()
    }
}

impl<E: Extent> Dual<E> {
    /// The independent variable number `i` with value `x` (derivative `e_i`).
    pub fn variable(x: f64, i: usize) -> Self {
        Dual {
            re: x,
            d: E::from_fn(|k| if k == i { 1.0 } else { 0.0 }),
        }
    }
    fn constant(x: f64) -> Self {
        Dual {
            re: x,
            d: E::from_fn(|_| 0.0),
        }
    }
    /// Applies `f` with `f(v) = f0` and `f'(v) = f1`.
    fn chain(self, f0: f64, f1: f64) -> Self {
        let d = self.d;
        Dual {
            re: f0,
            d: E::from_fn(|k| f1 * d.as_ref()[k]),
        }
    }
}

impl<E: Extent> Dual2<E> {
    /// The independent variable number `i` with value `x`.
    pub fn variable(x: f64, i: usize) -> Self {
        Dual2 {
            re: x,
            d: E::from_fn(|k| if k == i { 1.0 } else { 0.0 }),
            h: E::sym_from_fn(|_| 0.0),
        }
    }
    fn constant(x: f64) -> Self {
        Dual2 {
            re: x,
            d: E::from_fn(|_| 0.0),
            h: E::sym_from_fn(|_| 0.0),
        }
    }
    /// The second derivative `∂_c ∂_d`.
    pub fn second(&self, c: usize, d: usize) -> f64 {
        self.h.as_ref()[sym_index(c, d, E::N)]
    }
    /// Applies `f` with `f(v) = f0`, `f'(v) = f1`, `f''(v) = f2`:
    /// `∂_c∂_d f = f2 ∂_c v ∂_d v + f1 ∂_c∂_d v`.
    fn chain(self, f0: f64, f1: f64, f2: f64) -> Self {
        let n = E::N;
        let d = self.d;
        let h = self.h;
        let mut out_h = E::sym_from_fn(|_| 0.0);
        for c in 0..n {
            for e in c..n {
                let k = sym_index(c, e, n);
                out_h.as_mut()[k] = f2 * d.as_ref()[c] * d.as_ref()[e] + f1 * h.as_ref()[k];
            }
        }
        Dual2 {
            re: f0,
            d: E::from_fn(|k| f1 * d.as_ref()[k]),
            h: out_h,
        }
    }
}

// ---------------------------------------------------------------------------
// Arithmetic

impl<E: Extent> Add for Dual<E> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Dual {
            re: self.re + o.re,
            d: E::from_fn(|k| self.d.as_ref()[k] + o.d.as_ref()[k]),
        }
    }
}
impl<E: Extent> Sub for Dual<E> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Dual {
            re: self.re - o.re,
            d: E::from_fn(|k| self.d.as_ref()[k] - o.d.as_ref()[k]),
        }
    }
}
impl<E: Extent> Neg for Dual<E> {
    type Output = Self;
    fn neg(self) -> Self {
        Dual {
            re: -self.re,
            d: E::from_fn(|k| -self.d.as_ref()[k]),
        }
    }
}
impl<E: Extent> Mul for Dual<E> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Dual {
            re: self.re * o.re,
            d: E::from_fn(|k| self.d.as_ref()[k] * o.re + self.re * o.d.as_ref()[k]),
        }
    }
}
impl<E: Extent> Div for Dual<E> {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let inv = 1.0 / o.re;
        self * o.chain(inv, -inv * inv)
    }
}

impl<E: Extent> Add for Dual2<E> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Dual2 {
            re: self.re + o.re,
            d: E::from_fn(|k| self.d.as_ref()[k] + o.d.as_ref()[k]),
            h: E::sym_from_fn(|k| self.h.as_ref()[k] + o.h.as_ref()[k]),
        }
    }
}
impl<E: Extent> Sub for Dual2<E> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Dual2 {
            re: self.re - o.re,
            d: E::from_fn(|k| self.d.as_ref()[k] - o.d.as_ref()[k]),
            h: E::sym_from_fn(|k| self.h.as_ref()[k] - o.h.as_ref()[k]),
        }
    }
}
impl<E: Extent> Neg for Dual2<E> {
    type Output = Self;
    fn neg(self) -> Self {
        Dual2 {
            re: -self.re,
            d: E::from_fn(|k| -self.d.as_ref()[k]),
            h: E::sym_from_fn(|k| -self.h.as_ref()[k]),
        }
    }
}
impl<E: Extent> Mul for Dual2<E> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        let n = E::N;
        let (a, b) = (self, o);
        let mut h = E::sym_from_fn(|_| 0.0);
        for c in 0..n {
            for e in c..n {
                let k = sym_index(c, e, n);
                h.as_mut()[k] = a.h.as_ref()[k] * b.re
                    + a.re * b.h.as_ref()[k]
                    + a.d.as_ref()[c] * b.d.as_ref()[e]
                    + a.d.as_ref()[e] * b.d.as_ref()[c];
            }
        }
        Dual2 {
            re: a.re * b.re,
            d: E::from_fn(|k| a.d.as_ref()[k] * b.re + a.re * b.d.as_ref()[k]),
            h,
        }
    }
}
impl<E: Extent> Div for Dual2<E> {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let inv = 1.0 / o.re;
        self * o.chain(inv, -inv * inv, 2.0 * inv * inv * inv)
    }
}

// ---------------------------------------------------------------------------
// Elementary functions: f(v), f'(v) (and f''(v) for Dual2).

macro_rules! scalar_impl {
    ($t:ident, $chain:ident) => {
        impl<E: Extent> Scalar for $t<E> {
            fn from_f64(x: f64) -> Self {
                Self::constant(x)
            }
            fn re(self) -> f64 {
                self.re
            }
            fn sqrt(self) -> Self {
                let s = self.re.sqrt();
                $chain!(self, s, 0.5 / s, -0.25 / (s * self.re))
            }
            fn exp(self) -> Self {
                let e = self.re.exp();
                $chain!(self, e, e, e)
            }
            fn ln(self) -> Self {
                let x = self.re;
                $chain!(self, x.ln(), 1.0 / x, -1.0 / (x * x))
            }
            fn sin(self) -> Self {
                let (s, c) = self.re.sin_cos();
                $chain!(self, s, c, -s)
            }
            fn cos(self) -> Self {
                let (s, c) = self.re.sin_cos();
                $chain!(self, c, -s, -c)
            }
            fn tan(self) -> Self {
                let t = self.re.tan();
                let sec2 = 1.0 + t * t;
                $chain!(self, t, sec2, 2.0 * t * sec2)
            }
            fn sinh(self) -> Self {
                let (s, c) = (self.re.sinh(), self.re.cosh());
                $chain!(self, s, c, s)
            }
            fn cosh(self) -> Self {
                let (s, c) = (self.re.sinh(), self.re.cosh());
                $chain!(self, c, s, c)
            }
            fn tanh(self) -> Self {
                let t = self.re.tanh();
                // 1/cosh² instead of 1 − tanh² (which cancels badly for large |x|).
                let ch = self.re.cosh();
                let sech2 = 1.0 / (ch * ch);
                $chain!(self, t, sech2, -2.0 * t * sech2)
            }
            fn asin(self) -> Self {
                let x = self.re;
                let q = 1.0 - x * x;
                $chain!(self, x.asin(), 1.0 / q.sqrt(), x / (q * q.sqrt()))
            }
            fn acos(self) -> Self {
                let x = self.re;
                let q = 1.0 - x * x;
                $chain!(self, x.acos(), -1.0 / q.sqrt(), -x / (q * q.sqrt()))
            }
            fn atan(self) -> Self {
                let x = self.re;
                let q = 1.0 + x * x;
                $chain!(self, x.atan(), 1.0 / q, -2.0 * x / (q * q))
            }
            fn atan2(self, x: Self) -> Self {
                // Away from the origin, atan2(y, x) has the derivatives of atan(y/x)
                // (or of π/2 − atan(x/y)); only the value needs the quadrant.
                let value = self.re.atan2(x.re);
                let t = if x.re.abs() >= self.re.abs() {
                    (self / x).atan()
                } else {
                    -(x / self).atan()
                };
                let mut out = t;
                out.re = value;
                out
            }
            fn powi(self, k: i32) -> Self {
                // x^0 and x^1 directly: the general formulas multiply a zero
                // coefficient by 0^-1 = ∞ at x = 0, which gives NaN.
                match k {
                    0 => return Self::constant(1.0),
                    1 => return self,
                    _ => {}
                }
                let x = self.re;
                let kf = k as f64;
                // k − 1 and k − 2 do not fit in i32 near i32::MIN; powf takes the
                // same integral exponent there.
                let pow =
                    |j: i64| i32::try_from(j).map_or_else(|_| x.powf(j as f64), |j| x.powi(j));
                $chain!(
                    self,
                    x.powi(k),
                    kf * pow(i64::from(k) - 1),
                    kf * (kf - 1.0) * pow(i64::from(k) - 2)
                )
            }
            fn powf(self, p: f64) -> Self {
                // As for powi: exponents 0 and 1 would give 0 · ∞ at x = 0.
                if p == 0.0 {
                    return Self::constant(1.0);
                }
                if p == 1.0 {
                    return self;
                }
                let x = self.re;
                $chain!(
                    self,
                    x.powf(p),
                    p * x.powf(p - 1.0),
                    p * (p - 1.0) * x.powf(p - 2.0)
                )
            }
            fn abs(self) -> Self {
                if self.re < 0.0 { -self } else { self }
            }
        }
    };
}
macro_rules! chain1 {
    ($v:expr, $f0:expr, $f1:expr, $f2:expr) => {
        $v.chain($f0, $f1)
    };
}
macro_rules! chain2 {
    ($v:expr, $f0:expr, $f1:expr, $f2:expr) => {
        $v.chain($f0, $f1, $f2)
    };
}
scalar_impl!(Dual, chain1);
scalar_impl!(Dual2, chain2);
