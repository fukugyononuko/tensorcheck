//! Index labels for Einstein-convention expressions (API-08).
//!
//! Labels are compared by number. The labels here use numbers 0 to 39; define more
//! with [`label!`](crate::label) using numbers 40 to 63.

crate::label!(pub a = 0);
crate::label!(pub b = 1);
crate::label!(pub c = 2);
crate::label!(pub d = 3);
crate::label!(pub e = 4);
crate::label!(pub f = 5);
crate::label!(pub g = 6);
crate::label!(pub h = 7);
crate::label!(pub i = 8);
crate::label!(pub j = 9);
crate::label!(pub k = 10);
crate::label!(pub l = 11);
crate::label!(pub m = 12);
crate::label!(pub n = 13);
crate::label!(pub mu = 14);
crate::label!(pub nu = 15);
crate::label!(pub rho = 16);
crate::label!(pub sigma = 17);
crate::label!(pub alpha = 18);
crate::label!(pub beta = 19);
crate::label!(pub I = 20);
crate::label!(pub J = 21);
crate::label!(pub K = 22);
crate::label!(pub L = 23);
crate::label!(pub M = 24);
crate::label!(pub N = 25);
crate::label!(pub kappa_ = 26);
crate::label!(pub lambda = 27);
crate::label!(pub tau = 28);
crate::label!(pub gamma = 29);
crate::label!(pub delta = 30);
crate::label!(pub epsilon = 31);
crate::label!(pub p = 32);
crate::label!(pub q = 33);
crate::label!(pub r = 34);
crate::label!(pub s = 35);
crate::label!(pub P = 36);
crate::label!(pub Q = 37);
crate::label!(pub R = 38);
crate::label!(pub S = 39);

/// Numbers of the labels above, for the uniqueness test (TC-U13).
#[doc(hidden)]
pub const NUMBERS: [(&str, usize); 40] = [
    ("a", 0),
    ("b", 1),
    ("c", 2),
    ("d", 3),
    ("e", 4),
    ("f", 5),
    ("g", 6),
    ("h", 7),
    ("i", 8),
    ("j", 9),
    ("k", 10),
    ("l", 11),
    ("m", 12),
    ("n", 13),
    ("mu", 14),
    ("nu", 15),
    ("rho", 16),
    ("sigma", 17),
    ("alpha", 18),
    ("beta", 19),
    ("I", 20),
    ("J", 21),
    ("K", 22),
    ("L", 23),
    ("M", 24),
    ("N", 25),
    ("kappa_", 26),
    ("lambda", 27),
    ("tau", 28),
    ("gamma", 29),
    ("delta", 30),
    ("epsilon", 31),
    ("p", 32),
    ("q", 33),
    ("r", 34),
    ("s", 35),
    ("P", 36),
    ("Q", 37),
    ("R", 38),
    ("S", 39),
];
