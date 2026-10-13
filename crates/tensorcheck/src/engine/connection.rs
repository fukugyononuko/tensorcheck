//! The Levi-Civita connection and the Riemann tensor in a coordinate basis
//! (CMP-E09, DD-10).
//!
//! Pure arithmetic on flat arrays (`n` = dimension):
//! `ginv[a*n+d] = g^ad`, `dg[(c*n+a)*n+b] = ∂_c g_ab`,
//! `ddg[((c*n+d)*n+a)*n+b] = ∂_c∂_d g_ab`.

/// `Γ^a_bc = ½ g^ad (∂_b g_dc + ∂_c g_db − ∂_d g_bc)`, as `gamma[(a*n+b)*n+c]`.
pub fn christoffel(n: usize, ginv: &[f64], dg: &[f64]) -> Vec<f64> {
    let dgi = |c: usize, a: usize, b: usize| dg[(c * n + a) * n + b];
    let mut gamma = vec![0.0; n * n * n];
    for a in 0..n {
        for b in 0..n {
            for c in 0..n {
                let mut s = 0.0;
                for d in 0..n {
                    s += ginv[a * n + d] * (dgi(b, d, c) + dgi(c, d, b) - dgi(d, b, c));
                }
                gamma[(a * n + b) * n + c] = 0.5 * s;
            }
        }
    }
    gamma
}

/// `∂_e Γ^a_bc`, as `dgamma[((e*n+a)*n+b)*n+c]`, using
/// `∂_e g^ad = −g^af (∂_e g_fh) g^hd`.
pub fn christoffel_derivative(n: usize, ginv: &[f64], dg: &[f64], ddg: &[f64]) -> Vec<f64> {
    let dgi = |c: usize, a: usize, b: usize| dg[(c * n + a) * n + b];
    let ddgi = |c: usize, d: usize, a: usize, b: usize| ddg[((c * n + d) * n + a) * n + b];
    // ∂_e g^ad
    let mut dginv = vec![0.0; n * n * n];
    for e in 0..n {
        for a in 0..n {
            for d in 0..n {
                let mut s = 0.0;
                for f in 0..n {
                    for h in 0..n {
                        s += ginv[a * n + f] * dgi(e, f, h) * ginv[h * n + d];
                    }
                }
                dginv[(e * n + a) * n + d] = -s;
            }
        }
    }
    let mut out = vec![0.0; n * n * n * n];
    for e in 0..n {
        for a in 0..n {
            for b in 0..n {
                for c in 0..n {
                    let mut s = 0.0;
                    for d in 0..n {
                        let first = dgi(b, d, c) + dgi(c, d, b) - dgi(d, b, c);
                        let second = ddgi(e, b, d, c) + ddgi(e, c, d, b) - ddgi(e, d, b, c);
                        s += dginv[(e * n + a) * n + d] * first + ginv[a * n + d] * second;
                    }
                    out[((e * n + a) * n + b) * n + c] = 0.5 * s;
                }
            }
        }
    }
    out
}

/// `R^a_bcd = ∂_c Γ^a_db − ∂_d Γ^a_cb + Γ^a_ce Γ^e_db − Γ^a_de Γ^e_cb`, as
/// `riemann[((a*n+b)*n+c)*n+d]`.
pub fn riemann(n: usize, gamma: &[f64], dgamma: &[f64]) -> Vec<f64> {
    let g = |a: usize, b: usize, c: usize| gamma[(a * n + b) * n + c];
    let dg = |e: usize, a: usize, b: usize, c: usize| dgamma[((e * n + a) * n + b) * n + c];
    let mut r = vec![0.0; n * n * n * n];
    for a in 0..n {
        for b in 0..n {
            for c in 0..n {
                for d in 0..n {
                    let mut s = dg(c, a, d, b) - dg(d, a, c, b);
                    for e in 0..n {
                        s += g(a, c, e) * g(e, d, b) - g(a, d, e) * g(e, c, b);
                    }
                    r[((a * n + b) * n + c) * n + d] = s;
                }
            }
        }
    }
    r
}
