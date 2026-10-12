//! Phase 2: slots, storage, labels and contraction.
//! TC-U03, TC-U04, TC-U05, TC-U13, TC-I01, TC-R03.

use tensorcheck::labels::{I, J, a, b, c, d};
use tensorcheck::spaces::euclid3::{Cartesian, Euclid3};
use tensorcheck::spaces::plane::{Polar, R, Theta};
use tensorcheck::units::si::{
    Acceleration, Dimensionless, Force, Length, N, Velocity, kg, m, one, rad, s,
};
use tensorcheck::utility::typelevel::{Cons, Nil};
use tensorcheck::{Down, Quantity, Tensor, Up, Vector};

type T2<C1, C2> = Tensor<Dimensionless, Cons<C1, Cons<C2, Nil>>>;
type UpC = Up<Cartesian>;
type DnC = Down<Cartesian>;

/// Deterministic test values.
fn val(i: usize) -> f64 {
    ((i * 7 + 3) % 11) as f64 - 4.5
}

#[test]
fn storage_size_is_components_times_scalar() {
    // TC-U03 / NFR-03.
    assert_eq!(size_of::<Vector<Length, Cartesian>>(), 3 * 8);
    assert_eq!(size_of::<T2<UpC, DnC>>(), 9 * 8);
    tensorcheck::space!(Five, 5);
    tensorcheck::chart!(F5 on Five { X1: Length, X2: Length, X3: Length, X4: Length, X5: Length });
    type R4 =
        Tensor<Dimensionless, Cons<Up<F5>, Cons<Down<F5>, Cons<Down<F5>, Cons<Down<F5>, Nil>>>>>;
    assert_eq!(size_of::<R4>(), 625 * 8);
    type Mixed = Tensor<Dimensionless, Cons<Down<Cartesian>, Cons<Up<Polar>, Nil>>>;
    assert_eq!(size_of::<Mixed>(), 6 * 8);
}

#[test]
fn rank_two_components_round_trip() {
    let mut arr = [[0.0; 3]; 3];
    for (i, row) in arr.iter_mut().enumerate() {
        for (j, x) in row.iter_mut().enumerate() {
            *x = val(i * 3 + j);
        }
    }
    let t: T2<UpC, DnC> = Tensor::from_components(arr, one);
    assert_eq!(t.components_in(one), arr);
}

fn tensor2(off: usize) -> ([[f64; 3]; 3], T2<UpC, DnC>) {
    let arr = core::array::from_fn(|i| core::array::from_fn(|j| val(off + i * 3 + j)));
    (arr, Tensor::from_components(arr, one))
}

#[test]
fn contraction_matches_hand_loops() {
    // TC-U04: T^a_b V^b.
    let (t, tt) = tensor2(0);
    let v = [1.0, -2.0, 0.5];
    let vv: Vector<Dimensionless, Cartesian> = Vector::from_components(v, one);
    let w = (tt.idx((a, b)) * vv.idx(b)).into_tensor();
    let mut want = [0.0; 3];
    for (i, wi) in want.iter_mut().enumerate() {
        for j in 0..3 {
            *wi += t[i][j] * v[j];
        }
    }
    assert_eq!(w.components_in(one), want);

    // Two pairs: S^a_b T^b_a (a full contraction to a scalar).
    let (s2, ss) = tensor2(5);
    let full = (ss.idx((a, b)) * tt.idx((b, a))).into_tensor();
    let mut want = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            want += s2[i][j] * t[j][i];
        }
    }
    assert_eq!(full.value_in(one), want);

    // Free indices keep their order: S^a_b T^b_c -> (a, c).
    let prod = (ss.idx((a, b)) * tt.idx((b, c))).into_tensor();
    let got = prod.components_in(one);
    for i in 0..3 {
        for k in 0..3 {
            let want: f64 = (0..3).map(|j| s2[i][j] * t[j][k]).sum();
            assert_eq!(got[i][k], want);
        }
    }
}

#[test]
fn trace_and_outer_product() {
    let (t, tt) = tensor2(2);
    let tr = tt.idx((a, a)).into_tensor();
    assert_eq!(tr.value_in(one), t[0][0] + t[1][1] + t[2][2]);

    let u: Vector<Length, Cartesian> = Vector::from_components([1.0, 2.0, 3.0], m);
    let w: Vector<Length, Cartesian> = Vector::from_components([4.0, 5.0, 6.0], m);
    let o = u.outer(w);
    let oc = o.components_in(m * m);
    assert_eq!(oc[1][2], 2.0 * 6.0);
    // Contracting the outer product with a covector equals scaling.
    let k: tensorcheck::Covector<Dimensionless, Cartesian> =
        tensorcheck::Covector::from_components([1.0, 0.0, 1.0], one);
    let r = (o.idx((a, b)) * k.idx(b)).into_tensor();
    assert_eq!(r.components_in(m * m), [10.0, 20.0, 30.0]);
}

#[test]
fn chained_products() {
    // T^a_b T^b_c T^c_d contracted step by step equals the matrix cube.
    let (t, tt) = tensor2(1);
    let cube = (tt.idx((a, b)) * tt.idx((b, c)) * tt.idx((c, d))).into_tensor();
    let got = cube.components_in(one);
    for i in 0..3 {
        for l in 0..3 {
            let mut want = 0.0;
            for j in 0..3 {
                for k in 0..3 {
                    want += t[i][j] * t[j][k] * t[k][l];
                }
            }
            assert!((got[i][l] - want).abs() <= 1e-12 * want.abs().max(1.0));
        }
    }
}

#[test]
fn component_dimensions_in_polar_coordinates() {
    // TC-U05: V^r has dimension [V][r] = m/s, V^θ has [V][θ] = 1/s.
    let v: Vector<tensorcheck::Quot<Dimensionless, tensorcheck::units::si::Time>, Polar> =
        Vector::from_components([1.5, 0.25], one / s);
    let vr: Quantity<Velocity> = v.at((R,));
    let vt: Quantity<tensorcheck::units::si::Frequency> = v.at((Theta,));
    assert_eq!(vr.value_in(m / s), 1.5);
    assert_eq!(vt.value_in(rad / s), 0.25);
}

#[test]
fn labels_have_distinct_numbers() {
    // TC-U13.
    let mut seen = std::collections::HashSet::new();
    for (name, n) in tensorcheck::labels::NUMBERS {
        assert!(seen.insert(n), "label {name} reuses number {n}");
    }
}

#[test]
fn newtons_second_law_for_vectors() {
    // TC-I01 / CC-01.
    let mass = 2.0 * kg;
    let acc: Vector<Acceleration, Cartesian> =
        Vector::from_components([0.0, 0.0, -9.8], m / s.pow::<2>());
    let force: Vector<Force, Cartesian> = mass * acc;
    let f = force.components_in(N);
    assert_eq!(f, [0.0, 0.0, -19.6]);
}

// TC-R03 (VOL-03): a three-component field space defined by the user.
tensorcheck::space!(Fields3, 3);
tensorcheck::chart!(Phi3 on Fields3 { P1: Dimensionless, P2: Dimensionless, P3: Dimensionless });

#[test]
fn user_defined_field_space() {
    let g: Tensor<Dimensionless, Cons<Down<Phi3>, Cons<Down<Phi3>, Nil>>> =
        Tensor::from_components([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]], one);
    let dphi: Vector<Dimensionless, Phi3> = Vector::from_components([1.0, 1.0, 2.0], one);
    let kin = (g.idx((I, J)) * dphi.idx(I) * dphi.idx(J)).into_tensor();
    assert_eq!(kin.value_in(one), 1.0 + 2.0 + 12.0);
    let _ = Euclid3;
}
