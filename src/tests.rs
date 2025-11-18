use crate::{self as gp, Sign};

// Helper to initialize and terminate geogram
#[cfg(feature = "legacy")]
fn setup() {
    // initialize the C++ geogram library
    gp::initialize();
}

#[test]
fn test_geo_sgn() {
    assert_eq!(gp::geo_sign(42.0), Sign::Positive);
    assert_eq!(gp::geo_sign(-3.14), Sign::Negative);
    assert_eq!(gp::geo_sign(0.0), Sign::Zero);
}

#[test]
fn test_orient_2d() {
    // CCW
    let a = [0.0, 0.0];
    let b = [1.0, 0.0];
    let c = [0.0, 1.0];
    assert_eq!(gp::orient_2d(&a, &b, &c), 1);
    // CW
    assert_eq!(gp::orient_2d(&a, &c, &b), -1);
    // Co-linear
    let d = [2.0, 0.0];
    assert_eq!(gp::orient_2d(&a, &b, &d), 0);
}

#[test]
fn test_orient_3d() {
    let a = [0.0, 0.0, 0.0];
    let b = [1.0, 0.0, 0.0];
    let c = [0.0, 1.0, 0.0];
    let d = [0.0, 0.0, 1.0];
    // Right-handed tetrahedron
    assert_eq!(gp::orient_3d(&a, &b, &c, &d), 1);
    // Flip two vertices for negative orientation
    assert_eq!(gp::orient_3d(&a, &c, &b, &d), -1);
    // Flat
    let e = [0.5, 0.5, 0.0];
    assert_eq!(gp::orient_3d(&a, &b, &c, &e), 0);
}

#[test]
fn test_dot_3d() {
    let a = [0.0, 0.0, 0.0];
    let b = [1.0, 0.0, 0.0];
    let c = [0.0, 1.0, 0.0];
    // orthogonal
    assert_eq!(gp::dot_3d(&a, &b, &c), true);
    // same direction
    let d = [2.0, 0.0, 0.0];
    assert_eq!(gp::dot_3d(&a, &b, &d), true);
}

#[test]
fn test_points_identical() {
    let p1 = [4.0, 2.0];
    let p2 = [4.0, 2.0];
    let p3 = [4.0, 2.1];
    assert!(gp::points_are_identical_2d(&p1, &p2));
    assert!(!gp::points_are_identical_2d(&p1, &p3));

    let q1 = [1.0, 2.0, 3.0];
    let q2 = [1.0, 2.0, 3.0];
    let q3 = [1.0, 2.0, 3.1];
    assert!(gp::points_are_identical_3d(&q1, &q2));
    assert!(!gp::points_are_identical_3d(&q1, &q3));
}

#[test]
fn test_points_colinear_3d() {
    let p1 = [0.0, 0.0, 0.0];
    let p2 = [1.0, 1.0, 1.0];
    let p3 = [2.0, 2.0, 2.0];
    assert!(gp::points_are_colinear_3d(&p1, &p2, &p3));
    let p4 = [1.0, 0.0, 0.0];
    assert!(!gp::points_are_colinear_3d(&p1, &p2, &p4));
}

#[test]
fn test_in_circle_2d() {
    let a = [0.0, 0.0];
    let b = [1.0, 0.0];
    let c = [0.0, 1.0];
    let p_in = [0.1, 0.1];
    let p_out = [2.0, 2.0];

    assert_eq!(gp::in_circle_2d_sos::<false>(&a, &b, &c, &p_in), 1);
    assert_eq!(gp::in_circle_2d_sos::<true>(&a, &b, &c, &p_in), 1);
    assert_eq!(gp::in_circle_2d_sos::<false>(&a, &b, &c, &p_out), -1);
    assert_eq!(gp::in_circle_2d_sos::<true>(&a, &b, &c, &p_out), -1);
}

// On-border tests for incircle: PERTURB false -> -1, true -> +1
#[test]
fn test_in_circle_2d_on_border() {
    let a = [0.0, 0.0];
    let b = [1.0, 0.0];
    let c = [0.0, 1.0];
    // circumcenter = (0.5,0.5), radius = sqrt(0.5)
    let r = 0.5_f64.sqrt();
    // choose a point at angle π/4 around the circumcenter
    let p_on = [0.5 + r / 2.0_f64.sqrt(), 0.5 + r / 2.0_f64.sqrt()];
    // exactly on circle
    assert_eq!(gp::in_circle_2d_sos::<false>(&a, &b, &c, &p_on), -1);
    assert_eq!(gp::in_circle_2d_sos::<true>(&a, &b, &c, &p_on), 1);
}

#[test]
fn test_in_sphere_3d() {
    let a = [0.0, 0.0, 0.0];
    let b = [1.0, 0.0, 0.0];
    let c = [0.0, 1.0, 0.0];
    let d = [0.0, 0.0, 1.0];
    let p_in = [0.1, 0.1, 0.1];
    let p_out = [2.0, 2.0, 2.0];

    assert_eq!(gp::in_sphere_3d_sos::<false>(&a, &b, &c, &d, &p_in), -1);
    assert_eq!(gp::in_sphere_3d_sos::<true>(&a, &b, &c, &d, &p_in), -1);
    assert_eq!(gp::in_sphere_3d_sos::<false>(&a, &b, &c, &d, &p_out), 1);
    assert_eq!(gp::in_sphere_3d_sos::<true>(&a, &b, &c, &d, &p_out), 1);
}

// On-border tests for insphere: PERTURB false -> -1, true -> +1
#[test]
fn test_in_sphere_3d_on_border() {
    let a = [0.0, 0.0, 0.0];
    let b = [1.0, 0.0, 0.0];
    let c = [0.0, 1.0, 0.0];
    let d = [0.0, 0.0, 1.0];
    // circumsphere center = (0.25,0.25,0.25), radius = sqrt(3)/4
    let r = (3.0_f64).sqrt() / 4.0;
    let offset = r / (3.0_f64).sqrt();
    let p_on = [0.25 + offset, 0.25 + offset, 0.25 + offset];

    assert_eq!(gp::in_sphere_3d_sos::<false>(&a, &b, &c, &d, &p_on), -1);
    assert_eq!(gp::in_sphere_3d_sos::<true>(&a, &b, &c, &d, &p_on), -1);
}

#[cfg(feature = "legacy")]
#[test]
fn test_det_4d() {
    setup();
    let a4 = [1.0, 2.0, 3.0, 4.0];
    let b4 = [5.0, 6.0, 7.0, 8.0];
    let c4 = [9.0, 10.0, 11.0, 12.0];
    let d4 = [13.0, 14.0, 15.0, 16.0];
    assert_eq!(gp::det_4d(&a4, &b4, &c4, &d4), 0);
    gp::terminate();
}

#[cfg(feature = "legacy")]
#[test]
fn test_det_3d() {
    setup();
    // Linearly dependent vectors -> determinant should be 0
    let a3 = [1.0, 2.0, 3.0];
    let b3 = [4.0, 5.0, 6.0];
    let c3 = [7.0, 8.0, 9.0];
    assert_eq!(gp::det_3d(&a3, &b3, &c3), 0); // this passes with geogram

    // Standard basis vectors -> determinant should be 1
    let a3 = [1.0, 0.0, 0.0];
    let b3 = [0.0, 1.0, 0.0];
    let c3 = [0.0, 0.0, 1.0];
    assert_eq!(gp::det_3d(&a3, &b3, &c3), 1);

    // Negated basis vector -> determinant should be -1
    let a3 = [1.0, 0.0, 0.0];
    let b3 = [0.0, -1.0, 0.0];
    let c3 = [0.0, 0.0, 1.0];
    assert_eq!(gp::det_3d(&a3, &b3, &c3), -1);

    // Swapped rows should negate the determinant
    let a3 = [0.0, 1.0, 0.0];
    let b3 = [1.0, 0.0, 0.0];
    let c3 = [0.0, 0.0, 1.0];
    assert_eq!(gp::det_3d(&a3, &b3, &c3), -1);

    // Determinant of identity matrix should be 1
    let a3 = [1.0, 0.0, 0.0];
    let b3 = [0.0, 1.0, 0.0];
    let c3 = [0.0, 0.0, 1.0];
    assert_eq!(gp::det_3d(&a3, &b3, &c3), 1);

    let a = [1.0, 0.0, 0.0];
    let b = [0.0, 1.0, 0.0];
    let c = [1.0, 1.0, 0.0];
    assert_eq!(gp::det_3d(&a, &b, &c), 0);

    let a = [1e-10, 0.0,   0.0];
    let b = [0.0,   1e-10, 0.0];
    let c = [0.0,   0.0,  1e-10];
    assert_eq!(gp::det_3d(&a, &b, &c), 1);

    let a = [1e25, 0.0,  0.0];
    let b = [0.0,  1e25, 0.0];
    let c = [0.0,  0.0, -1e25];
    assert_eq!(gp::det_3d(&a, &b, &c), -1);
    gp::terminate();
}

#[test]
fn test_orient_2dlifted() {
    // 2D lifted: trivial case with zero weights equals incircle
    let a2 = [0.0, 0.0];
    let b2 = [1.0, 0.0];
    let c2 = [0.0, 1.0];
    let p2 = [0.1, 0.1];
    let h = 0.0;
    let res = gp::orient_2dlifted_sos(&a2, &b2, &c2, &p2, [h, h, h, h]);
    assert_eq!(res, gp::in_circle_2d_sos::<false>(&a2, &b2, &c2, &p2));
}

#[test]
fn test_orient_3dlifted() {
    // 3D lifted: trivial with zero weights equals insphere
    let a3 = [0.0, 0.0, 0.0];
    let b3 = [1.0, 0.0, 0.0];
    let c3 = [0.0, 1.0, 0.0];
    let d3 = [0.0, 0.0, 1.0];
    let p3 = [0.1, 0.1, 0.1];
    let h3 = 0.0;
    let res3 = gp::orient_3dlifted_sos(&a3, &b3, &c3, &d3, &p3, [h3, h3, h3, h3, h3]);
    assert_eq!(res3, 1 /* res from c++ */);

    assert_eq!(res3, gp::in_sphere_3d_sos::<false>(&a3, &b3, &c3, &d3, &p3));

    let a4 = [0.0, 3.0, 0.0];
    let res4 = gp::orient_3dlifted_sos(&a4, &b3, &c3, &d3, &p3, [h3, h3, h3, h3, h3]);
    assert_eq!(res4, 1 /* res from c++ */);

    let a = [0.0, 3.0, 0.0];
    let b = [-1.0, 0.0, 0.0];
    let c = [0.0, -1.0, 0.0];
    let d = [-9.0, 0.0, -1.0];
    let p = [-0.1, -0.1, -0.1];
    let res5 = gp::orient_3dlifted_sos(
        &a, &b, &c, &d, &p,
        [2.0, 0.1, -0.1, -3.0, 0.0]
    );
    assert_eq!(res5, 1 /* res from c++ */);

    // Case: zero weights, point inside the sphere
    let a = [1.0, 0.0, 0.0];
    let b = [0.0, 1.0, 0.0];
    let c = [0.0, 0.0, 1.0];
    let d = [1.0, 1.0, 1.0];
    let p = [0.5, 0.5, 0.5];
    let weights = [0.0; 5];
    let res = gp::orient_3dlifted_sos(&a, &b, &c, &d, &p, weights);
    assert_eq!(res, 1);
    assert_eq!(res, gp::in_sphere_3d_sos::<false>(&a, &b, &c, &d, &p));

    // Case: non-zero weights, point outside sphere
    let weights = [1.0, 2.0, 3.0, 4.0, 10.0];
    let res = gp::orient_3dlifted_sos(&a, &b, &c, &d, &p, weights);
    assert_eq!(res, -1);

    // Case: point lies on the lifted sphere (weights chosen for coplanarity)
    let weights = [1.0, 1.0, 1.0, 1.0, 1.0];
    let res = gp::orient_3dlifted_sos(&a, &b, &c, &d, &p, weights);
    assert_eq!(res, 0);

    // Degenerate case: repeated point
    let res = gp::orient_3dlifted_sos(&a, &a, &c, &d, &p, weights);
    assert_eq!(res, 0); // Expect zero or handled degeneracy

    // Degenerate case: flat tetrahedron (coplanar points)
    let a = [0.0, 0.0, 0.0];
    let b = [1.0, 0.0, 0.0];
    let c = [0.0, 1.0, 0.0];
    let d = [1.0, 1.0, 0.0];
    let p = [0.5, 0.5, 0.0];
    let weights = [0.0, 0.0, 0.0, 0.0, 0.0];
    let res = gp::orient_3dlifted_sos(&a, &b, &c, &d, &p, weights);
    assert_eq!(res, 0); // Planar configuration should be degenerate

    // Case: large weights, test numerical stability
    let weights = [1e9, 1e9, 1e9, 1e9, 1e9];
    let res = gp::orient_3dlifted_sos(&a, &b, &c, &d, &p, weights);
    // Result depends on actual implementation; this checks for crash/stability
    assert!((res as i8).abs() <= 1);
}

#[test]
fn test_orient_3d_inexact() {
    let a = [0.0, 0.0, 0.0];
    let b = [1.0, 0.0, 0.0];
    let c = [0.0, 1.0, 0.0];
    let d = [0.0, 0.0, 1.0];

    assert_eq!(gp::orient_3d_inexact(&a, &b, &c, &d), Sign::Positive);
}
