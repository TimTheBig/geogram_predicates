//! # Geogram Predicates
//!
//! A crate for rust interoperability with `geogram`s _robust predicates_; via `cxx`.
#![no_std]

#[cfg(test)]
mod tests;

use core::cmp::Ordering;
pub use geogram_ffi::*;
use nalgebra::Point3;
use robust::{Coord, Coord3D};

mod expansion;
pub use expansion::Expansion;

pub type Point3d = [f64; 3];
pub type Point2d = [f64; 2];

/// A helper for functions that return signs
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone)]
#[repr(i8)]
pub enum Sign {
    Positive = 1,
    Zero = 0,
    Negative = -1,
}

impl From<Sign> for i8 {
    fn from(sign: Sign) -> i8 {
        sign as i8
    }
}

impl PartialEq<i8> for Sign {
    fn eq(&self, other: &i8) -> bool {
        (self.clone() as i8) == *other
    }
}

impl PartialEq<Sign> for i8 {
    fn eq(&self, other: &Sign) -> bool {
        (other.clone() as i8) == *self
    }
}

impl PartialOrd<i8> for Sign {
    fn partial_cmp(&self, other: &i8) -> Option<Ordering> {
        (self.clone() as i8).partial_cmp(other)
    }
}

/// Gets the sign of a value.
///
/// ### Parameters
/// - `x` value to test
///
/// # Example
/// ```
/// use geogram_predicates::{Sign, geo_sign};
///
/// let a = 42.0;
/// let b = -42.0;
/// let c = 0.0;
///
/// assert_eq!(geo_sign(a), Sign::Positive);
/// assert_eq!(geo_sign(b), Sign::Negative);
/// assert_eq!(geo_sign(c), Sign::Zero);
/// ```
#[inline]
pub const fn geo_sign(value: f64) -> Sign {
    if value > 0.0 {
        Sign::Positive
    } else if value < 0.0 {
        Sign::Negative
    } else {
        Sign::Zero
    }
}

/// Computes the orientation predicate in 3d.
///
/// Computes the sign of the signed volume of the tetrahedron `a`, `b`, `c`, `d`.
///
/// ### Parameters
/// - `a`, `b`, `c`, `d` vertices of the tetrahedron
///
/// ### Return values
/// * `+1` - if the tetrahedron is oriented positively
/// * `0` - if the tetrahedron is flat
/// * `-1` - if the tetrahedron is oriented negatively
///
/// # Example
/// ```
/// use geogram_predicates::orient_3d;
///
/// // Define four points that form a tetrahedron
/// let a = [0.0, 0.0, 0.0];
/// let b = [2.0, 0.0, 0.0];
/// let c = [0.0, 2.0, 0.0];
/// let d = [0.75, 0.75, 1.0];
///
/// assert_eq!(1, orient_3d(&a, &b, &c, &d));
///```
pub fn orient_3d(a: &Point3d, b: &Point3d, c: &Point3d, d: &Point3d) -> i8 {
    let orientation = robust::orient3d(
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*a) },
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*b) },
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*c) },
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*d) },
    );
    if orientation > 0.0 {
        -1
    } else if orientation < 0.0 {
        1
    } else {
        0
    }
}

/// Computes the orientation predicate in 2d.
///
/// Computes the sign of the signed area of the triangle `a`, `b`, `c`.
///
/// ### Parameters
/// - `a`, `b`, `c` vertices of the triangle
///
/// ### Return values
/// * `+1` - if the triangle is oriented counter-clockwise
/// * `0` - if the triangle is flat
/// * `-1` - if the triangle is oriented clockwise
///
/// # Example
/// ```
/// use geogram_predicates::orient_2d;
///
/// // Define three points that form a triangle
/// let a = [0.0, 0.0];
/// let b = [2.0, 0.0];
/// let c = [1.0, 1.0];
///
/// let orientation = orient_2d(&a, &b, &c);
/// assert_eq!(1, orientation);
/// ```
pub fn orient_2d(a: &Point2d, b: &Point2d, c: &Point2d) -> i8 {
    let orientation = robust::orient2d(
        unsafe { core::mem::transmute::<Point2d, Coord<f64>>(*a) },
        unsafe { core::mem::transmute::<Point2d, Coord<f64>>(*b) },
        unsafe { core::mem::transmute::<Point2d, Coord<f64>>(*c) },
    );
    if orientation > 0.0 {
        1
    } else if orientation < 0.0 {
        -1
    } else {
        0
    }
}

/// Computes the sign of the dot product between two vectors.
///
/// ### Parameters
/// - `a`, `b`, `c`, three 3d points
///
/// ### Returns
/// - the sign of the dot product between the vectors `ab` and `ac`, true is positive
///
/// # Example
/// ```
/// use geogram_predicates::dot_3d;
///
/// // Define four points that form a matrix
/// let a = [0.0, 0.0, 0.0];
/// let b = [1.0, 0.0, 0.0];
/// let c = [0.0, 1.0, 0.0];
///
/// let dot_sign = dot_3d(&a, &b, &c); // should be orthogonal
/// assert_eq!(dot_sign, true);
/// ```
pub fn dot_3d(a: &Point3d, b: &Point3d, c: &Point3d) -> bool {
    let ab_diff = Into::<Point3<_>>::into(*a) - Into::<Point3<_>>::into(*b);
    let ab_distance = ab_diff.x.hypot(ab_diff.y).hypot(ab_diff.z);
    let ac_diff = Into::<Point3<_>>::into(*a) - Into::<Point3<_>>::into(*c);
    let ac_distance = ac_diff.x.hypot(ac_diff.y).hypot(ac_diff.z);

    nalgebra::vector![ab_distance].dot(&nalgebra::vector![ac_distance]) >= 0.0
}

/// Tests whether a point is in the circum-circle of a triangle.
///
/// If the triangle `a` , `b` , `c` is oriented clockwise instead of counter-clockwise, then the result is inversed.
///
/// ### Parameters
/// - `a`, `b`, `c` vertices of the triangle
/// - `p` point to test
/// - `PERTURB` (const) - true is `+1` false is `-1` consistent perturbation
///
/// ### Return values
/// * `+1` - if `p` is inside the circum-circle of `a`, `b`, `c`
/// * `-1` - if `p` is outside the circum-circle of `a`, `b`, `c`
/// * `PERTURB` - if `p` is exactly on the circum-circle of the triangle `a`, `b`, `c`, where `perturb()` denotes a consistent perturbation, that returns either `+1` or `-1`
///
/// # Example
/// ```
/// use geogram_predicates as gp;
///
/// // Define three points that form a triangle
/// let a = [0.0, 0.0];
/// let b = [2.0, 0.0];
/// let c = [1.0, 1.0];
///
/// // Define two points, to test against the triangles circum-circle
/// let p_in = [1.0, -0.4];
/// let p_out = [1.0, -1.2];
///
/// let is_in_circle_p_in = gp::in_circle_2d_sos::<false>(&a, &b, &c, &p_in);
/// assert_eq!(1, is_in_circle_p_in);
/// # let is_in_circle_p_in = gp::in_circle_2d_sos::<true>(&a, &b, &c, &p_in);
/// # assert_eq!(1, is_in_circle_p_in);
///
/// let is_in_circle_p_out = gp::in_circle_2d_sos::<true>(&a, &b, &c, &p_out);
/// assert_eq!(-1, is_in_circle_p_out);
/// # let is_in_circle_p_out = gp::in_circle_2d_sos::<false>(&a, &b, &c, &p_out);
/// # assert_eq!(-1, is_in_circle_p_out);
/// ```
pub fn in_circle_2d_sos<const PERTURB: bool>(a: &Point2d, b: &Point2d, c: &Point2d, p: &Point2d) -> i8 {
    let incircle = robust::incircle(
        unsafe { core::mem::transmute::<Point2d, Coord<f64>>(*a) },
        unsafe { core::mem::transmute::<Point2d, Coord<f64>>(*b) },
        unsafe { core::mem::transmute::<Point2d, Coord<f64>>(*c) },
        unsafe { core::mem::transmute::<Point2d, Coord<f64>>(*p) },
    );

    if incircle > 0.0 {
        1
    } else if incircle < 0.0 {
        -1
    } else {
        const { if PERTURB { 1 } else { -1 } }
    }
}

/// Tests whether a point is in the circum-sphere of a tetrahedron.
///
/// ### Parameters
/// - `a`, `b`, `c`, `d` vertices of the tetrahedron
/// - `p` point to test
/// - `PERTURB` (const) - what it should be if `p` is exactly on the circum-sphere, true is `+1`, false is `-1`
///
/// ### Return values
/// * `+1` - if `p` is inside the circum-sphere of `a`, `b`, `c`, `d`
/// * `-1` - if `p` is outside the circum-sphere of `a`, `b`, `c`, `d`
/// * `PERTURB` - if `p` is exactly on the circum-sphere of the tetrahedron `a`, `b`, `c`, `d`, where `perturb()` denotes a consistent perturbation, that returns either `+1` or `-1`
///
/// # Example
/// ```
/// use geogram_predicates as gp;
///
/// // Define four points that form a tetrahedron
/// let a = [0.0, 0.0, 0.0];
/// let b = [2.0, 0.0, 0.0];
/// let c = [0.0, 2.0, 0.0];
/// let d = [0.75, 0.75, 1.0];
///
/// // Define two points, to test against the tetrahedrons circum-sphere
/// let p_in = [0.75, 0.75, 0.5];
/// assert_eq!(-1, gp::in_sphere_3d_sos::<false>(&a, &b, &c, &d, &p_in));
/// # assert_eq!(-1, gp::in_sphere_3d_sos::<true>(&a, &b, &c, &d, &p_in));
///
/// let p_out = [0.75, 0.75, 1.5];
/// assert_eq!(1, gp::in_sphere_3d_sos::<true>(&a, &b, &c, &d, &p_out));
/// # assert_eq!(1, gp::in_sphere_3d_sos::<false>(&a, &b, &c, &d, &p_out));
/// ```
pub fn in_sphere_3d_sos<const PERTURB: bool>(
    a: &Point3d,
    b: &Point3d,
    c: &Point3d,
    d: &Point3d,
    p: &Point3d,
) -> i8 {
    // let insphere = if orient_3d(a, b, c, d) == 1 {
    //     robust::insphere(
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*a) },
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*b) },
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*c) },
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*d) },
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*p) },
    //     )
    // } else {
    //     robust::insphere(
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*d) },
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*c) },
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*b) },
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*a) },
    //         unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*p) },
    //     )
    // };
    let insphere = robust::insphere(
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*a) },
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*b) },
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*c) },
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*d) },
        unsafe { core::mem::transmute::<Point3d, Coord3D<f64>>(*p) },
    );

    if insphere > 0.0 {
        1
    } else if insphere < 0.0 {
        -1
    } else {
        const { if PERTURB { 1 } else { -1 } }
    }
}

/// Computes the sign of the determinant of a 3x3 matrix formed by three 3D points.
///
/// ### Parameters
/// - `a`, `b`, `c` the three points that form the matrix
///
/// ### Returns
/// - the sign of the determinant of the matrix
///
/// # Example
/// ```
/// use geogram_predicates as gp;
///
/// // Define three points that form a matrix
/// let a = [1.0, 2.0, 3.0];
/// let b = [4.0, 5.0, 6.0];
/// let c = [7.0, 8.0, 9.0];
///
/// let det = gp::det_3d(&a, &b, &c);
/// assert_eq!(det, 0);
/// ```
pub const fn det_3d(a: &Point3d, b: &Point3d, c: &Point3d) -> i8 {
    det_3d_filter(a, b, c)
}

#[inline]
const fn det_3d_filter(p0: &Point3d, p1: &Point3d, p2: &Point3d) -> i8 {
    const FPG_UNCERTAIN_VALUE: i8 = 0;

    let delta = ((p0[0] * ((p1[1] * p2[2]) - (p1[2] * p2[1]))) - (p1[0] * ((p0[1] * p2[2]) - (p0[2] * p2[1]))))
        + (p2[0] * ((p0[1] * p1[2]) - (p0[2] * p1[1])));

    let max1 = p0[0].abs().max(p1[0].abs()).max(p2[0].abs());
    let max2 = p0[1].abs().max(p0[2].abs()).max(p1[1].abs()).max(p1[2].abs());
    let max3 = p1[1].abs().max(p1[2].abs()).max(p2[1].abs()).max(p2[2].abs());

    let mut lower_bound_1 = max1;
    let mut upper_bound_1 = max1;

    if max2 < lower_bound_1 {
        lower_bound_1 = max2;
    } else if max2 > upper_bound_1 {
        upper_bound_1 = max2;
    }
    if max3 < lower_bound_1 {
        lower_bound_1 = max3;
    } else if max3 > upper_bound_1 {
        upper_bound_1 = max3;
    }

    if lower_bound_1 < 1.92663387981871579179e-98 || upper_bound_1 > 1.11987237108890185662e+102 {
        FPG_UNCERTAIN_VALUE
    } else {
        let eps = 3.11133555671680765034e-15 * ((max2 * max3) * max1);
        if delta > eps {
            1
        } else if delta < -eps {
            -1
        } else {
            FPG_UNCERTAIN_VALUE
        }
    }
}

/// Computes the sign of the determinant of a 3x3
/// matrix formed by three 3d points using exact arithmetics.
///
/// ### Parameters
/// p0 , p1 , p2 the three points
///
/// ### Returns
/// The sign of the determinant of the matrix.
// fn det_3d_exact(p0: &Point3d, p1: &Point3d, p2: &Point3d) -> i8 {
//     let p0_0 = Expansion::from(p0[0]);
//     let p0_1 = Expansion::from(p0[1]);
//     let p0_2 = Expansion::from(p0[2]);

//     let p1_0 = Expansion::from(p1[0]);
//     let p1_1 = Expansion::from(p1[1]);
//     let p1_2 = Expansion::from(p1[2]);

//     let p2_0 = Expansion::from(p2[0]);
//     let p2_1 = Expansion::from(p2[1]);
//     let p2_2 = Expansion::from(p2[2]);

//     let Delta = expansion_det3x3!(
//         p0_0, p0_1, p0_2,
//         p1_0, p1_1, p1_2,
//         p2_0, p2_1, p2_2
//     );

//     return Delta.sign().into();
// }

#[cxx::bridge(namespace = "GEOGRAM")]
mod geogram_ffi {
    // Shared structs with fields visible to both languages.
    // ...

    // Rust types and signatures exposed to C++.
    // ...

    // C++ types and signatures exposed to Rust.
    unsafe extern "C++" {
        include!("geogram_predicates/include/geogram_ffi.h");

        /// Computes the sign of the determinant of a 4x4 matrix formed by four 4D points.
        ///
        /// ### Parameters
        /// - `a`, `b`, `c`, `d` the four points that form the matrix
        ///
        /// ### Returns
        /// - the sign of the determinant of the matrix
        ///
        /// # Example
        /// ```
        /// use geogram_predicates as gp;
        ///
        /// // Define four points that form a matrix
        /// let a = [1.0, 2.0, 3.0, 4.0];
        /// let b = [5.0, 6.0, 7.0, 8.0];
        /// let c = [9.0, 10.0, 11.0, 12.0];
        /// let d = [13.0, 14.0, 15.0, 16.0];
        ///
        /// let det = gp::det_4d(&a, &b, &c, &d);
        /// assert_eq!(det, 0);
        /// ```
        fn det_4d(a: &[f64; 4], b: &[f64; 4], c: &[f64; 4], d: &[f64; 4]) -> i16;

        /// Needs to be called before using any predicate.
        fn initialize();

        /// Computes the 3d orientation test with lifted points, i.e the regularity test for 2d.
        ///
        /// Given three lifted points `a'`, `b'`, `c'` in R^3, tests if the lifted point `p'` in R^3 lies below or above the plane passing through the three points `a'`, `b'`, `c'`.
        ///
        /// The coordinates and the heights are specified in separate arguments for each vertex.
        ///
        /// Note: if `w_i` = 0 this is equal to the in-circle test for the triangle `a`, `b`, `c` w.r.t `p`.
        ///
        /// ### Parameters
        /// - `a` ,`b`, `c`	vertices of the triangle
        /// - `p` point to test
        /// - `h_a` ,`h_b` ,`h_c` the heights of the lifted points, e.g. `a' = a.x**2 + a.y**2 - a.w`
        /// - `h_p` the height of the lifted point `p`
        ///
        /// ### Return values
        /// - `+1` - if p3' lies below the plane
        /// - `-1` - if p3' lies above the plane
        /// - perturb()	- if `p'` lies exactly on the hyperplane, where perturb() denotes a globally consistent perturbation, that returns either `+1` or `-1`
        ///
        /// # Example
        /// For a graphical representation see this [geogebra example](https://www.geogebra.org/m/etyzj96t) of the code below.
        /// ```
        /// use geogram_predicates as gp;
        ///
        /// // Define three points that form a triangle
        /// let a: [f64; 2] = [0.0, 0.0];
        /// let b: [f64; 2] = [2.0, 0.0];
        /// let c: [f64; 2] = [0.0, 2.0];
        ///
        /// // Additionally in this scenario, each point is associated with a weight w_i
        /// // And the height of a point is defined as h_i = x_i**2 + y_i**2 - w_i
        /// // One can interpret the height as the z-coordinate of a point lifted to R^3
        /// let h_a = a[0].powf(2.0) + a[1].powf(2.0) + 2.0;  // i.e. w_a = -2.0
        /// let h_b = b[0].powf(2.0) + b[1].powf(2.0) - 1.0;  // i.e. w_b = 1.0
        /// let h_c = c[0].powf(2.0) + c[1].powf(2.0) - 0.5;  // i.e. w_c = 0.5
        ///
        /// // Define weighted points, to test against the plane, that contains the lifted triangle
        /// let p_below: [f64; 2] = [0.6, 0.6];
        /// let h_p_below = p_below[0].powf(2.0) + p_below[1].powf(2.0) + 1.28;  // i.e. w_p_below = -1.28
        ///
        /// let p_above: [f64; 2] = [0.6, 0.6];
        /// let h_p_above = p_above[0].powf(2.0) + p_above[1].powf(2.0) + 2.78;  // i.e. w_p_above = -2.78
        ///
        /// let orientation_below = gp::orient_2dlifted_SOS(&a, &b, &c, &p_below, h_a, h_b, h_c, h_p_below);
        /// assert_eq!(1, orientation_below);
        ///
        /// let orientation_above = gp::orient_2dlifted_SOS(&a, &b, &c, &p_above, h_a, h_b, h_c, h_p_above);
        /// assert_eq!(-1, orientation_above);
        /// ```
        #[allow(clippy::too_many_arguments)]
        fn orient_2dlifted_SOS(
            a: &[f64; 2],
            b: &[f64; 2],
            c: &[f64; 2],
            p: &[f64; 2],
            h_a: f64,
            h_b: f64,
            h_c: f64,
            h_p: f64,
        ) -> i16;

        /// Computes the 4d orientation test with lifted points, i.e the regularity test for 3d.
        ///
        /// Given four lifted points `a'`, `b'`, `c'`, `d'` in R^4, tests if the lifted point `p'` in R^4 lies below or above the hyperplane passing through the four points `a'`, `b'`, `c'`, `d'`.
        ///
        /// Symbolic perturbation is applied, whenever the 5 vertices are not linearly independent.
        ///
        /// The coordinates and the heights are specified in separate arguments for each vertex.
        ///
        /// Note: if `w_i` = 0 this is equal to the in-sphere test for the tetrahedron `a`, `b`, `c`, `d` w.r.t `p`.
        ///
        /// ### Parameters
        /// - `a` ,`b`, `c`, `d` vertices of the tetrahedron
        /// - `p` point to test
        /// - `h_a` ,`h_b` ,`h_c`, `h_d` the heights of the lifted points, e.g. `h_a' = a.x**2 + a.y**2 - a.w`
        /// - `h_p` the height of the lifted point `p'`
        ///
        /// ### Return values
        /// - `+1` - if `p'` lies below the plane
        /// - `-1` - if `p'` lies above the plane
        /// - perturb()	- if `p'` lies exactly on the hyperplane, where perturb() denotes a globally consistent perturbation, that returns either `+1` or `-1`
        ///
        /// # Example
        /// ```
        /// use geogram_predicates as gp;
        /// // Define four points that form a tetrahedron
        /// let a: [f64; 3] = [0.0, 0.0, 0.0];
        /// let b: [f64; 3] = [2.0, 0.0, 0.0];
        /// let c: [f64; 3] = [0.0, 2.0, 0.0];
        /// let d: [f64; 3] = [0.75, 0.75, 1.0];
        ///
        /// // Additionally in this scenario, each point is associated with a weight w_i
        /// // And the height of a point is defined as h_i = x_i**2 + y_i**2 + z_i**2 - w_i
        /// // One can interpret the height as the 4th-coordinate of a point lifted to R^4
        /// let h_a = a[0].powf(2.0) + a[1].powf(2.0) + a[2].powf(2.0) + 2.0;  // i.e. w_a = -2.0
        /// let h_b = b[0].powf(2.0) + b[1].powf(2.0) + b[2].powf(2.0) - 1.0;  // i.e. w_b = 1.0
        /// let h_c = c[0].powf(2.0) + c[1].powf(2.0) + c[2].powf(2.0) - 0.5;  // i.e. w_c = 0.5
        /// let h_d = d[0].powf(2.0) + d[1].powf(2.0) + d[2].powf(2.0) - 0.5;  // i.e. w_c = 0.5
        ///
        /// // Define weighted points, to test against the hyperplane, that contains the lifted tetrahedron
        /// let p_below: [f64; 3] = [0.6, 0.6, 0.6];
        /// let h_p_below = p_below[0].powf(2.0) + p_below[1].powf(2.0) + 0.28;  // i.e. w_p_below = -0.28
        ///
        /// let p_above: [f64; 3] = [0.6, 0.6, 0.6];
        /// let h_p_above = p_above[0].powf(2.0) + p_above[1].powf(2.0) + 2.78;  // i.e. w_p_above = -2.78
        ///
        /// let orientation_below = gp::orient_3dlifted_SOS(&a, &b, &c, &d, &p_below, h_a, h_b, h_c, h_d, h_p_below);
        /// assert_eq!(1, orientation_below);
        ///
        /// let orientation_above = gp::orient_3dlifted_SOS(&a, &b, &c, &d, &p_above, h_a, h_b, h_c, h_d, h_p_above);
        /// assert_eq!(-1, orientation_above);
        ///
        /// ```
        #[allow(clippy::too_many_arguments)]
        fn orient_3dlifted_SOS(
            a: &[f64; 3],
            b: &[f64; 3],
            c: &[f64; 3],
            d: &[f64; 3],
            p: &[f64; 3],
            h_a: f64,
            h_b: f64,
            h_c: f64,
            h_d: f64,
            h_p: f64,
        ) -> i16;

        /// Tests whether three 3D points are colinear.
        ///
        /// ### Parameters
        /// - `p1` first point
        /// - `p2` second point
        /// - `p3` third point
        ///
        /// ### Return values
        /// - `true` - if `p1`, `p2` and `p3` are colinear
        /// - `false` - otherwise
        ///
        /// # Example
        /// ```
        /// use geogram_predicates as gp;
        ///
        /// // Define three points on a line
        /// let p1 = [0.0, 0.0, 0.0];
        /// let p2 = [0.0, 0.0, 1.0];
        /// let p3 = [0.0, 0.0, 2.0];
        ///
        /// assert!(gp::points_are_colinear_3d(&p1, &p2, &p3));
        /// ```
        fn points_are_colinear_3d(p1: &[f64; 3], p2: &[f64; 3], p3: &[f64; 3]) -> bool;

        /// Displays some statistics about predicates, including the number of calls, the number of exact arithmetics calls, and the number of Simulation of Simplicity calls.
        fn show_stats();

        /// Needs to be called at the end of the program.
        fn terminate();
    }
}

/// Tests whether two 2d points are identical.
///
/// ### Parameters
/// - `p1` first point
/// - `p2` second point
///
/// ### Return values
/// - `true` - if `p1` and `p2` have exactly the same coordinates
/// - `false` - otherwise
///
/// # Example
/// ```
/// use geogram_predicates::points_are_identical_2d;
///
/// let p1 = [4.0, 2.0];
/// let p2 = [4.0, 2.0];
///
/// assert!(points_are_identical_2d(&p1, &p2));
/// ```
pub const fn points_are_identical_2d(p1: &Point2d, p2: &Point2d) -> bool {
    p1[0] == p2[0] && p1[1] == p2[1]
}

/// Tests whether two 3d points are identical.
///
/// ### Parameters
/// - `p1` first point
/// - `p2` second point
///
/// ### Return values
/// - `true` - if `p1` and `p2` have exactly the same coordinates
/// - `false` - otherwise
///
/// # Example
/// ```
/// use geogram_predicates::points_are_identical_3d;
///
/// let p1 = [4.0, 2.0, 0.42];
/// let p2 = [4.0, 2.0, 0.42];
///
/// assert!(points_are_identical_3d(&p1, &p2));
/// ```
pub const fn points_are_identical_3d(p1: &Point3d, p2: &Point3d) -> bool {
    p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2]
}

/// Computes the (approximate) orientation predicate in 3d.
///
/// Computes the sign of the (approximate) signed volume of the tetrahedron `a`, `b`, `c`, `d`.
///
/// ### Parameters
/// - `a`, `b`, `c`, `d` vertices of the tetrahedron
///
/// ### Return values
/// * `+1` - if the tetrahedron is oriented positively
/// * `0` - if the tetrahedron is flat
/// * `-1` - if the tetrahedron is oriented negatively
///
/// # Example
/// ```
/// use geogram_predicates as gp;
///
/// // Define four points that form a tetrahedron
/// let a = [0.0, 0.0, 0.0];
/// let b = [2.0, 0.0, 0.0];
/// let c = [0.0, 2.0, 0.0];
/// let d = [0.75, 0.75, 1.0];
///
/// assert_eq!(1, gp::orient_3d_inexact(&a, &b, &c, &d));
///```
pub const fn orient_3d_inexact(a: &Point3d, b: &Point3d, c: &Point3d, d: &Point3d) -> Sign {
    let a11 = b[0] - a[0];
    let a12 = b[1] - a[1];
    let a13 = b[2] - a[2];

    let a21 = c[0] - a[0];
    let a22 = c[1] - a[1];
    let a23 = c[2] - a[2];

    let a31 = d[0] - a[0];
    let a32 = d[1] - a[1];
    let a33 = d[2] - a[2];

    let delta =  det3x3(
        [a11, a12, a13],
        [a21, a22, a23],
        [a31, a32, a33],
    );

    geo_sign(delta)
}

/// Computes the determinant of a 3x3 matrix given by coefficients.
#[inline]
const fn det3x3(
    [a00, a01, a02]: Point3d,
    [a10, a11, a12]: Point3d,
    [a20, a21, a22]: Point3d,
) -> f64 {
    let m01 = a00 * a11 - a10 * a01;
    let m02 = a00 * a21 - a20 * a01;
    let m12 = a10 * a21 - a20 * a11;

    m01 * a22 - m02 * a12 + m12 * a02
}
