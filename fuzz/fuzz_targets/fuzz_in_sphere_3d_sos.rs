#![no_main]

use libfuzzer_sys::fuzz_target;
use geogram_predicates::{Point3d, Sign};

fuzz_target!(|data: (Point3d, Point3d, Point3d, Point3d, Point3d)| {
    let (a, b, c, d, p) = data;
    match geogram_predicates::in_sphere_3d_sos::<true>(&a, &b, &c, &d, &p) {
        Sign::Zero => panic!("can't be zero"),
        _ => ()
    }

    match geogram_predicates::in_sphere_3d_sos::<false>(&a, &b, &c, &d, &p) {
        Sign::Zero => panic!("can't be zero"),
        _ => ()
    }
});
