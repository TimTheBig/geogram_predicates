#![no_main]

use libfuzzer_sys::fuzz_target;
use geogram_predicates::Point3d;

fuzz_target!(|data: (Point3d, Point3d, Point3d, Point3d, Point3d, [f64; 5])| {
    let (a, b, c, d, p, heights) = data;

    match geogram_predicates::orient_3dlifted_sos(&a, &b, &c, &d, &p, heights) {
        _ => ()
    }
});
