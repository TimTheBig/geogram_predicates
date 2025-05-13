#![no_main]

use libfuzzer_sys::fuzz_target;
use geogram_predicates::Point3d;

fuzz_target!(|data: (Point3d, Point3d, Point3d, Point3d)| {
    let (a, b, c, d) = data;
    assert_eq!(geogram_predicates::orient_3d(&a, &b, &c, &d), geogram_predicates_legacy::orient_3d(&a, &b, &c, &d) as i8);
    // match geogram_predicates::orient_3d(&a, &b, &c, &d) {
    //     _ => ()
    // }
});
