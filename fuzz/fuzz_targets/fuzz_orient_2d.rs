#![no_main]

use libfuzzer_sys::fuzz_target;
use geogram_predicates::Point2d;

fuzz_target!(|data: (Point2d, Point2d, Point2d)| {
    let (a, b, c) = data;
    assert_eq!(geogram_predicates::orient_2d(&a, &b, &c), geogram_predicates_legacy::orient_2d(&a, &b, &c) as i8);
    // match geogram_predicates::orient_2d(&a, &b, &c) {
    //     _ => ()
    // }
});
