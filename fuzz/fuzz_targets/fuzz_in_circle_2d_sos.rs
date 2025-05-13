#![no_main]

use libfuzzer_sys::fuzz_target;
use geogram_predicates::{Point2d, Sign};

fuzz_target!(|data: (Point2d, Point2d, Point2d, Point2d)| {
    let (a, b, c, p) = data;
    match geogram_predicates::in_circle_2d_sos::<true>(&a, &b, &c, &p) {
        Sign::Zero => panic!("can't be zero"),
        _ => ()
    }

    match geogram_predicates::in_circle_2d_sos::<false>(&a, &b, &c, &p) {
        Sign::Zero => panic!("can't be zero"),
        _ => ()
    }
});
