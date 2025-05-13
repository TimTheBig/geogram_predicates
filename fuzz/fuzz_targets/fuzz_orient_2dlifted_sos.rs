#![no_main]

use libfuzzer_sys::fuzz_target;
use geogram_predicates::Point2d;

fuzz_target!(|data: (Point2d, Point2d, Point2d, Point2d, [f64; 4])| {
    let (a, b, c, p, heights) = data;

    // fails if input contains nan, geogram is failing not us
    assert_eq!(
        geogram_predicates::orient_2dlifted_sos(&a, &b, &c, &p, heights),
        geogram_predicates_legacy::orient_2dlifted_SOS(&a, &b, &c, &p, heights[0], heights[1], heights[2], heights[3]) as i8
    )
});
