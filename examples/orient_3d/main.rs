use geogram_predicates::orient_3d;
use std::{fs, path::Path};
use test_utils::{predicate_3d_test, write_to_png};

fn main() {
    let args = std::env::args().collect::<Vec<_>>();

    if args.len() != 2 {
        usage()
    }

    let mode = args[1].as_str();

    let p1 = [12., 12., 0.];
    let p2 = [24., 24., 12.];
    let p3 = [24., 24., 24.];

    let predicate: Box<dyn Fn([f64; 3]) -> f64> = match mode {
        "naive" => Box::new(|p| naive_orient_3d(&p1, &p, &p2, &p3)),
        "robust" => Box::new(|p| (orient_3d(&p1, &p, &p2, &p3) as i8).into()),
        "clean" => {
            let _ = fs::remove_file("out_naive_orient_3d.png");
            let _ = fs::remove_file("out_robust_orient_3d.png");
            println!("example images removed");
            std::process::exit(1);
        }
        "help" | _ => usage(),
    };

    let predicate_results = predicate_3d_test(predicate, [0.5, 0.5], 256, 256);

    let out_path = format!("out_{}_orient_3d.png", mode);

    write_to_png(&predicate_results, Path::new(&out_path), 256, 256);
}

// Directly evaluate the orient3d determinant (signed volume of a tetrahedron).
// Refer: https://www.cs.cmu.edu/~quake/robust.html
fn naive_orient_3d(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3], d: &[f64; 3]) -> f64 {
    let adx = a[0] - d[0];
    let ady = a[1] - d[1];
    let adz = a[2] - d[2];

    let bdx = b[0] - d[0];
    let bdy = b[1] - d[1];
    let bdz = b[2] - d[2];

    let cdx = c[0] - d[0];
    let cdy = c[1] - d[1];
    let cdz = c[2] - d[2];

    adx * (bdy * cdz - bdz * cdy)
  - ady * (bdx * cdz - bdz * cdx)
  + adz * (bdx * cdy - bdy * cdx)
}

fn usage() -> ! {
    eprintln!("
    Usage:
        orient_3d [option]

        MODES:
        naive - output an image showing the output of a naive orient_3d implementation
        robust - output an image showing the output of the robust orient_3d implementation
        OTHER:
        help - show this help message
        clean - remove the example output images
    ");
    std::process::exit(1);
}
