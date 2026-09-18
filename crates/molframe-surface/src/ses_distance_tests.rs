use super::distance_from_exterior;
use crate::cavity::{EMPTY, EXTERIOR, Grid};
use crate::numeric::{f64_to_f32, usize_to_f64};

#[test]
fn transform_returns_exact_axis_and_diagonal_distances() {
    let grid = Grid {
        origin: [0.0; 3],
        step: 1.0,
        dims: [3, 3, 3],
    };
    let mut state = vec![EMPTY; 27];
    state[grid.index(1, 1, 1)] = EXTERIOR;
    let distance = distance_from_exterior(&grid, state).expect("the small transform must fit");

    assert!((distance[grid.index(2, 1, 1)] - 1.0).abs() < f32::EPSILON);
    assert!((distance[grid.index(2, 2, 1)] - 2.0_f32.sqrt()).abs() < 1.0e-6);
    assert!((distance[grid.index(2, 2, 2)] - 3.0_f32.sqrt()).abs() < 1.0e-6);
}

#[test]
fn transform_preserves_zero_at_every_exterior_cell() {
    let grid = Grid {
        origin: [0.0; 3],
        step: 0.25,
        dims: [3, 3, 3],
    };
    let state = vec![EXTERIOR; 27];
    let distance = distance_from_exterior(&grid, state).expect("the small transform must fit");
    assert!(distance.iter().all(|value| value.abs() < f32::EPSILON));
}

#[test]
fn separable_transform_matches_brute_force_for_multiple_sources() {
    let grid = Grid {
        origin: [0.0; 3],
        step: 0.5,
        dims: [4, 3, 2],
    };
    let sources = [[0, 0, 0], [3, 2, 1]];
    let mut state = vec![EMPTY; 24];
    for [x, y, z] in sources {
        state[grid.index(x, y, z)] = EXTERIOR;
    }
    let distance = distance_from_exterior(&grid, state).expect("the small transform must fit");

    for z in 0..grid.dims[2] {
        for y in 0..grid.dims[1] {
            for x in 0..grid.dims[0] {
                let expected_squared = sources
                    .iter()
                    .map(|source| {
                        let dx = x.abs_diff(source[0]);
                        let dy = y.abs_diff(source[1]);
                        let dz = z.abs_diff(source[2]);
                        f64_to_f32(usize_to_f64(dx * dx + dy * dy + dz * dz))
                    })
                    .fold(f32::INFINITY, f32::min);
                let expected = expected_squared.sqrt() * 0.5;
                assert!((distance[grid.index(x, y, z)] - expected).abs() < 1.0e-6);
            }
        }
    }
}
