use super::{EMPTY, EXTERIOR, FloodWorkspace, Grid, SOLID, cavities, cavities_with_options};
use crate::{SasaError, SurfaceGridOptions};
use crate::{numeric::f64_to_f32, sampling::fibonacci_sphere};

#[test]
fn a_lone_atom_encloses_nothing() {
    let found = cavities(&[[0.0, 0.0, 0.0]], &[2.0], 0.0, 0.5);
    let Ok(found) = found else { panic!("valid") };
    assert!(found.is_empty());
}

#[test]
fn a_sealed_shell_of_atoms_encloses_one_cavity() {
    // A dense shell on a sphere of radius 5, each atom large enough that its
    // neighbours overlap, seals a hollow interior.
    let positions: Vec<[f32; 3]> = fibonacci_sphere(120)
        .into_iter()
        .map(|direction| {
            [
                f64_to_f32(direction[0] * 5.0),
                f64_to_f32(direction[1] * 5.0),
                f64_to_f32(direction[2] * 5.0),
            ]
        })
        .collect();
    let radii = vec![2.2f32; positions.len()];
    let Ok(found) = cavities(&positions, &radii, 0.0, 0.6) else {
        panic!("valid");
    };
    assert!(!found.is_empty(), "the hollow interior should be a cavity");
    assert!(found[0].volume > 0.0);
    // The largest cavity's representative point sits near the centre.
    let point = found[0].representative;
    let distance =
        (f64::from(point[0]).powi(2) + f64::from(point[1]).powi(2) + f64::from(point[2]).powi(2))
            .sqrt();
    assert!(
        distance < 3.0,
        "cavity point should be interior, was {distance}"
    );
}

#[test]
fn a_non_positive_resolution_is_rejected() {
    assert!(cavities(&[[0.0, 0.0, 0.0]], &[2.0], 0.0, 0.0).is_err());
}

#[test]
fn allocation_ceiling_is_caller_controlled() {
    let error = cavities_with_options(
        &[[0.0, 0.0, 0.0]],
        &[2.0],
        0.0,
        SurfaceGridOptions {
            resolution: 0.5,
            max_cells: 1,
            max_workspace_bytes: SurfaceGridOptions::STANDARD_WORKSPACE_BYTES,
        },
    )
    .expect_err("one cell cannot hold the bounded grid");
    assert!(matches!(error, SasaError::GridTooLarge { .. }));
    assert!(matches!(
        cavities_with_options(
            &[[0.0, 0.0, 0.0]],
            &[2.0],
            0.0,
            SurfaceGridOptions {
                resolution: 0.5,
                max_cells: 0,
                max_workspace_bytes: SurfaceGridOptions::STANDARD_WORKSPACE_BYTES,
            },
        ),
        Err(SasaError::InvalidGridOptions)
    ));
}

#[test]
fn the_stress_fixture_has_exactly_one_million_cells() {
    let options = SurfaceGridOptions {
        resolution: 0.1,
        max_cells: 1_000_000,
        max_workspace_bytes: SurfaceGridOptions::STANDARD_WORKSPACE_BYTES,
    };
    let grid = Grid::new(&[[0.0; 3]], &[4.85], options).expect("the exact grid must fit");
    assert_eq!(grid.dims, [100, 100, 100]);
    assert_eq!(grid.cell_count(), 1_000_000);
}

#[test]
fn scanline_flood_matches_a_cell_reference_on_fragmented_grids() {
    let grid = Grid {
        origin: [0.0; 3],
        step: 1.0,
        dims: [4, 4, 4],
    };
    let mut seed = 0x9e37_79b9_u32;
    for _ in 0..128 {
        let mut state = vec![EMPTY; grid.cell_count()];
        for cell in &mut state {
            seed = seed.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            if seed.trailing_zeros() >= 2 {
                *cell = SOLID;
            }
        }
        let mut expected = state.clone();
        reference_exterior(&grid, &mut expected);
        let mut flood = FloodWorkspace::new(grid.cell_count()).expect("small workspace must fit");
        grid.flood_exterior(&mut state, &mut flood);
        assert_eq!(state, expected);
    }
}

fn reference_exterior(grid: &Grid, state: &mut [u8]) {
    let mut stack = Vec::new();
    for index in 0..state.len() {
        let plane = grid.dims[0] * grid.dims[1];
        let z = index / plane;
        let within_plane = index % plane;
        let y = within_plane / grid.dims[0];
        let x = within_plane % grid.dims[0];
        if x == 0
            || x + 1 == grid.dims[0]
            || y == 0
            || y + 1 == grid.dims[1]
            || z == 0
            || z + 1 == grid.dims[2]
        {
            reference_mark(state, &mut stack, index);
        }
    }
    while let Some(index) = stack.pop() {
        let nx = grid.dims[0];
        let plane = nx * grid.dims[1];
        let row = index % plane;
        let x = row % nx;
        if x > 0 {
            reference_mark(state, &mut stack, index - 1);
        }
        if x + 1 < nx {
            reference_mark(state, &mut stack, index + 1);
        }
        if row >= nx {
            reference_mark(state, &mut stack, index - nx);
        }
        if row + nx < plane {
            reference_mark(state, &mut stack, index + nx);
        }
        if index >= plane {
            reference_mark(state, &mut stack, index - plane);
        }
        if index + plane < state.len() {
            reference_mark(state, &mut stack, index + plane);
        }
    }
}

fn reference_mark(state: &mut [u8], stack: &mut Vec<usize>, index: usize) {
    if state[index] == EMPTY {
        state[index] = EXTERIOR;
        stack.push(index);
    }
}
