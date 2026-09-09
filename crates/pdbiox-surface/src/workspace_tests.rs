use super::{ensure_cavity_grid, ensure_ses_grid, ensure_ses_output};
use crate::cavity::Grid;
use crate::{
    SasaError, SurfaceGridOptions, cavities_with_options, solvent_excluded_surface_with_options,
};

/// The allowance the named standard profile grants.
const LIMIT: usize = SurfaceGridOptions::STANDARD_WORKSPACE_BYTES;

#[test]
fn the_standard_eight_million_cell_grid_fits_the_standard_allowance() {
    assert!(ensure_cavity_grid(1, 8_000_000, LIMIT).is_ok());
    assert!(ensure_ses_grid(1, 8_000_000, 200, LIMIT).is_ok());
}

#[test]
fn ses_output_capacity_is_counted_before_triangle_allocation() {
    let mut triangles = 1usize;
    while ensure_ses_output(1, 8_000_000, triangles, LIMIT).is_ok() {
        triangles = triangles.saturating_mul(2);
        assert_ne!(
            triangles,
            usize::MAX,
            "an unbounded surface output was accepted"
        );
    }

    let error = ensure_ses_output(1, 8_000_000, triangles, LIMIT)
        .expect_err("an oversized output must be rejected before allocation");
    assert!(matches!(
        error,
        SasaError::WorkspaceTooLarge { limit: LIMIT, .. }
    ));
}

#[test]
fn custom_cell_ceilings_cannot_bypass_the_absolute_byte_limit() {
    let error = ensure_cavity_grid(1, 20_000_000, LIMIT)
        .expect_err("the complete worst-case cavity output exceeds 500 MB");
    assert!(matches!(error, SasaError::WorkspaceTooLarge { .. }));
}

#[test]
#[ignore = "eight-million-cell RSS stress"]
fn standard_eight_million_cell_cavity_and_ses_complete() {
    let positions = [[0.0_f32; 3]];
    let radii = [9.8499_f32];
    let options = SurfaceGridOptions::standard(0.1);
    let grid = Grid::new(&positions, &[f64::from(radii[0])], options)
        .expect("the standard stress grid must fit");
    assert_eq!(grid.dims, [200, 200, 200]);
    assert_eq!(grid.cell_count(), 8_000_000);
    let cavities = cavities_with_options(&positions, &radii, 0.0, options)
        .expect("the standard cavity ceiling must complete");
    assert!(cavities.is_empty());

    let surface = solvent_excluded_surface_with_options(&positions, &radii, 0.0, options)
        .expect("the standard SES ceiling must complete");
    assert!(!surface.triangles.is_empty());
}
