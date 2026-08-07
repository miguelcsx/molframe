use super::*;

#[test]
fn orthorhombic_minimum_image_wraps_across_a_boundary() {
    let periodic = PeriodicBox::from_cell(UnitCell {
        lengths: [10.0; 3],
        angles: [90.0; 3],
    });
    let periodic = match periodic {
        Ok(periodic) => periodic,
        Err(error) => panic!("cell failed: {error}"),
    };
    assert!((periodic.distance_squared([0.2, 0.0, 0.0], [9.8, 0.0, 0.0]) - 0.16).abs() < 1e-5);
}

#[test]
fn a_skewed_cell_checks_neighbouring_lattice_images() {
    let periodic = PeriodicBox::from_cell(UnitCell {
        lengths: [8.0, 9.0, 10.0],
        angles: [70.0, 80.0, 65.0],
    });
    let periodic = match periodic {
        Ok(periodic) => periodic,
        Err(error) => panic!("cell failed: {error}"),
    };
    let forward = periodic.distance_squared([0.1, 0.2, 0.3], [7.9, 0.2, 0.3]);
    let reverse = periodic.distance_squared([7.9, 0.2, 0.3], [0.1, 0.2, 0.3]);
    assert!((forward - reverse).abs() < f32::EPSILON);
    assert!(forward < 4.0);
}

#[test]
fn a_degenerate_cell_is_refused() {
    assert_eq!(
        PeriodicBox::from_cell(UnitCell {
            lengths: [10.0, 10.0, 0.0],
            angles: [90.0; 3],
        }),
        Err(SpatialError::InvalidCell)
    );
}
