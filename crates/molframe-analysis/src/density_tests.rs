use super::{CartesianAxis, DensityGridSpec, LinearDensityOptions, density_map, linear_density};

#[test]
fn linear_density_preserves_weight_and_units() {
    let positions = [[0.25, 0.0, 0.0], [1.25, 0.0, 0.0], [5.0, 0.0, 0.0]];
    let Ok(bins) = linear_density(
        &positions,
        &[2.0, 4.0, 8.0],
        LinearDensityOptions {
            axis: CartesianAxis::X,
            minimum: 0.0,
            maximum: 2.0,
            bins: 2,
        },
    ) else {
        panic!("valid density profile");
    };
    assert!((bins[0].weight - 2.0).abs() < f64::EPSILON);
    assert!((bins[1].weight - 4.0).abs() < f64::EPSILON);
    assert!((bins[0].density - 2.0).abs() < f64::EPSILON);
}

#[test]
fn density_map_reports_outside_weight() {
    let spec = DensityGridSpec {
        origin: [0.0; 3],
        spacing: [1.0; 3],
        shape: [2; 3],
    };
    let Ok(grid) = density_map(&[[0.5, 0.5, 0.5], [2.0, 0.5, 0.5]], &[3.0, 7.0], spec) else {
        panic!("valid density grid");
    };
    assert!((grid.density[0] - 3.0).abs() < f64::EPSILON);
    assert!((grid.excluded_weight - 7.0).abs() < f64::EPSILON);
}

#[test]
fn density_map_normalises_by_the_voxel_volume() {
    let spec = DensityGridSpec {
        origin: [0.0; 3],
        spacing: [2.0; 3],
        shape: [1; 3],
    };
    let Ok(grid) = density_map(&[[1.0, 1.0, 1.0]], &[8.0], spec) else {
        panic!("valid density grid");
    };
    assert!((grid.density[0] - 1.0).abs() < f64::EPSILON);
}

#[test]
fn density_requires_explicit_valid_grid() {
    let result = density_map(
        &[],
        &[],
        DensityGridSpec {
            origin: [0.0; 3],
            spacing: [0.0; 3],
            shape: [1; 3],
        },
    );
    assert!(result.is_err());
}
