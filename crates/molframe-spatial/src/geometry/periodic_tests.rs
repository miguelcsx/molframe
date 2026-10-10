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
    assert_eq!(
        periodic
            .minimum_image([0.2, 0.0, 0.0], [9.8, 0.0, 0.0])
            .lattice_shift,
        [1, 0, 0]
    );
}

#[test]
fn shortest_interpolation_crosses_the_periodic_boundary() {
    let periodic = PeriodicBox::from_cell(UnitCell {
        lengths: [10.0; 3],
        angles: [90.0; 3],
    });
    let Ok(periodic) = periodic else {
        panic!("valid cell rejected");
    };
    let middle = periodic.interpolate([9.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5);
    assert!(middle.is_some_and(|point| point[0].abs() < 1e-6));
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

#[test]
fn wrapping_uses_fractional_primary_cell_for_a_skewed_box() {
    let periodic = PeriodicBox::from_cell(UnitCell {
        lengths: [8.0, 9.0, 10.0],
        angles: [70.0, 80.0, 65.0],
    })
    .unwrap_or_else(|error| panic!("cell failed: {error}"));
    let position = periodic.cartesian([1.2, -0.25, 2.75]);
    let wrapped = periodic.fractional(periodic.wrap(position));
    assert!(
        wrapped
            .iter()
            .zip([0.2, 0.75, 0.75])
            .all(|(value, expected)| (value - expected).abs() < 1.0e-6)
    );
}

fn exhaustive_min_squared(periodic: &PeriodicBox, fractional: [f64; 3], range: i32) -> f64 {
    let mut best = f64::INFINITY;
    for i in -range..=range {
        for j in -range..=range {
            for k in -range..=range {
                let image = [
                    fractional[0] - f64::from(i),
                    fractional[1] - f64::from(j),
                    fractional[2] - f64::from(k),
                ];
                let cart = multiply(periodic.basis, image);
                best = best.min(squared(cart));
            }
        }
    }
    best
}

fn cell(lengths: [f64; 3], angles: [f64; 3]) -> Option<PeriodicBox> {
    PeriodicBox::from_cell(UnitCell { lengths, angles }).ok()
}

#[test]
fn a_highly_skewed_cell_finds_the_true_minimum_image() {
    let Some(periodic) = cell([10.0, 1.0, 10.0], [90.0, 90.0, 1.0]) else {
        panic!("valid cell rejected");
    };
    let fractional = [0.4, 0.4, 0.0];
    let left = [0.0_f64; 3];
    let right = multiply(periodic.basis, fractional);
    let found = periodic.displacement_f64(left, right);
    let expected = exhaustive_min_squared(&periodic, fractional, 12);
    assert!((squared(found) - expected).abs() < 1e-9, "{found:?}");
    assert!(expected.sqrt() < 0.41);
}

#[test]
fn reported_lattice_shift_reproduces_the_displacement() {
    let Some(periodic) = cell([10.0, 1.0, 10.0], [90.0, 90.0, 1.0]) else {
        panic!("valid cell rejected");
    };
    let right = multiply(periodic.basis, [0.4, 0.4, 0.0]).map(f64_f32);
    let image = periodic.minimum_image([0.0; 3], right);
    let shift = multiply(
        periodic.basis,
        image.lattice_shift.map(|v| i64_f64(v).unwrap_or(0.0)),
    );
    for axis in 0..3 {
        let expected = f64::from(right[axis]) - shift[axis];
        assert!((expected - f64::from(image.displacement[axis])).abs() < 1e-4);
    }
}

mod property {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn matches_exhaustive_search_for_random_cells(
            a in 1.0_f64..20.0, b in 1.0_f64..20.0, c in 1.0_f64..20.0,
            alpha in 1.0_f64..179.0, beta in 1.0_f64..179.0, gamma in 1.0_f64..179.0,
            f0 in -3.0_f64..3.0, f1 in -3.0_f64..3.0, f2 in -3.0_f64..3.0,
        ) {
            let Some(periodic) = cell([a, b, c], [alpha, beta, gamma]) else {
                return Ok(());
            };
            let fractional = [f0, f1, f2];
            let right = multiply(periodic.basis, fractional);
            let found = squared(periodic.displacement_f64([0.0; 3], right));
            let expected = exhaustive_min_squared(&periodic, fractional, 10);
            prop_assert!(found <= expected + 1e-9 * (1.0 + expected), "{found} vs {expected}");
        }
    }
}
