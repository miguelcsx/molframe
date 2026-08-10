use super::*;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

#[test]
fn the_centre_of_a_symmetric_pair_is_the_point_between_them() {
    let pair = [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let centre = centroid(&pair);
    assert!(centre.is_some_and(|c| close(c[0], 0.0) && close(c[1], 0.0)));
}

#[test]
fn the_centre_of_nothing_is_not_the_origin() {
    assert_eq!(centroid(&[]), None);
    assert_eq!(radius_of_gyration(&[], &[]), None);
    assert_eq!(inertia_tensor(&[], &[]), None);
}

#[test]
fn mass_moves_the_centre_towards_the_heavier_atom() {
    let pair = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
    let Some(geometric) = centroid(&pair) else {
        panic!("expected a centre")
    };
    let Some(weighted) = centre_of_mass(&pair, &[1.0, 9.0]) else {
        panic!("expected a centre")
    };
    assert!(close(geometric[0], 5.0));
    assert!(close(weighted[0], 9.0), "got {}", weighted[0]);
}

#[test]
fn masses_shorter_than_the_positions_are_refused_rather_than_padded() {
    let three = [[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    assert_eq!(centre_of_mass(&three, &[1.0, 1.0]), None);
}

#[test]
fn the_radius_of_gyration_of_a_symmetric_pair_is_its_half_separation() {
    let pair = [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    assert!(radius_of_gyration(&pair, &[]).is_some_and(|rg| close(rg, 1.0)));
}

#[test]
fn moving_a_set_does_not_change_its_radius_of_gyration() {
    let here = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
    let there = [[100.0, 50.0, 7.0], [101.0, 50.0, 7.0], [100.0, 52.0, 7.0]];
    let (Some(first), Some(second)) = (
        radius_of_gyration(&here, &[]),
        radius_of_gyration(&there, &[]),
    ) else {
        panic!("expected both")
    };
    assert!((first - second).abs() < 1e-6, "{first} vs {second}");
}

#[test]
fn the_inertia_tensor_is_symmetric_by_construction() {
    let points = [
        [1.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.0, 0.0, 3.0],
        [1.0, 1.0, 1.0],
    ];
    let Some(tensor) = inertia_tensor(&points, &[]) else {
        panic!("expected a tensor")
    };
    for (row, source) in tensor.iter().enumerate() {
        for (column, value) in source.iter().enumerate() {
            assert!(close(*value, tensor[column][row]), "({row}, {column})");
        }
    }
}

#[test]
fn the_long_axis_of_an_elongated_set_points_along_its_length() {
    let rod: Vec<[f32; 3]> = (-10_i16..=10)
        .map(|coordinate| [f32::from(coordinate), 0.0, 0.0])
        .collect();
    let Ok(Some(axes)) = gyration_axes(&rod) else {
        panic!("expected axes")
    };
    let long = axes.dominant();
    assert!(
        long[0].abs() > 0.999,
        "the long axis should be x, got {long:?}"
    );
}

#[test]
fn a_line_is_maximally_aspherical_and_a_symmetric_cloud_is_not() {
    let rod: Vec<[f32; 3]> = (-10_i16..=10)
        .map(|coordinate| [f32::from(coordinate), 0.0, 0.0])
        .collect();
    let Ok(Some(elongated)) = asphericity(&rod) else {
        panic!("expected a value")
    };
    assert!(elongated > 0.9, "a rod should be near one, got {elongated}");

    let cube = [
        [-1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [1.0, -1.0, 1.0],
        [-1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0],
    ];
    let Ok(Some(round)) = asphericity(&cube) else {
        panic!("expected a value")
    };
    assert!(
        round.abs() < 1e-9,
        "a cube should be near zero, got {round}"
    );
}

#[test]
fn principal_axes_come_back_ordered_by_their_moments() {
    let points = [
        [-3.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    let Ok(Some(axes)) = principal_axes(&points, &[]) else {
        panic!("expected axes")
    };
    assert!(axes.values[0] >= axes.values[1]);
    assert!(axes.values[1] >= axes.values[2]);
}
