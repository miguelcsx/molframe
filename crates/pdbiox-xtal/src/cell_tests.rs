use super::*;

#[test]
fn orthorhombic_axes_scale_independently() {
    let transform = CellTransform::new(&UnitCell {
        lengths: [10.0, 20.0, 30.0],
        angles: [90.0; 3],
    })
    .expect("cell is valid");
    let cartesian = transform.to_cartesian([0.5, 0.25, 0.1]);
    for (actual, expected) in cartesian.into_iter().zip([5.0, 5.0, 3.0]) {
        assert!((actual - expected).abs() < 1.0e-12);
    }
}

#[test]
fn triclinic_round_trip_recovers_fractional_coordinates() {
    let transform = CellTransform::new(&UnitCell {
        lengths: [43.1, 51.7, 62.3],
        angles: [73.0, 81.0, 67.0],
    })
    .expect("cell is valid");
    let expected = [0.125, -0.25, 1.75];
    let actual = transform.to_fractional(transform.to_cartesian(expected));
    for (actual, expected) in actual.into_iter().zip(expected) {
        assert!((actual - expected).abs() < 1.0e-12);
    }
}

#[test]
fn degenerate_cells_are_refused_before_division() {
    let result = CellTransform::new(&UnitCell {
        lengths: [1.0, 1.0, 0.0],
        angles: [90.0; 3],
    });
    assert_eq!(result.err().map(|error| error.code()), Some(Code::E5004));
}
