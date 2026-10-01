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

#[test]
fn reciprocal_vectors_are_dual_to_triclinic_direct_axes() {
    let transform = CellTransform::new(&UnitCell {
        lengths: [43.1, 51.7, 62.3],
        angles: [73.0, 81.0, 67.0],
    })
    .expect("valid cell");
    let hkl = [2, -3, 5];
    let reciprocal = transform.reciprocal_vector(hkl);
    for axis in 0..3 {
        let mut unit = [0.0; 3];
        unit[axis] = 1.0;
        let direct = transform.to_cartesian(unit);
        let dot: f64 = direct.into_iter().zip(reciprocal).map(|(a, b)| a * b).sum();
        assert!((dot - f64::from(hkl[axis])).abs() < 1e-12);
    }
    assert!(transform.d_spacing([0; 3]).is_infinite());
    assert_eq!(
        transform.reciprocal_spacing_squared([0; 3]).to_bits(),
        0.0_f64.to_bits()
    );
    assert!(
        (transform.d_spacing(hkl).powi(-2) - transform.reciprocal_spacing_squared(hkl)).abs()
            < 1e-14
    );
}
