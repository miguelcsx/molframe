use super::*;

#[test]
fn borrowed_samples_preserve_endpoints_and_handle_extreme_finite_coordinates() {
    let start = [[f32::MAX, -f32::MAX, -0.0], [1.0, 2.0, 3.0]];
    let end = [[-f32::MAX, f32::MAX, 0.0], [3.0, 4.0, 5.0]];
    let sample = CoordinateSample::new(&start, &end, 0.5).expect("aligned columns");
    assert_eq!(sample.start.as_ptr(), start.as_ptr());
    assert_eq!(sample.end.as_ptr(), end.as_ptr());
    assert_eq!(sample.len(), 2);
    assert!(!sample.is_empty());
    assert_eq!(sample.position(0), Some([0.0; 3]));
    assert_eq!(sample.position(1), Some([2.0, 3.0, 4.0]));
    assert!(sample.position(2).is_none());
    assert_eq!(
        CoordinateSample::new(&start, &end, 0.0)
            .expect("start")
            .position(0)
            .expect("row")[2]
            .to_bits(),
        (-0.0f32).to_bits()
    );
    assert_eq!(
        CoordinateSample::new(&start, &end, 1.0)
            .expect("end")
            .position(0)
            .expect("row")[0]
            .to_bits(),
        (-f32::MAX).to_bits()
    );
}

#[test]
fn malformed_metadata_and_nonfinite_atom_rows_are_refused() {
    let start = [[0.0; 3]];
    let end = [[f32::NAN; 3]];
    assert!(CoordinateSample::new(&start, &[], 0.5).is_err());
    for weight in [-0.1, 1.1, f32::NAN, f32::INFINITY] {
        assert!(CoordinateSample::new(&start, &end, weight).is_err());
    }
    assert!(
        CoordinateSample::new(&start, &end, 0.5)
            .expect("metadata")
            .position(0)
            .is_none()
    );
    assert!(
        CoordinateSample::new(&[], &[], 0.5)
            .expect("empty columns")
            .is_empty()
    );
}
