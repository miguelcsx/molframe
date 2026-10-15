use molframe_core::structure::UnitCell;

use super::*;

#[test]
fn round_trip_requires_and_applies_units() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("trajectory.gsd");
    let frames = vec![Timestep {
        frame: 17,
        positions: vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
        cell: Some(UnitCell {
            lengths: [20.0, 30.0, 40.0],
            angles: [90.0, 90.0, 90.0],
        }),
        ..Timestep::default()
    }];
    let options = GsdOptions::new(10.0).expect("explicit nm conversion");
    write_gsd(&path, &frames, options).expect("write GSD");
    let decoded = parse_gsd(&path, options).expect("read GSD");

    assert_eq!(decoded.steps, [17]);
    assert_eq!(decoded.frames[0].positions, frames[0].positions);
    let cell = decoded.frames[0].cell.expect("periodic cell");
    for (actual, expected) in cell.lengths.iter().zip([20.0, 30.0, 40.0]) {
        assert!((actual - expected).abs() < 0.000_01);
    }
}

#[test]
fn refuses_implicit_units() {
    assert_eq!(GsdOptions::new(0.0), Err(GsdError::InvalidUnits));
    assert_eq!(GsdOptions::new(f64::NAN), Err(GsdError::InvalidUnits));
}
