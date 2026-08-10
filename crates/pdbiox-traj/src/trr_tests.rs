use super::*;
use crate::{TrrWriteOptions, write_trr};
use pdbiox_core::structure::UnitCell;

#[test]
fn single_and_double_precision_round_trip_all_blocks() {
    for precision in [TrrPrecision::Single, TrrPrecision::Double] {
        let frames = vec![frame(0.0), frame(2.0)];
        let bytes = write_trr(&frames, TrrWriteOptions { precision })
            .unwrap_or_else(|error| panic!("TRR write failed: {error}"));
        let parsed = parse_trr(&bytes).unwrap_or_else(|error| panic!("TRR parse failed: {error}"));
        assert_eq!(parsed.precision, vec![precision, precision]);
        assert_eq!(parsed.frames.len(), 2);
        assert_close(parsed.frames[1].positions[1], [-4.0, 5.0, 6.0], 1.0e-5);
        assert_close(
            parsed.frames[0]
                .velocities
                .as_ref()
                .map(|values| values[0])
                .unwrap_or_default(),
            [0.1, 0.2, 0.3],
            1.0e-5,
        );
        assert_close(
            parsed.frames[0]
                .forces
                .as_ref()
                .map(|values| values[1])
                .unwrap_or_default(),
            [-7.0, 8.0, 9.0],
            1.0e-5,
        );
        assert!(parsed.frames[0].cell.is_some());
        assert!(parsed.frames.iter().all(|frame| frame.dt == Some(2.0)));
    }
}

#[test]
fn truncated_xdr_is_rejected() {
    let bytes = write_trr(&[frame(0.0)], TrrWriteOptions::default())
        .unwrap_or_else(|error| panic!("TRR write failed: {error}"));
    assert!(matches!(
        parse_trr(&bytes[..bytes.len() - 1]),
        Err(TrrError::Truncated { .. })
    ));
}

#[test]
fn external_trr_fixture_when_configured() {
    let Some(path) = std::env::var_os("PDBIOX_TRR_FIXTURE") else {
        return;
    };
    let bytes = std::fs::read(path).unwrap_or_else(|error| panic!("fixture read failed: {error}"));
    let parsed = parse_trr(&bytes).unwrap_or_else(|error| panic!("fixture parse failed: {error}"));
    assert!(!parsed.frames.is_empty());
    assert!(
        parsed
            .frames
            .iter()
            .all(|frame| !frame.positions.is_empty())
    );
}

fn frame(time: f64) -> Timestep {
    Timestep {
        time: Some(time),
        positions: vec![[1.0, 2.0, 3.0], [-4.0, 5.0, 6.0]],
        velocities: Some(vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]),
        forces: Some(vec![[1.0, 2.0, 3.0], [-7.0, 8.0, 9.0]]),
        cell: Some(UnitCell {
            lengths: [10.0, 11.0, 12.0],
            angles: [70.0, 80.0, 75.0],
        }),
        ..Timestep::default()
    }
}

fn assert_close(actual: [f32; 3], expected: [f32; 3], tolerance: f32) {
    assert!(
        actual
            .iter()
            .zip(expected)
            .all(|(left, right)| (*left - right).abs() < tolerance)
    );
}
