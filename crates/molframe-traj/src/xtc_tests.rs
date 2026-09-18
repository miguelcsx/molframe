use molframe_core::structure::UnitCell;

use super::*;

fn frame(index: u16, offset: f32) -> Timestep {
    Timestep {
        frame: usize::from(index),
        time: Some(f64::from(index) * 0.5),
        positions: (0_i16..12)
            .map(|atom| [f32::from(atom) + offset, offset, -f32::from(atom)])
            .collect(),
        cell: Some(UnitCell {
            lengths: [20.0, 21.0, 22.0],
            angles: [90.0, 90.0, 90.0],
        }),
        ..Timestep::default()
    }
}

#[test]
fn round_trip_preserves_frames_and_metadata() {
    let input = vec![frame(4, 0.0), frame(5, 0.25)];
    let bytes = write_xtc(&input, XtcWriteOptions::default()).expect("encode XTC");
    let decoded = parse_xtc(&bytes).expect("decode XTC");

    assert_eq!(decoded.steps, [4, 5]);
    assert_eq!(decoded.precision, [1_000.0, 1_000.0]);
    assert_eq!(decoded.frames.len(), 2);
    assert_eq!(decoded.frames[0].dt, Some(0.5));
    for (actual, expected) in decoded.frames.iter().zip(input) {
        assert_eq!(actual.frame, expected.frame - 4);
        assert_eq!(actual.time, expected.time);
        for (left, right) in actual.positions.iter().zip(expected.positions) {
            for (value, wanted) in left.iter().zip(right) {
                assert!((value - wanted).abs() < 0.002);
            }
        }
        let cell = actual.cell.expect("cell");
        for (length, expected) in cell.lengths.iter().zip([20.0, 21.0, 22.0]) {
            assert!((length - expected).abs() < 0.000_01);
        }
    }
}

#[test]
fn refuses_truncation_and_shape_changes() {
    let input = vec![frame(0, 0.0)];
    let mut bytes = write_xtc(&input, XtcWriteOptions::default()).expect("encode XTC");
    bytes.pop();
    assert!(parse_xtc(&bytes).is_err());

    let mut mismatched = frame(1, 0.0);
    mismatched.positions.pop();
    assert!(write_xtc(&[frame(0, 0.0), mismatched], XtcWriteOptions::default()).is_err());
}
