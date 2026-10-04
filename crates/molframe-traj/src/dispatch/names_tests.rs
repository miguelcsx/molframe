use super::*;
use crate::Timestep;
use std::path::Path;

#[test]
fn every_named_format_round_trips_and_extension_inference_agrees() {
    for name in TrajectoryFormat::names() {
        let format: TrajectoryFormat = name.parse().expect("a listed name parses");
        assert_eq!(format.name(), name);
    }
    assert_eq!("netcdf".parse(), Ok(TrajectoryFormat::AmberNetcdf));
    assert_eq!("lammps_dump".parse(), Ok(TrajectoryFormat::LammpsDump));
    assert_eq!(
        TrajectoryFormat::infer(Path::new("run.xtc")),
        Some("xtc".parse().expect("xtc parses"))
    );
    assert!("mdcrd".parse::<TrajectoryFormat>().is_err());
}

#[test]
fn dense_positions_are_one_frame_major_block_with_missing_times_marked() {
    let frame = |frame: usize, time: Option<f64>, x: f32, atoms: usize| Timestep {
        frame,
        time,
        dt: None,
        positions: vec![[x, x + 1.0, x + 2.0]; atoms],
        velocities: None,
        forces: None,
        cell: None,
        data: std::collections::BTreeMap::new(),
    };
    let mut data = TrajectoryData {
        format: TrajectoryFormat::Xyz,
        frames: vec![frame(0, Some(0.5), 1.0, 2), frame(1, None, 4.0, 2)],
        metadata: TrajectoryMetadata {
            steps: None,
            format: FormatMetadata::None,
        },
    };
    let dense = data.dense_positions().expect("rectangular");
    assert_eq!((dense.frames, dense.atoms), (2, 2));
    assert_eq!(dense.coordinates.len(), 12);
    assert_eq!(&dense.coordinates[6..9], &[4.0, 5.0, 6.0]);
    assert!((dense.times[0] - 0.5).abs() < f64::EPSILON);
    assert!(dense.times[1].is_nan());
    data.frames.push(frame(2, None, 0.0, 3));
    assert!(matches!(
        data.dense_positions(),
        Err(crate::TrajectoryError::AtomCountMismatch {
            expected: 2,
            found: 3
        })
    ));
}
