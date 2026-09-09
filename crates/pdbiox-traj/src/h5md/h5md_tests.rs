use pdbiox_core::structure::UnitCell;

use crate::Timestep;

use super::{
    H5mdError, H5mdMetadata, H5mdOptions, H5mdUnitSystem, parse_h5md,
    parse_h5md_record_with_options, parse_h5md_with_options, write_h5md, write_h5md_with_metadata,
};

fn frames() -> Vec<Timestep> {
    vec![
        Timestep {
            frame: 10,
            time: Some(2.0),
            positions: vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            velocities: Some(vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]),
            forces: Some(vec![[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]]),
            cell: Some(UnitCell {
                lengths: [10.0, 11.0, 12.0],
                angles: [90.0, 90.0, 90.0],
            }),
            ..Timestep::default()
        },
        Timestep {
            frame: 20,
            time: Some(2.5),
            positions: vec![[1.5, 2.5, 3.5], [4.5, 5.5, 6.5]],
            velocities: Some(vec![[0.2, 0.3, 0.4], [0.5, 0.6, 0.7]]),
            forces: Some(vec![[11.0, 21.0, 31.0], [41.0, 51.0, 61.0]]),
            cell: Some(UnitCell {
                lengths: [10.5, 11.5, 12.5],
                angles: [90.0, 90.0, 90.0],
            }),
            ..Timestep::default()
        },
    ]
}

#[test]
fn roundtrips_canonical_h5md_with_all_streams() {
    let bytes = write_h5md(&frames(), &H5mdOptions::default()).expect("encode H5MD");
    assert_eq!(&bytes[..8], b"\x89HDF\r\n\x1a\n");
    let decoded = parse_h5md(&bytes).expect("decode H5MD");
    assert_eq!(decoded.len(), 2);
    assert_vector_close(decoded[0].positions[1], [4.0, 5.0, 6.0]);
    assert_close(
        decoded[1].velocities.as_ref().expect("velocities")[1][2],
        0.7,
    );
    assert_close(decoded[1].forces.as_ref().expect("forces")[0][0], 11.0);
    assert_eq!(
        decoded.iter().map(|frame| frame.frame).collect::<Vec<_>>(),
        [10, 20]
    );
    assert_eq!(decoded[1].time, Some(2.5));
    assert_eq!(decoded[1].dt, Some(0.5));
    let cell = decoded[1].cell.expect("cell");
    assert!((cell.lengths[2] - 12.5).abs() < 1.0e-10);
}

#[test]
fn applies_only_explicit_noncanonical_units() {
    let options = H5mdOptions {
        particle_group: "solvent".into(),
        units: H5mdUnitSystem {
            length_unit: "nm".into(),
            length_to_angstrom: 10.0,
            time_unit: "fs".into(),
            time_to_picosecond: 0.001,
            velocity_unit: "nm/fs".into(),
            velocity_to_angstrom_per_picosecond: 10_000.0,
            force_unit: "kJ/mol/nm".into(),
            force_to_kilojoule_per_mole_angstrom: 0.1,
        },
    };
    let bytes = write_h5md(&frames(), &options).expect("encode scaled H5MD");
    assert!(matches!(
        parse_h5md(&bytes),
        Err(H5mdError::InvalidUnits { .. })
    ));
    let decoded = parse_h5md_with_options(&bytes, &options).expect("decode scaled H5MD");
    assert_vector_close(decoded[1].positions[1], [4.5, 5.5, 6.5]);
    assert_eq!(decoded[0].time, Some(2.0));
    assert_close(decoded[0].forces.as_ref().expect("forces")[1][1], 50.0);
}

fn assert_close(actual: f32, expected: f32) {
    assert!((actual - expected).abs() < 1.0e-5);
}

fn assert_vector_close(actual: [f32; 3], expected: [f32; 3]) {
    actual
        .into_iter()
        .zip(expected)
        .for_each(|(actual, expected)| assert_close(actual, expected));
}

#[test]
fn rejects_invalid_unit_policy_and_mixed_schema() {
    let mut options = H5mdOptions::default();
    options.units.length_to_angstrom = 0.0;
    assert_eq!(
        write_h5md(&frames(), &options),
        Err(H5mdError::InvalidValue)
    );

    let mut input = frames();
    input[1].cell = None;
    assert_eq!(
        write_h5md(&input, &H5mdOptions::default()),
        Err(H5mdError::InconsistentFrames)
    );
}

#[test]
fn reads_external_nanometre_h5md_fixture_when_configured() {
    let Ok(path) = std::env::var("PDBIOX_H5MD_NANOMETRE_FIXTURE") else {
        return;
    };
    let Ok(particle_group) = std::env::var("PDBIOX_H5MD_PARTICLE_GROUP") else {
        return;
    };
    let bytes = std::fs::read(path).expect("read external H5MD fixture");
    let options = H5mdOptions {
        particle_group,
        units: H5mdUnitSystem {
            length_unit: "nm".into(),
            length_to_angstrom: 10.0,
            time_unit: "ps".into(),
            time_to_picosecond: 1.0,
            velocity_unit: "nm ps-1".into(),
            velocity_to_angstrom_per_picosecond: 10.0,
            force_unit: "kJ mol-1 nm-1".into(),
            force_to_kilojoule_per_mole_angstrom: 0.1,
        },
    };
    let frames = parse_h5md_with_options(&bytes, &options).expect("parse external H5MD fixture");
    assert!(!frames.is_empty());
    assert!(!frames[0].positions.is_empty());
}

#[test]
fn declarative_dispatch_preserves_creator_metadata() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let source_path = directory.path().join("source.h5md");
    let rewritten_path = directory.path().join("rewritten.h5md");
    let metadata = H5mdMetadata {
        options: H5mdOptions::default(),
        creator_name: "source-engine".into(),
        creator_version: "4.2".into(),
    };
    std::fs::write(
        &source_path,
        write_h5md_with_metadata(&frames(), &metadata).expect("source encode"),
    )
    .expect("source file");
    let trajectory =
        crate::read_trajectory_materialized(&source_path, &crate::TrajectoryReadOptions::default())
            .expect("dispatch read");
    crate::write_trajectory(
        &rewritten_path,
        &trajectory,
        &crate::TrajectoryWriteOptions::default(),
    )
    .expect("dispatch rewrite");
    let rewritten = parse_h5md_record_with_options(
        &std::fs::read(&rewritten_path).expect("rewritten bytes"),
        &H5mdOptions::default(),
    )
    .expect("rewritten decode");
    assert_eq!(rewritten.metadata.creator_name, "source-engine");
    assert_eq!(rewritten.metadata.creator_version, "4.2");
}
