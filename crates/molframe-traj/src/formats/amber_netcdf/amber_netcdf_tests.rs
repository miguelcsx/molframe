use molframe_core::structure::UnitCell;

use crate::Timestep;

use super::{
    AmberNetcdfError, AmberNetcdfPrecision, AmberNetcdfWriteOptions, parse_amber_netcdf,
    parse_amber_netcdf_record, write_amber_netcdf,
};

fn frames() -> Vec<Timestep> {
    vec![
        Timestep {
            frame: 0,
            time: Some(1.0),
            positions: vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            velocities: Some(vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]),
            forces: Some(vec![[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]]),
            cell: Some(UnitCell {
                lengths: [10.0, 11.0, 12.0],
                angles: [90.0, 91.0, 92.0],
            }),
            ..Timestep::default()
        },
        Timestep {
            frame: 1,
            time: Some(1.5),
            positions: vec![[1.5, 2.5, 3.5], [4.5, 5.5, 6.5]],
            velocities: Some(vec![[0.2, 0.3, 0.4], [0.5, 0.6, 0.7]]),
            forces: Some(vec![[11.0, 21.0, 31.0], [41.0, 51.0, 61.0]]),
            cell: Some(UnitCell {
                lengths: [10.5, 11.5, 12.5],
                angles: [90.5, 91.5, 92.5],
            }),
            ..Timestep::default()
        },
    ]
}

#[test]
fn roundtrips_all_supported_streams_in_both_precisions() {
    for precision in [AmberNetcdfPrecision::Single, AmberNetcdfPrecision::Double] {
        let bytes = write_amber_netcdf(
            &frames(),
            AmberNetcdfWriteOptions {
                precision,
                ..AmberNetcdfWriteOptions::default()
            },
        )
        .expect("encode AMBER NetCDF");
        assert_eq!(&bytes[..4], b"CDF\x02");
        let decoded = parse_amber_netcdf(&bytes).expect("decode AMBER NetCDF");
        assert_eq!(decoded.len(), 2);
        assert_vector_close(decoded[1].positions[1], [4.5, 5.5, 6.5]);
        assert_eq!(decoded[1].dt, Some(0.5));
        assert_close(
            decoded[0].velocities.as_ref().expect("velocities")[0][1],
            0.2,
        );
        assert_close(decoded[1].forces.as_ref().expect("forces")[1][2], 61.0);
        let angles = decoded[1].cell.expect("cell").angles;
        assert!(
            angles
                .into_iter()
                .zip([90.5, 91.5, 92.5])
                .all(|(a, b)| (a - b).abs() < 1.0e-6)
        );
    }
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
fn rejects_inconsistent_optional_streams() {
    let mut input = frames();
    input[1].velocities = None;
    assert_eq!(
        write_amber_netcdf(&input, AmberNetcdfWriteOptions::default()),
        Err(AmberNetcdfError::InconsistentFrames)
    );
}

#[test]
fn rejects_missing_unit_metadata() {
    let mut builder = netcdf_writer::NcFileBuilder::new();
    builder
        .add_attribute(
            "Conventions",
            netcdf_writer::NcAttrValue::Chars("AMBER".into()),
        )
        .expect("convention");
    builder
        .add_attribute(
            "ConventionVersion",
            netcdf_writer::NcAttrValue::Chars("1.0".into()),
        )
        .expect("version");
    let frame = builder.add_dimension("frame", 1).expect("frame");
    let atom = builder.add_dimension("atom", 1).expect("atom");
    let spatial = builder.add_dimension("spatial", 3).expect("spatial");
    let coordinates = builder
        .add_variable::<f32>("coordinates", &[frame, atom, spatial])
        .expect("coordinates");
    builder
        .write_variable(coordinates, &[1.0_f32, 2.0, 3.0])
        .expect("values");
    let (_, bytes) = builder
        .to_vec(netcdf_writer::NcWriteOptions::classic())
        .expect("file");
    assert!(matches!(
        parse_amber_netcdf(&bytes),
        Err(AmberNetcdfError::InvalidUnits { .. })
    ));
}

#[test]
fn reads_external_amber_netcdf_fixture_when_configured() {
    let Ok(path) = std::env::var("MOLFRAME_AMBER_NETCDF_FIXTURE") else {
        return;
    };
    let bytes = std::fs::read(path).expect("read external AMBER NetCDF fixture");
    let frames = parse_amber_netcdf(&bytes).expect("parse external AMBER NetCDF fixture");
    assert!(!frames.is_empty());
    assert!(!frames[0].positions.is_empty());
}

#[test]
fn declarative_dispatch_preserves_precision_and_producer_metadata() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let source_path = directory.path().join("source.ncdf");
    let rewritten_path = directory.path().join("rewritten.ncdf");
    let source_frames = frames();
    let source_options = AmberNetcdfWriteOptions {
        precision: AmberNetcdfPrecision::Double,
        program: "source-engine".into(),
        program_version: "9.4".into(),
    };
    std::fs::write(
        &source_path,
        write_amber_netcdf(&source_frames, source_options.clone()).expect("source encode"),
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
    let rewritten =
        parse_amber_netcdf_record(&std::fs::read(&rewritten_path).expect("rewritten file bytes"))
            .expect("rewritten decode");
    assert_eq!(rewritten.metadata.precision, source_options.precision);
    assert_eq!(rewritten.metadata.program.as_deref(), Some("source-engine"));
    assert_eq!(rewritten.metadata.program_version.as_deref(), Some("9.4"));
}
