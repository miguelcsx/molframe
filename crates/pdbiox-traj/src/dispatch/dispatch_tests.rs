use crate::{Timestep, XtcWriteOptions};
use std::path::Path;

use super::*;

#[test]
fn xtc_path_round_trip_preserves_steps() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("frames.xtc");
    let data = TrajectoryData {
        format: TrajectoryFormat::Xtc,
        frames: vec![
            Timestep {
                frame: 0,
                time: Some(0.0),
                positions: vec![[0.0, 1.0, 2.0]],
                ..Timestep::default()
            },
            Timestep {
                frame: 1,
                time: Some(1.0),
                positions: vec![[3.0, 4.0, 5.0]],
                ..Timestep::default()
            },
        ],
        metadata: TrajectoryMetadata {
            steps: Some(vec![10, 20]),
            format: FormatMetadata::Xtc {
                precision: vec![1_000.0; 2],
            },
        },
    };
    write_trajectory(
        &path,
        &data,
        &TrajectoryWriteOptions {
            xtc: Some(XtcWriteOptions::default()),
            ..TrajectoryWriteOptions::default()
        },
    )
    .expect("write trajectory");
    let decoded =
        read_trajectory(&path, &TrajectoryReadOptions::default()).expect("read trajectory");
    assert_eq!(decoded.metadata.steps, Some(vec![10, 20]));
}

#[test]
fn xtc_rewrite_preserves_per_frame_precision() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("variable.xtc");
    let data = two_frame_data(
        TrajectoryFormat::Xtc,
        FormatMetadata::Xtc {
            precision: vec![100.0, 10_000.0],
        },
    );
    write_trajectory(&path, &data, &TrajectoryWriteOptions::default()).expect("write trajectory");
    let decoded =
        read_trajectory(&path, &TrajectoryReadOptions::default()).expect("read trajectory");
    assert_eq!(
        decoded.metadata.format,
        FormatMetadata::Xtc {
            precision: vec![100.0, 10_000.0]
        }
    );
}

#[test]
fn trr_rewrite_preserves_steps_and_per_frame_precision() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("variable.trr");
    let data = two_frame_data(
        TrajectoryFormat::Trr,
        FormatMetadata::Trr {
            precision: vec![crate::TrrPrecision::Single, crate::TrrPrecision::Double],
        },
    );
    write_trajectory(&path, &data, &TrajectoryWriteOptions::default()).expect("write trajectory");
    let decoded =
        read_trajectory(&path, &TrajectoryReadOptions::default()).expect("read trajectory");
    assert_eq!(decoded.metadata.steps, Some(vec![10, 20]));
    assert_eq!(decoded.metadata.format, data.metadata.format);
}

fn two_frame_data(format: TrajectoryFormat, metadata: FormatMetadata) -> TrajectoryData {
    let first = (0_u16..10)
        .map(|index| [f32::from(index), 1.0, 2.0])
        .collect();
    let second = (0_u16..10)
        .map(|index| [f32::from(index) + 3.0, 4.0, 5.0])
        .collect();
    TrajectoryData {
        format,
        frames: vec![
            Timestep {
                time: Some(0.0),
                positions: first,
                ..Timestep::default()
            },
            Timestep {
                time: Some(1.0),
                positions: second,
                ..Timestep::default()
            },
        ],
        metadata: TrajectoryMetadata {
            steps: Some(vec![10, 20]),
            format: metadata,
        },
    }
}

#[test]
fn frame_selection_reorders_repeats_and_aligns_metadata() {
    let data = two_frame_data(
        TrajectoryFormat::Xtc,
        FormatMetadata::Xtc {
            precision: vec![100.0, 200.0],
        },
    );
    let selected = data
        .select_frames(&[1, 0, 1])
        .expect("valid selection should succeed");
    assert_eq!(selected.frames.len(), 3);
    assert_eq!(selected.metadata.steps, Some(vec![20, 10, 20]));
    assert_eq!(
        selected.metadata.format,
        FormatMetadata::Xtc {
            precision: vec![200.0, 100.0, 200.0]
        }
    );
}

#[test]
fn frame_selection_rejects_out_of_range_indices() {
    let data = two_frame_data(TrajectoryFormat::Trr, FormatMetadata::None);
    assert!(matches!(
        data.select_frames(&[2]),
        Err(TrajectoryIoError::InvalidMetadata)
    ));
}

#[test]
fn frame_selection_rejects_preexisting_metadata_misalignment() {
    let mut data = two_frame_data(
        TrajectoryFormat::Xtc,
        FormatMetadata::Xtc {
            precision: vec![100.0],
        },
    );
    data.metadata.steps = Some(vec![10]);
    assert!(matches!(
        data.select_frames(&[0]),
        Err(TrajectoryIoError::InvalidMetadata)
    ));
}

#[test]
fn frame_selection_keeps_text_format_records_aligned() {
    let records =
        crate::parse_xyz("1\nfirst\nC 1 2 3\n1\nsecond\nO 4 5 6\n").expect("valid XYZ records");
    let data = TrajectoryData {
        format: TrajectoryFormat::Xyz,
        frames: records
            .iter()
            .enumerate()
            .map(|(frame, record)| Timestep {
                frame,
                positions: record.atoms.iter().map(|atom| atom.position).collect(),
                ..Timestep::default()
            })
            .collect(),
        metadata: TrajectoryMetadata {
            steps: None,
            format: FormatMetadata::Xyz(records),
        },
    };
    let selected = data.select_frames(&[1, 0, 1]).expect("select frames");
    let FormatMetadata::Xyz(records) = selected.metadata.format else {
        panic!("XYZ metadata expected");
    };
    assert_eq!(records[0].comment, "second");
    assert_eq!(records[1].comment, "first");
    assert_eq!(records[2].comment, "second");
}

#[test]
fn gsd_refuses_implicit_units() {
    let path = Path::new("trajectory.gsd");
    let error = read_trajectory(path, &TrajectoryReadOptions::default())
        .expect_err("GSD units must be explicit");
    assert!(matches!(error, TrajectoryIoError::MissingUnitOptions));
}

#[test]
fn namd_path_roundtrip_preserves_endian_and_requires_one_frame() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("snapshot.coor");
    let data = TrajectoryData {
        format: TrajectoryFormat::Namd,
        frames: vec![Timestep {
            positions: vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            ..Timestep::default()
        }],
        metadata: TrajectoryMetadata {
            steps: None,
            format: FormatMetadata::Namd(crate::NamdEndian::Big),
        },
    };
    write_trajectory(&path, &data, &TrajectoryWriteOptions::default()).expect("write NAMD");
    let decoded = read_trajectory(&path, &TrajectoryReadOptions::default()).expect("read NAMD");
    assert_eq!(decoded.frames, data.frames);
    assert_eq!(decoded.metadata.format, data.metadata.format);
}

#[test]
fn gro_dispatch_preserves_topology_velocities_and_box() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("frame.gro");
    let source = "title\n1\n    1WAT     OW    1   1.000   2.000   3.000  0.1000  0.2000  0.3000\n   5.00000   6.00000   7.00000\n";
    std::fs::write(&path, source).expect("fixture write");
    let data = read_trajectory(&path, &TrajectoryReadOptions::default()).expect("read GRO");
    write_trajectory(&path, &data, &TrajectoryWriteOptions::default()).expect("write GRO");
    let decoded = read_trajectory(&path, &TrajectoryReadOptions::default()).expect("re-read GRO");
    assert_eq!(decoded, data);
}

#[test]
fn xyz_dispatch_preserves_elements_and_comments() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("frames.xyz");
    std::fs::write(&path, "2\ncomment\nC 0 1 2\nO 3 4 5\n").expect("fixture write");
    let data = read_trajectory(&path, &TrajectoryReadOptions::default()).expect("read XYZ");
    write_trajectory(&path, &data, &TrajectoryWriteOptions::default()).expect("write XYZ");
    let decoded = read_trajectory(&path, &TrajectoryReadOptions::default()).expect("re-read XYZ");
    assert_eq!(decoded, data);
}

#[test]
fn txyz_dispatch_preserves_types_and_connectivity() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("molecule.txyz");
    std::fs::write(&path, "2 water\n1 O 0 0 0 10 2\n2 H 1 0 0 11 1\n").expect("fixture write");
    let data = read_trajectory(&path, &TrajectoryReadOptions::default()).expect("read TXYZ");
    write_trajectory(&path, &data, &TrajectoryWriteOptions::default()).expect("write TXYZ");
    assert_eq!(
        read_trajectory(&path, &TrajectoryReadOptions::default()).expect("re-read TXYZ"),
        data
    );
}

#[test]
fn aims_dispatch_preserves_species_and_lattice() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("geometry.in");
    let source =
        "lattice_vector 10 0 0\nlattice_vector 0 10 0\nlattice_vector 0 0 10\natom 1 2 3 C\n";
    std::fs::write(&path, source).expect("fixture write");
    let data = read_trajectory(&path, &TrajectoryReadOptions::default()).expect("read aims");
    write_trajectory(&path, &data, &TrajectoryWriteOptions::default()).expect("write aims");
    assert_eq!(
        read_trajectory(&path, &TrajectoryReadOptions::default()).expect("re-read aims"),
        data
    );
}

#[test]
fn dlpoly_config_dispatch_preserves_identity_and_controls() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("state.config");
    std::fs::write(&path, "title\n0 0 1\nC 7\n1.0 2.0 3.0\n").expect("fixture write");
    let data = read_trajectory(&path, &TrajectoryReadOptions::default()).expect("read CONFIG");
    write_trajectory(&path, &data, &TrajectoryWriteOptions::default()).expect("write CONFIG");
    assert_eq!(
        read_trajectory(&path, &TrajectoryReadOptions::default()).expect("re-read CONFIG"),
        data
    );
}

#[test]
fn dlpoly_history_dispatch_preserves_topology_steps_and_dynamics() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("trajectory.history");
    let source = "trajectory\n2 3 1\n\
timestep 10 1 2 3 0.5\n10 0 0\n0 10 0\n0 0 10\n\
C 7 12.011 -0.2\n1 2 3\n0.1 0.2 0.3\n100 200 300\n\
timestep 20 1 2 3 1.0\n10 0 0\n0 10 0\n0 0 10\n\
C 7 12.011 -0.2\n2 3 4\n0.2 0.3 0.4\n200 300 400\n";
    std::fs::write(&path, source).expect("fixture write");
    let data = read_trajectory(&path, &TrajectoryReadOptions::default()).expect("read HISTORY");
    write_trajectory(&path, &data, &TrajectoryWriteOptions::default()).expect("write HISTORY");
    assert_eq!(
        read_trajectory(&path, &TrajectoryReadOptions::default()).expect("re-read HISTORY"),
        data
    );
}

#[test]
fn charmm_card_dispatch_preserves_all_atom_identity_fields() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("coordinates.crd");
    let source = "* coordinate title\n* generated test\n    1\n    1    1 MOL  CA  -123.45600   2.50000  99.00000 SYS  1      0.00000\n";
    std::fs::write(&path, source).expect("fixture write");
    let data = read_trajectory(&path, &TrajectoryReadOptions::default()).expect("read CARD");
    write_trajectory(&path, &data, &TrajectoryWriteOptions::default()).expect("write CARD");
    assert_eq!(
        read_trajectory(&path, &TrajectoryReadOptions::default()).expect("re-read CARD"),
        data
    );
}

#[test]
fn documented_read_only_formats_refuse_writes_explicitly() {
    let directory = tempfile::tempdir().expect("temporary directory");
    for (name, format, metadata) in [
        (
            "calculation.gms",
            TrajectoryFormat::Gamess,
            FormatMetadata::None,
        ),
        (
            "trajectory.lammpsdump",
            TrajectoryFormat::LammpsDump,
            FormatMetadata::None,
        ),
        (
            "trajectory.trc",
            TrajectoryFormat::Gromos11,
            FormatMetadata::None,
        ),
    ] {
        let data = TrajectoryData {
            format,
            frames: Vec::new(),
            metadata: TrajectoryMetadata {
                steps: None,
                format: metadata,
            },
        };
        let error = write_trajectory(
            &directory.path().join(name),
            &data,
            &TrajectoryWriteOptions::default(),
        )
        .expect_err("read-only format must refuse writes");
        assert!(matches!(error, TrajectoryIoError::ReadOnlyFormat));
    }
}

#[test]
fn amber_ascii_dispatch_requires_and_uses_explicit_topology() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("trajectory.mdcrd");
    std::fs::write(
        &path,
        "AMBER trajectory\n   1.000   2.000   3.000   4.000   5.000   6.000\n  10.000  11.000  12.000\n",
    )
    .expect("fixture write");
    let error = read_trajectory(&path, &TrajectoryReadOptions::default())
        .expect_err("headerless trajectory needs its topology");
    assert!(matches!(error, TrajectoryIoError::MissingTopologyOptions));

    let data = read_trajectory(
        &path,
        &TrajectoryReadOptions {
            amber_ascii: Some(AmberAsciiReadOptions {
                atom_count: 2,
                periodic_box: true,
            }),
            ..TrajectoryReadOptions::default()
        },
    )
    .expect("read AMBER ASCII");
    assert_eq!(data.frames.len(), 1);
    assert_eq!(data.frames[0].positions.len(), 2);
    assert!(matches!(
        data.metadata.format,
        FormatMetadata::AmberAscii {
            atom_count: 2,
            periodic_box: true,
            ..
        }
    ));
}

#[test]
fn dms_dispatch_retains_complete_topology_in_format_metadata() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("system.dms");
    let source = crate::DmsSystem {
        topology: crate::DmsTopology {
            particles: vec![crate::DmsParticle {
                atomic_number: Some(6),
                residue_name: Some("LIG".into()),
                name: Some("C1".into()),
                ..crate::DmsParticle::default()
            }],
            bonds: Vec::new(),
        },
        frame: crate::DmsFrame {
            positions: vec![[1.0, 2.0, 3.0]],
            velocities: vec![Some([0.1, 0.2, 0.3])],
            cell: None,
        },
        version: None,
    };
    crate::write_dms(&path, &source).expect("write DMS fixture");
    let data = read_trajectory(&path, &TrajectoryReadOptions::default()).expect("dispatch DMS");
    assert_eq!(data.frames[0].positions, vec![[1.0, 2.0, 3.0]]);
    assert!(matches!(
        data.metadata.format,
        FormatMetadata::Dms(system) if *system == source
    ));
}

#[test]
fn every_documented_trajectory_suffix_has_declarative_dispatch() {
    for (name, expected) in [
        ("trajectory.ncdf", TrajectoryFormat::AmberNetcdf),
        ("trajectory.trj", TrajectoryFormat::AmberAscii),
        ("trajectory.crdbox", TrajectoryFormat::AmberAscii),
        ("coordinates.restrt", TrajectoryFormat::AmberRestart),
        ("coordinates.arc", TrajectoryFormat::Txyz),
        ("system.dms", TrajectoryFormat::Dms),
        ("geometry.in", TrajectoryFormat::Aims),
        ("trajectory.history", TrajectoryFormat::DlPolyHistory),
        ("calculation.gms", TrajectoryFormat::Gamess),
    ] {
        assert_eq!(TrajectoryFormat::infer(Path::new(name)), Some(expected));
    }
}

#[test]
fn tng_dispatch_preserves_lossy_codec_and_precision() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let source_path = directory.path().join("source.tng");
    let rewritten_path = directory.path().join("rewritten.tng");
    let frames = vec![
        Timestep {
            frame: 0,
            time: Some(0.0),
            dt: Some(1.0),
            positions: vec![[1.0, 2.0, 3.0], [2.0, 3.0, 4.0]],
            ..Timestep::default()
        },
        Timestep {
            frame: 5,
            time: Some(1.0),
            dt: Some(1.0),
            positions: vec![[4.0, 5.0, 6.0], [5.0, 6.0, 7.0]],
            ..Timestep::default()
        },
    ];
    crate::write_tng(
        &source_path,
        &frames,
        crate::TngWriteOptions {
            compression: crate::TngCompression::Lossy { precision: 400.0 },
            ..crate::TngWriteOptions::default()
        },
    )
    .expect("write source TNG");
    let data = read_trajectory(&source_path, &TrajectoryReadOptions::default()).expect("read TNG");
    write_trajectory(&rewritten_path, &data, &TrajectoryWriteOptions::default())
        .expect("rewrite TNG");
    let rewritten = read_trajectory(&rewritten_path, &TrajectoryReadOptions::default())
        .expect("read rewritten TNG");
    assert!(matches!(
        rewritten.metadata.format,
        FormatMetadata::Tng {
            compression: crate::TngCompression::Lossy { precision },
            ..
        } if (precision - 400.0).abs() < f64::EPSILON
    ));
}
