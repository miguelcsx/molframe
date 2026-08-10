use super::{AmberRestartLayout, parse_amber_ascii_trajectory, parse_amber_restart};
use crate::{parse_amber_restart_record, write_amber_restart};
use std::fmt::Write;

fn fixed(values: &[f32], width: usize, per_line: usize) -> String {
    let mut output = String::new();
    for (index, value) in values.iter().enumerate() {
        if width == 8 {
            let _ = write!(output, "{value:>8.3}");
        } else {
            let _ = write!(output, "{value:>12.7}");
        }
        if (index + 1) % per_line == 0 {
            output.push('\n');
        }
    }
    if !output.ends_with('\n') {
        output.push('\n');
    }
    output
}

#[test]
fn restart_reads_velocities_time_and_six_value_cell() {
    let mut source = String::from("restart\n     2  1.5000000\n");
    let values = [
        0.0, 1.0, 2.0, 3.0, 4.0, 5.0, // coordinates
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, // velocities
        10.0, 11.0, 12.0, 80.0, 90.0, 100.0,
    ];
    source.push_str(&fixed(&values, 12, 6));
    let frame = parse_amber_restart(&source, AmberRestartLayout::Auto)
        .unwrap_or_else(|error| panic!("restart failed: {error}"));
    assert_eq!(frame.time, Some(1.5));
    assert_eq!(frame.positions.len(), 2);
    assert!(
        frame
            .velocities
            .as_ref()
            .is_some_and(|velocities| { (velocities[0][0] - 2.0455).abs() < 1.0e-4 })
    );
    let Some(cell) = frame.cell else {
        panic!("cell expected");
    };
    assert!((cell.angles[0] - 80.0).abs() < f64::EPSILON);
}

#[test]
fn one_atom_tail_requires_an_explicit_layout() {
    let source = format!(
        "inpcrd\n1\n{}",
        fixed(&[0.0, 0.0, 0.0, 10.0, 10.0, 10.0], 12, 6)
    );
    assert!(parse_amber_restart(&source, AmberRestartLayout::Auto).is_err());
    let frame = parse_amber_restart(&source, AmberRestartLayout::CoordinatesBox3)
        .unwrap_or_else(|error| panic!("explicit layout failed: {error}"));
    assert!(frame.cell.is_some());
}

#[test]
fn mdcrd_uses_explicit_atom_count_and_periodic_tail() {
    let mut source = String::from("trajectory\n");
    let values = [
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 10.0, 10.0, 10.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 11.0, 11.0,
        11.0,
    ];
    source.push_str(&fixed(&values, 8, 10));
    let frames = parse_amber_ascii_trajectory(&source, 2, true)
        .unwrap_or_else(|error| panic!("trajectory failed: {error}"));
    assert_eq!(frames.len(), 2);
    assert!(
        frames[1].positions[1]
            .iter()
            .zip([1.0, 1.0, 0.0])
            .all(|(observed, expected)| (observed - expected).abs() < f32::EPSILON)
    );
    assert!(frames.iter().all(|frame| frame.cell.is_some()));
}

#[test]
fn restart_roundtrips_title_units_velocities_and_cell() {
    let mut source = String::from("restart title\n     2  1.5000000\n");
    source.push_str(&fixed(
        &[
            0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 10.0, 11.0, 12.0,
        ],
        12,
        6,
    ));
    let record = parse_amber_restart_record(&source, AmberRestartLayout::CoordinatesVelocitiesBox3)
        .unwrap_or_else(|error| panic!("restart failed: {error}"));
    let encoded = write_amber_restart(&record)
        .unwrap_or_else(|error| panic!("restart write failed: {error}"));
    let reparsed =
        parse_amber_restart_record(&encoded, AmberRestartLayout::CoordinatesVelocitiesBox3)
            .unwrap_or_else(|error| panic!("restart reparse failed: {error}"));
    assert_eq!(reparsed.title, record.title);
    assert_eq!(reparsed.layout, record.layout);
    assert_eq!(reparsed.timestep.positions, record.timestep.positions);
    assert_eq!(reparsed.timestep.cell, record.timestep.cell);
    let original = record.timestep.velocities.as_ref().expect("velocities");
    let observed = reparsed.timestep.velocities.as_ref().expect("velocities");
    assert!(
        original
            .iter()
            .flatten()
            .zip(observed.iter().flatten())
            .all(|(left, right)| (left - right).abs() < 1.0e-5)
    );
}
