use super::*;
use crate::{DcdWriteOptions, write_dcd};

#[test]
fn little_and_big_endian_charmm_cells_round_trip() {
    for endian in [DcdEndian::Little, DcdEndian::Big] {
        let frames = vec![Timestep {
            positions: vec![[1.0, 2.0, 3.0], [-4.0, 5.0, -6.0]],
            cell: Some(UnitCell {
                lengths: [10.0, 11.0, 12.0],
                angles: [70.0, 80.0, 75.0],
            }),
            ..Timestep::default()
        }];
        let bytes = write_dcd(
            &frames,
            &DcdWriteOptions {
                endian,
                ..DcdWriteOptions::default()
            },
        )
        .unwrap_or_else(|error| panic!("DCD write failed: {error}"));
        let parsed = parse_dcd(&bytes).unwrap_or_else(|error| panic!("DCD parse failed: {error}"));
        assert_eq!(parsed.frames[0].positions, frames[0].positions);
        let cell = parsed.frames[0]
            .cell
            .unwrap_or_else(|| panic!("cell absent"));
        assert!(
            cell.angles
                .iter()
                .zip([70.0, 80.0, 75.0])
                .all(|(left, right)| (left - right).abs() < 1.0e-10)
        );
    }
}

#[test]
fn fixed_atoms_are_carried_from_the_first_frame() {
    let endian = DcdEndian::Little;
    let mut bytes = Vec::new();
    let mut header = vec![0u8; 84];
    header[0..4].copy_from_slice(b"CORD");
    put(&mut header, 4, 2, endian);
    put(&mut header, 8, 0, endian);
    put(&mut header, 12, 1, endian);
    put(&mut header, 36, 1, endian);
    header[40..44].copy_from_slice(&1.0f32.to_le_bytes());
    put(&mut header, 80, 24, endian);
    add_record(&mut bytes, &header, endian);
    let mut title = vec![b' '; 84];
    put(&mut title, 0, 1, endian);
    add_record(&mut bytes, &title, endian);
    add_record(&mut bytes, &3i32.to_le_bytes(), endian);
    add_record(
        &mut bytes,
        &[2i32.to_le_bytes(), 3i32.to_le_bytes()].concat(),
        endian,
    );
    for axis in [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]] {
        add_record(&mut bytes, &floats(&axis, endian), endian);
    }
    for axis in [[20.0, 30.0], [50.0, 60.0], [80.0, 90.0]] {
        add_record(&mut bytes, &floats(&axis, endian), endian);
    }
    let parsed = parse_dcd(&bytes).unwrap_or_else(|error| panic!("fixed DCD failed: {error}"));
    assert_position(parsed.frames[1].positions[0], [1.0, 4.0, 7.0]);
    assert_position(parsed.frames[1].positions[2], [30.0, 60.0, 90.0]);
}

#[test]
fn xplor_lammps_style_header_without_cell_is_supported() {
    let endian = DcdEndian::Little;
    let mut bytes = Vec::new();
    let mut header = vec![0u8; 84];
    header[0..4].copy_from_slice(b"CORD");
    put(&mut header, 4, 1, endian);
    put(&mut header, 12, 1, endian);
    header[40..48].copy_from_slice(&1.0f64.to_le_bytes());
    add_record(&mut bytes, &header, endian);
    let mut title = vec![b' '; 84];
    put(&mut title, 0, 1, endian);
    add_record(&mut bytes, &title, endian);
    add_record(&mut bytes, &1i32.to_le_bytes(), endian);
    for value in [1.0f32, 2.0, 3.0] {
        add_record(&mut bytes, &value.to_le_bytes(), endian);
    }
    let parsed = parse_dcd(&bytes).unwrap_or_else(|error| panic!("X-PLOR DCD failed: {error}"));
    assert!(!parsed.header.charmm);
    assert_eq!(parsed.frames[0].positions, [[1.0, 2.0, 3.0]]);
}

fn put(bytes: &mut [u8], offset: usize, value: i32, endian: DcdEndian) {
    let value = match endian {
        DcdEndian::Little => value.to_le_bytes(),
        DcdEndian::Big => value.to_be_bytes(),
    };
    bytes[offset..offset + 4].copy_from_slice(&value);
}

fn add_record(output: &mut Vec<u8>, payload: &[u8], endian: DcdEndian) {
    let length = i32::try_from(payload.len())
        .unwrap_or_else(|_| panic!("test record exceeds the DCD marker range"));
    let marker = match endian {
        DcdEndian::Little => length.to_le_bytes(),
        DcdEndian::Big => length.to_be_bytes(),
    };
    output.extend(marker);
    output.extend(payload);
    output.extend(marker);
}

fn assert_position(actual: [f32; 3], expected: [f32; 3]) {
    assert!(
        actual
            .iter()
            .zip(expected)
            .all(|(left, right)| (*left - right).abs() < f32::EPSILON)
    );
}

fn floats(values: &[f32], endian: DcdEndian) -> Vec<u8> {
    values
        .iter()
        .flat_map(|value| match endian {
            DcdEndian::Little => value.to_le_bytes(),
            DcdEndian::Big => value.to_be_bytes(),
        })
        .collect()
}
