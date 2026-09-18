use super::{parse_gro_records, write_gro};

/// Builds a GRO atom line in the fixed `%5d%-5s%5s%5d%8.3f%8.3f%8.3f` layout.
fn atom_line(resnum: i32, resname: &str, atom: &str, atomnum: i32, xyz: [f32; 3]) -> String {
    format!(
        "{resnum:>5}{resname:<5}{atom:>5}{atomnum:>5}{:>8.3}{:>8.3}{:>8.3}",
        xyz[0], xyz[1], xyz[2]
    )
}

#[test]
fn coordinates_are_read_and_converted_to_angstrom() {
    let text = format!(
        "water\n2\n{}\n{}\n   5.00000   5.00000   5.00000\n",
        atom_line(1, "WAT", "OW", 1, [1.0, 2.0, 3.0]),
        atom_line(1, "WAT", "HW", 2, [1.1, 2.0, 3.0]),
    );
    let frames = parse_gro_records(&text).expect("valid GRO");
    assert_eq!(frames.len(), 1);
    assert_eq!(frames[0].atoms.len(), 2);
    // 1.0 nm becomes 10 Å, and so on.
    let first = frames[0].atoms[0].position;
    assert!((first[0] - 10.0).abs() < 1e-3, "x {}", first[0]);
    assert!((first[1] - 20.0).abs() < 1e-3);
    assert!((first[2] - 30.0).abs() < 1e-3);
}

#[test]
fn two_stacked_frames_are_both_read() {
    let frame = format!(
        "t\n1\n{}\n   5.0   5.0   5.0\n",
        atom_line(1, "WAT", "OW", 1, [0.0, 0.0, 0.0]),
    );
    let frames = parse_gro_records(&format!("{frame}{frame}")).expect("valid");
    assert_eq!(frames.len(), 2);
}

#[test]
fn a_missing_box_line_fails_to_parse() {
    let text = format!("t\n1\n{}\n", atom_line(1, "WAT", "OW", 1, [0.0, 0.0, 0.0]));
    assert!(parse_gro_records(&text).is_err());
}

#[test]
fn written_frames_round_trip_coordinates_and_units() {
    let source = format!(
        "title\n1\n{}\n   5.0   5.0   5.0\n",
        atom_line(1, "WAT", "OW", 1, [1.0, 2.0, 3.0]),
    );
    let records = parse_gro_records(&source).expect("read complete GRO");
    let encoded = write_gro(&records).expect("write complete GRO");
    assert_eq!(parse_gro_records(&encoded).expect("reparse GRO"), records);
}
