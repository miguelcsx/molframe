use super::helix_class;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::{SecondarySource, SecondaryStructure as Ss};
use std::fmt::Write;

#[test]
fn helix_classes_name_alpha_pi_and_three_ten() {
    assert_eq!(helix_class(Some(1)), Ss::AlphaHelix);
    assert_eq!(helix_class(Some(3)), Ss::PiHelix);
    assert_eq!(helix_class(Some(5)), Ss::ThreeTenHelix);
    assert_eq!(helix_class(Some(7)), Ss::OtherHelix);
    assert_eq!(helix_class(None), Ss::OtherHelix);
}

/// Places `text` at 1-based inclusive PDB columns `from..=to`, right-aligned.
fn put(line: &mut [u8; 80], from: usize, to: usize, text: &str) {
    let width = to - from + 1;
    let padded = format!("{text:>width$}");
    line[from - 1..to].copy_from_slice(padded.as_bytes());
}

fn helix_record(serial: u32, first: u32, last: u32, class: u32) -> String {
    let mut line = [b' '; 80];
    put(&mut line, 1, 5, "HELIX");
    put(&mut line, 8, 10, &serial.to_string());
    put(&mut line, 12, 14, &serial.to_string());
    put(&mut line, 16, 18, "GLY");
    put(&mut line, 20, 20, "A");
    put(&mut line, 22, 25, &first.to_string());
    put(&mut line, 28, 30, "GLY");
    put(&mut line, 32, 32, "A");
    put(&mut line, 34, 37, &last.to_string());
    put(&mut line, 39, 40, &class.to_string());
    put(&mut line, 72, 76, &(last - first + 1).to_string());
    String::from_utf8(line.to_vec()).expect("ascii record")
}

/// One Cα per residue, so the HELIX class is all that names each range.
fn file_with_helix_classes() -> String {
    let mut text = String::new();
    for (serial, (first, last, class)) in (1..).zip([(1, 2, 1), (3, 4, 3), (5, 6, 5), (7, 8, 7)]) {
        writeln!(text, "{}", helix_record(serial, first, last, class)).expect("string write");
    }
    for residue in 1..=9_u32 {
        writeln!(
            text,
            "ATOM  {residue:>5}  CA  GLY A{residue:>4}    {:>8.3}   0.000   0.000  1.00 10.00           C",
            f64::from(residue) * 3.8
        )
        .expect("string write");
    }
    text.push_str("END\n");
    text
}

#[test]
fn a_helix_record_reaches_the_structure_with_its_class_and_a_file_source() {
    let input = InputBuffer::from_bytes(file_with_helix_classes().into_bytes());
    let (structure, _) = crate::read(&input, &ReadOptions::new()).expect("fixture reads");
    let states = structure.data().secondary_structure.to_vec();
    assert_eq!(
        states,
        [
            Ss::AlphaHelix,
            Ss::AlphaHelix,
            Ss::PiHelix,
            Ss::PiHelix,
            Ss::ThreeTenHelix,
            Ss::ThreeTenHelix,
            Ss::OtherHelix,
            Ss::OtherHelix,
            Ss::Unknown,
        ]
    );
    let sources = structure.data().secondary_source.to_vec();
    assert!(
        sources[..8]
            .iter()
            .all(|source| *source == SecondarySource::File)
    );
    assert_eq!(sources[8], SecondarySource::None);
}
