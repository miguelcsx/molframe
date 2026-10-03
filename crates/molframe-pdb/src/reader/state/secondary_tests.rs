use super::helix_class;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::{SecondarySource, SecondaryStructure as Ss};
use std::fmt::Write;

#[test]
fn sheet_file_records_keep_their_precedence_over_overlapping_helix_records() {
    let source = format!(
        "{}\n{}\n{}\n{}\n{}\n",
        identity_record(true, "A", ("10", "B"), ("10", "B"), 0),
        identity_record(false, "A", ("10", "A"), ("10", "B"), 10),
        identity_atom("A", "10", "", 1),
        identity_atom("A", "10", "A", 2),
        identity_atom("A", "10", "B", 3)
    );
    assert_eq!(
        read_identity(source).secondary_structure(),
        &[Ss::Unknown, Ss::PolyProline, Ss::Strand]
    );
}

#[test]
fn secondary_endpoints_do_not_connect_separate_chains_reusing_a_label() {
    let source = format!(
        "{}\n{}\nTER\n{}\n",
        identity_record(false, "A", ("1", ""), ("2", ""), 10),
        identity_atom("A", "1", "", 1),
        identity_atom("A", "2", "", 2)
    );
    assert_eq!(
        read_identity(source).secondary_structure(),
        &[Ss::Unknown, Ss::Unknown]
    );
}

#[test]
fn helix_classes_name_alpha_pi_and_three_ten() {
    assert_eq!(helix_class(Some(1)), Ss::AlphaHelix);
    assert_eq!(helix_class(Some(3)), Ss::PiHelix);
    assert_eq!(helix_class(Some(5)), Ss::ThreeTenHelix);
    assert_eq!(helix_class(Some(10)), Ss::PolyProline);
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

#[test]
fn excluding_helix_records_leaves_the_residues_unassigned() {
    use molframe_core::io::CategoryFilter;
    let input = InputBuffer::from_bytes(file_with_helix_classes().into_bytes());
    let options = ReadOptions::new().categories(CategoryFilter::except(["helix"]));
    let (structure, _) = crate::read(&input, &options).expect("fixture reads");
    assert_eq!(structure.atom_count(), 9);
    assert!(
        structure
            .data()
            .secondary_structure
            .iter()
            .all(|state| *state == Ss::Unknown)
    );
}

fn identity_record(
    sheet: bool,
    chain: &str,
    begin: (&str, &str),
    end: (&str, &str),
    class: u32,
) -> String {
    let mut line = [b' '; 80];
    put(&mut line, 1, 6, if sheet { "SHEET" } else { "HELIX" });
    let columns = if sheet {
        [22, 23, 26, 27, 33, 34, 37, 38]
    } else {
        [20, 22, 25, 26, 32, 34, 37, 38]
    };
    put(&mut line, columns[0], columns[0], chain);
    put(&mut line, columns[1], columns[2], begin.0);
    put(&mut line, columns[3], columns[3], begin.1);
    put(&mut line, columns[4], columns[4], chain);
    put(&mut line, columns[5], columns[6], end.0);
    put(&mut line, columns[7], columns[7], end.1);
    if !sheet {
        put(&mut line, 39, 40, &class.to_string());
    }
    String::from_utf8(line.to_vec()).expect("ASCII")
}

fn identity_atom(chain: &str, sequence: &str, insertion: &str, serial: u32) -> String {
    format!(
        "ATOM  {serial:>5}  CA  GLY {chain:1}{sequence:>4}{insertion:1}   {:>8.3}   0.000   0.000  1.00 10.00           C",
        f64::from(serial)
    )
}

fn read_identity(source: String) -> molframe_core::Structure {
    crate::read(
        &InputBuffer::from_bytes(source.into_bytes()),
        &ReadOptions::new(),
    )
    .expect("identity fixture reads")
    .0
}

#[test]
fn insertion_endpoints_bound_helix_and_sheet_without_touching_neighboring_insertions() {
    for (sheet, kind) in [(false, Ss::PolyProline), (true, Ss::Strand)] {
        let mut source = identity_record(sheet, "A", ("10", "A"), ("10", "B"), 10);
        source.push('\n');
        for (serial, insertion) in (1..).zip(["", "A", "B", "C"]) {
            writeln!(source, "{}", identity_atom("A", "10", insertion, serial)).expect("String");
        }
        writeln!(source, "{}", identity_atom("B", "10", "A", 5)).expect("String");
        let structure = read_identity(source);
        assert_eq!(
            structure.secondary_structure(),
            &[Ss::Unknown, kind, kind, Ss::Unknown, Ss::Unknown]
        );
        assert_eq!(
            &*structure.data().secondary_source,
            &[
                SecondarySource::None,
                SecondarySource::File,
                SecondarySource::File,
                SecondarySource::None,
                SecondarySource::None
            ]
        );
    }
}

#[test]
fn negative_and_hybrid36_secondary_endpoints_use_the_atom_number_decoder() {
    for sheet in [false, true] {
        for sequences in [
            ["-2", "-1", "0"],
            ["9999", "A000", "A001"],
            ["ZZZZ", "a000", "a001"],
        ] {
            let mut source = identity_record(sheet, "A", (sequences[1], ""), (sequences[2], ""), 7);
            source.push('\n');
            for (serial, sequence) in (1..).zip(sequences) {
                writeln!(source, "{}", identity_atom("A", sequence, "", serial)).expect("String");
            }
            let structure = read_identity(source);
            let kind = if sheet { Ss::Strand } else { Ss::OtherHelix };
            assert_eq!(structure.secondary_structure(), &[Ss::Unknown, kind, kind]);
        }
    }
}

#[test]
fn missing_insertion_endpoint_does_not_guess_a_numeric_range() {
    let source = format!(
        "{}
{}
{}
",
        identity_record(false, "A", ("10", "A"), ("10", "B"), 1),
        identity_atom("A", "10", "", 1),
        identity_atom("A", "10", "C", 2)
    );
    assert_eq!(
        read_identity(source).secondary_structure(),
        &[Ss::Unknown, Ss::Unknown]
    );
}

#[test]
fn a_blank_chain_is_a_real_secondary_identity() {
    let source = format!(
        "{}
{}
",
        identity_record(false, "", ("-1", ""), ("-1", ""), 10),
        identity_atom("", "-1", "", 1)
    );
    assert_eq!(
        read_identity(source).secondary_structure(),
        &[Ss::PolyProline]
    );
}

#[test]
fn later_model_records_do_not_overwrite_the_shared_first_model_secondary_column() {
    let atom = identity_atom("A", "1", "", 1);
    let source = format!(
        "MODEL        1
{}
{atom}
ENDMDL
MODEL        2
{}
{atom}
ENDMDL
",
        identity_record(false, "A", ("1", ""), ("1", ""), 10),
        identity_record(false, "A", ("1", ""), ("1", ""), 7)
    );
    let structure = read_identity(source);
    assert_eq!(structure.model_count(), 2);
    assert_eq!(structure.secondary_structure(), &[Ss::PolyProline]);
}

#[test]
fn ragged_model_secondary_records_are_local_but_file_headers_are_shared() {
    let source = format!(
        "{}
MODEL        1
{}
{}
{}
ENDMDL
MODEL        2
{}
{}
{}
{}
ENDMDL
",
        identity_record(true, "A", ("1", ""), ("1", ""), 0),
        identity_record(false, "A", ("2", ""), ("2", ""), 10),
        identity_atom("A", "1", "", 1),
        identity_atom("A", "2", "", 2),
        identity_record(false, "A", ("2", ""), ("2", ""), 7),
        identity_atom("A", "1", "", 1),
        identity_atom("A", "2", "", 2),
        identity_atom("A", "3", "", 3)
    );
    let structure = read_identity(source);
    let models = structure.ragged_models().expect("different topology");
    assert_eq!(
        models[0].secondary_structure(),
        &[Ss::Strand, Ss::PolyProline]
    );
    assert_eq!(
        models[1].secondary_structure(),
        &[Ss::Strand, Ss::OtherHelix, Ss::Unknown]
    );
}
