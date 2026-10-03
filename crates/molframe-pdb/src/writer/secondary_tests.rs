use super::super::pdb::{PdbOptions, write};
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::{SecondaryStructure as Ss, Structure};
use std::fmt::Write as _;

#[test]
fn author_numbering_need_not_increase_along_a_secondary_range() {
    let structure = identity_fixture(
        &[("A", "20", ""), ("A", "10", ""), ("A", "15", "")],
        &[Ss::PolyProline, Ss::PolyProline, Ss::Unknown],
    );
    let output = write(
        &structure,
        &PdbOptions::new().namespace(super::super::pdb::PdbIdentifierNamespace::Auth),
    )
    .expect("author projection");
    let (roundtrip, _) = crate::read(
        &InputBuffer::from_bytes(output.into_bytes()),
        &ReadOptions::new(),
    )
    .expect("read output");
    assert_eq!(
        roundtrip.secondary_structure(),
        structure.secondary_structure()
    );
}

#[test]
fn selecting_out_a_residue_splits_secondary_records_without_extending_endpoints() {
    struct WithoutMiddle;
    impl molframe_core::io::Select for WithoutMiddle {
        fn accept_residue(&self, residue: molframe_core::index::ResidueIndex) -> bool {
            residue.get() != 1
        }
    }
    let structure = identity_fixture(
        &[("A", "10", "A"), ("A", "10", "B"), ("A", "10", "C")],
        &[Ss::PolyProline; 3],
    );
    let output = crate::write_selected(&structure, &PdbOptions::new(), &WithoutMiddle)
        .expect("selected writer");
    assert_eq!(
        output
            .lines()
            .filter(|line| line.starts_with("HELIX"))
            .count(),
        2
    );
    let (roundtrip, _) = crate::read(
        &InputBuffer::from_bytes(output.into_bytes()),
        &ReadOptions::new(),
    )
    .expect("read selected output");
    assert_eq!(roundtrip.secondary_structure(), &[Ss::PolyProline; 2]);
    assert_eq!(
        roundtrip
            .data()
            .residues()
            .filter_map(molframe_core::structure::ResidueRef::ins_code)
            .collect::<Vec<_>>(),
        ["A", "C"]
    );
}

#[test]
fn a_missing_secondary_component_refuses_before_streaming() {
    let source = write(&fixture(), &PdbOptions::new())
        .expect("fixture writer")
        .replace("GLY", "   ");
    let (structure, _) = crate::read(
        &InputBuffer::from_bytes(source.into_bytes()),
        &ReadOptions::new(),
    )
    .expect("blank component fixture");
    let mut output = Vec::new();
    assert!(crate::write_to(&structure, &PdbOptions::new(), &mut output).is_err());
    assert!(output.is_empty());
}

fn fixture() -> Structure {
    let mut source = String::new();
    for residue in 1..=7 {
        writeln!(source, "ATOM  {residue:>5}  CA  GLY A{residue:>4}    {:>8.3}   0.000   0.000  1.00 10.00           C", f64::from(residue) * 3.8).expect("write into String");
    }
    source.push_str(
        "END
",
    );
    let input = InputBuffer::from_bytes(source.into_bytes());
    let (raw, _) = crate::read(&input, &ReadOptions::new()).expect("fixture reads");
    let mut data = raw.data().clone();
    data.secondary_structure = vec![
        Ss::AlphaHelix,
        Ss::ThreeTenHelix,
        Ss::PiHelix,
        Ss::PolyProline,
        Ss::OtherHelix,
        Ss::Strand,
        Ss::Unknown,
    ]
    .into();
    Structure::new(data)
}

#[test]
fn generated_helix_classes_and_sheets_round_trip() {
    let structure = fixture();
    let text = write(&structure, &PdbOptions::new()).expect("write");
    let input = InputBuffer::from_bytes(text.into_bytes());
    let (roundtrip, _) = crate::read(&input, &ReadOptions::new()).expect("read output");
    assert_eq!(
        structure.secondary_structure(),
        roundtrip.secondary_structure()
    );
}

#[test]
fn renamed_chains_are_used_in_secondary_records_too() {
    let structure = fixture();
    let text = write(&structure, &PdbOptions::new().chain_map("A", "Z")).expect("write");
    for record in text.lines().filter(|line| line.starts_with("HELIX")) {
        assert_eq!(&record[19..20], "Z");
        assert_eq!(&record[31..32], "Z");
    }
    let input = InputBuffer::from_bytes(text.into_bytes());
    let (roundtrip, _) = crate::read(&input, &ReadOptions::new()).expect("read output");
    assert_eq!(
        structure.secondary_structure(),
        roundtrip.secondary_structure()
    );
}

fn identity_fixture(rows: &[(&str, &str, &str)], states: &[Ss]) -> Structure {
    let mut source = String::new();
    for (serial, &(chain, sequence, insertion)) in (1..).zip(rows) {
        writeln!(source, "ATOM  {serial:>5}  CA  GLY {chain:1}{sequence:>4}{insertion:1}      0.000   0.000   0.000  1.00 10.00           C").expect("String");
    }
    let (raw, _) = crate::read(
        &InputBuffer::from_bytes(source.into_bytes()),
        &ReadOptions::new(),
    )
    .expect("read fixture");
    let mut data = raw.data().clone();
    data.secondary_structure = states.to_vec().into();
    Structure::new(data)
}

#[test]
fn secondary_writer_round_trips_insertions_negative_and_both_hybrid36_cases() {
    let rows = [
        ("A", "-2", ""),
        ("A", "-1", ""),
        ("A", "10", ""),
        ("A", "10", "A"),
        ("A", "10", "B"),
        ("A", "10", "C"),
        ("A", "A000", ""),
        ("A", "A001", ""),
        ("A", "a000", ""),
        ("A", "a001", ""),
        ("B", "10", "A"),
    ];
    let states = [
        Ss::OtherHelix,
        Ss::OtherHelix,
        Ss::Unknown,
        Ss::PolyProline,
        Ss::PolyProline,
        Ss::Unknown,
        Ss::Strand,
        Ss::Strand,
        Ss::ThreeTenHelix,
        Ss::ThreeTenHelix,
        Ss::PiHelix,
    ];
    let structure = identity_fixture(&rows, &states);
    let options = PdbOptions::new().hybrid36(true);
    let mut output = Vec::new();
    crate::write_to(&structure, &options, &mut output).expect("stream real writer");
    let text = str::from_utf8(&output).expect("ASCII");
    let pp = text
        .lines()
        .find(|line| line.starts_with("HELIX") && &line[38..40] == "10")
        .expect("PPII record");
    assert_eq!(&pp[21..26], "  10A");
    assert_eq!(&pp[33..38], "  10B");
    let (roundtrip, _) = crate::read(&InputBuffer::from_bytes(output), &ReadOptions::new())
        .expect("read actual writer");
    assert_eq!(roundtrip.secondary_structure(), states);
    let ids: Vec<_> = roundtrip
        .data()
        .residues()
        .map(|residue| (residue.auth_seq_id(), residue.ins_code()))
        .collect();
    let original: Vec<_> = structure
        .data()
        .residues()
        .map(|residue| (residue.auth_seq_id(), residue.ins_code()))
        .collect();
    assert_eq!(ids, original);
}

#[test]
fn a_chain_map_collision_refuses_secondary_output_before_streaming() {
    let structure = identity_fixture(
        &[("A", "1", ""), ("B", "1", "")],
        &[Ss::PolyProline, Ss::OtherHelix],
    );
    let mut output = Vec::new();
    let result = crate::write_to(
        &structure,
        &PdbOptions::new().chain_map("B", "A"),
        &mut output,
    );
    assert!(result.is_err());
    assert!(output.is_empty());
}

#[test]
fn secondary_sequence_overflow_refuses_even_with_hybrid36_enabled() {
    let raw = fixture();
    let mut data = raw.data().clone();
    data.topology.residues = molframe_core::topology::ResidueTable::default();
    for residue in raw.data().residues() {
        let record = molframe_core::topology::ResidueRecord {
            label_comp_id: raw
                .data()
                .topology
                .residues
                .label_comp_id(residue.index())
                .expect("component"),
            auth_comp_id: molframe_core::optional::OptionalSymbol::NONE,
            label_seq_id: molframe_core::optional::OptionalI32::some(i32::MAX),
            auth_seq_id: molframe_core::optional::OptionalI32::NONE,
            ins_code: molframe_core::optional::OptionalSymbol::NONE,
            het: false,
        };
        data.topology
            .residues
            .push(record, residue.index().get()..residue.index().get() + 1)
            .expect("residue");
    }
    let mut output = Vec::new();
    assert!(
        crate::write_to(
            &Structure::new(data),
            &PdbOptions::new().hybrid36(true),
            &mut output
        )
        .is_err()
    );
    assert!(output.is_empty());
}

#[test]
fn unrepresentable_secondary_states_are_not_relabelled_as_helix_or_strand() {
    let structure = identity_fixture(
        &[
            ("A", "1", ""),
            ("A", "2", ""),
            ("A", "3", ""),
            ("A", "4", ""),
        ],
        &[Ss::Coil, Ss::Bend, Ss::BetaBridge, Ss::Turn],
    );
    let text = write(&structure, &PdbOptions::new()).expect("write atoms");
    assert!(
        !text
            .lines()
            .any(|line| line.starts_with("HELIX") || line.starts_with("SHEET"))
    );
}
