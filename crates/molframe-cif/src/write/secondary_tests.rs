use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::{SecondaryStructure as Ss, Structure};
use std::fmt::Write as _;

fn fixture() -> Structure {
    let mut source = String::from(
        "data_secondary
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
",
    );
    for index in 1..=8 {
        writeln!(source, "ATOM {index} C CA GLY A {index} {index} 0 0 1")
            .expect("write into String");
    }
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
        Ss::Turn,
        Ss::Unknown,
    ]
    .into();
    Structure::new(data)
}

#[test]
fn canonical_secondary_records_round_trip_all_representable_states() {
    let structure = fixture();
    let text = crate::write_canonical_with_options(
        &structure,
        &crate::CifWriteOptions::new().with_block_id("secondary"),
    )
    .expect("write");
    assert!(text.contains("HELX_P C4 10"), "{text}");
    let input = InputBuffer::from_bytes(text.into_bytes());
    let (roundtrip, _) = crate::read(&input, &ReadOptions::new()).expect("read output");
    assert_eq!(
        structure.secondary_structure(),
        roundtrip.secondary_structure()
    );
}

#[test]
fn missing_secondary_sequence_refuses_before_streaming() {
    let input = InputBuffer::from_bytes(b"data_missing\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\nATOM 1 C CA GLY A . 0 0 0 1\n".to_vec());
    let (raw, _) =
        crate::read(&input, &ReadOptions::new()).expect("missing label sequence is readable");
    let mut data = raw.data().clone();
    data.secondary_structure = vec![Ss::PolyProline].into();
    let structure = Structure::new(data);
    let mut output = Vec::new();
    let result = crate::write_canonical_to(
        &structure,
        &crate::CifWriteOptions::new().with_block_id("secondary"),
        &mut output,
    );
    assert!(result.is_err());
    assert!(output.is_empty());
}

fn identity_fixture(rows: &str, states: &[Ss]) -> Structure {
    let source = format!(
        "data_identity\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.label_seq_id\n_atom_site.auth_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\n{rows}"
    );
    let (raw, _) = crate::read(
        &InputBuffer::from_bytes(source.into_bytes()),
        &ReadOptions::new(),
    )
    .expect("read identity fixture");
    let mut data = raw.data().clone();
    data.secondary_structure = states.to_vec().into();
    Structure::new(data)
}

#[test]
fn canonical_secondary_preserves_label_namespace_chain_order_and_author_insertions() {
    let states = [
        Ss::PolyProline,
        Ss::PolyProline,
        Ss::Unknown,
        Ss::OtherHelix,
        Ss::Strand,
    ];
    let structure = identity_fixture(
        "ATOM 1 C CA GLY A B 5 -2 A 0 0 0 1\nATOM 2 C CA GLY A B 1 -2 B 1 0 0 1\nATOM 3 C CA GLY A B 3 10000 . 2 0 0 1\nATOM 4 C CA GLY B A 5 -2 A 3 0 0 1\nATOM 5 C CA GLY B A 1 -2 B 4 0 0 1\n",
        &states,
    );
    let mut output = Vec::new();
    crate::write_canonical_to(
        &structure,
        &crate::CifWriteOptions::new().with_block_id("identity"),
        &mut output,
    )
    .expect("actual canonical writer");
    let (roundtrip, _) = crate::read(&InputBuffer::from_bytes(output), &ReadOptions::new())
        .expect("actual readback");
    assert_eq!(roundtrip.secondary_structure(), states);
    let ids: Vec<_> = roundtrip
        .data()
        .residues()
        .map(|r| (r.label_seq_id(), r.auth_seq_id(), r.ins_code()))
        .collect();
    let original: Vec<_> = structure
        .data()
        .residues()
        .map(|r| (r.label_seq_id(), r.auth_seq_id(), r.ins_code()))
        .collect();
    assert_eq!(ids, original);
}

#[test]
fn duplicate_secondary_label_identities_refuse_before_streaming() {
    let structure = identity_fixture(
        "ATOM 1 C CA GLY A A 1 10 A 0 0 0 1\nATOM 2 C CA GLY A A 1 10 B 1 0 0 1\n",
        &[Ss::PolyProline, Ss::OtherHelix],
    );
    let mut output = Vec::new();
    let result = crate::write_canonical_to(
        &structure,
        &crate::CifWriteOptions::new().with_block_id("identity"),
        &mut output,
    );
    assert!(result.is_err());
    assert!(output.is_empty());
}

#[test]
fn independent_model_secondary_annotations_refuse_before_streaming() {
    let first = fixture();
    let mut data = first.data().clone();
    data.coords = molframe_core::structure::CoordinateStore::Ragged {
        models: vec![first],
    };
    let structure = Structure::new(data);
    let mut output = Vec::new();
    let result = crate::write_canonical_to(
        &structure,
        &crate::CifWriteOptions::new().with_block_id("identity"),
        &mut output,
    );
    assert!(result.is_err());
    assert!(output.is_empty());
}

#[test]
fn coil_bend_and_isolated_bridge_are_not_projected_as_other_states() {
    let structure = identity_fixture(
        "ATOM 1 C CA GLY A A 1 1 . 0 0 0 1\nATOM 2 C CA GLY A A 2 2 . 1 0 0 1\nATOM 3 C CA GLY A A 3 3 . 2 0 0 1\n",
        &[Ss::Coil, Ss::Bend, Ss::BetaBridge],
    );
    let output = crate::write_canonical_with_options(
        &structure,
        &crate::CifWriteOptions::new().with_block_id("identity"),
    )
    .expect("atom output");
    assert!(!output.contains("_struct_conf."));
    assert!(!output.contains("_struct_sheet_range."));
}
