use super::write_mmtf;
use crate::mmtf::codec::{encode_f32, encode_i32, encode_strings};
use crate::mmtf::schema::{Entity, File, Group};
use molframe_core::{InputBuffer, ReadOptions};

#[test]
fn canonical_mmtf_round_trip_preserves_hierarchy_coordinates_and_bonds() {
    let input = InputBuffer::from_bytes(minimal_mmtf());
    let (structure, _) = crate::read_mmtf(&input, &ReadOptions::new()).expect("MMTF fixture");
    let encoded = write_mmtf(&structure).expect("MMTF write");
    let input = InputBuffer::from_bytes(encoded);
    let (decoded, findings) = crate::read_mmtf(&input, &ReadOptions::new()).expect("MMTF read");
    assert!(findings.is_empty());
    assert_eq!(decoded.atom_count(), 2);
    assert_eq!(decoded.residue_count(), 1);
    assert_eq!(decoded.chain_count(), 1);
    assert_eq!(decoded.data().bonds.len(), 1);
    assert_eq!(decoded.positions(), structure.positions());
}

#[test]
fn repeated_chemistry_occupies_one_group_dictionary_entry() {
    // Three glycines share one chemistry. MMTF stores that chemistry once and
    // refers to it by index, so the dictionary must not grow with the residue
    // count, and every residue must still resolve back to glycine.
    let input = InputBuffer::from_bytes(repeated_glycine_mmtf());
    let (structure, _) = crate::read_mmtf(&input, &ReadOptions::new()).expect("MMTF fixture");
    assert_eq!(structure.residue_count(), 3);

    let encoded = write_mmtf(&structure).expect("MMTF write");
    let written: File = rmp_serde::from_slice(&encoded).expect("written MessagePack");
    assert_eq!(
        written.group_list.len(),
        1,
        "one chemistry should occupy one dictionary entry"
    );
    assert_eq!(written.num_groups, 3, "three residues are still reported");

    let input = InputBuffer::from_bytes(encoded);
    let (decoded, findings) = crate::read_mmtf(&input, &ReadOptions::new()).expect("MMTF read");
    assert!(findings.is_empty());
    assert_eq!(decoded.residue_count(), 3);
    assert_eq!(decoded.atom_count(), 6);
    assert_eq!(decoded.positions(), structure.positions());
}

#[test]
fn writer_refuses_structure_without_explicit_mmtf_chemistry() {
    let source =
        b"ATOM      1  N   GLY A   1      11.104  13.207   9.301  1.00 20.00           N  \nEND\n";
    let input = InputBuffer::from_bytes(source.to_vec());
    let (structure, _) = crate::read(&input, &ReadOptions::new()).expect("PDB fixture");
    assert!(write_mmtf(&structure).is_err());
}

#[test]
fn reader_refuses_unknown_entity_types() {
    let mut file: File = rmp_serde::from_slice(&minimal_mmtf()).expect("fixture MessagePack");
    file.entity_list.as_mut().expect("fixture entity list")[0].kind = "compatibility-water".into();
    let bytes = rmp_serde::to_vec_named(&file).expect("fixture MessagePack");
    let input = InputBuffer::from_bytes(bytes);

    let diagnostics = crate::read_mmtf(&input, &ReadOptions::new())
        .expect_err("unknown entity types must not silently become unknown chemistry");

    assert!(diagnostics.iter().any(|diagnostic| {
        diagnostic
            .message()
            .contains("unsupported MMTF entity type")
    }));
}

fn minimal_mmtf() -> Vec<u8> {
    let file = File {
        mmtf_version: "1.0.0".into(),
        mmtf_producer: "fixture".into(),
        structure_id: Some("TEST".into()),
        title: Some("faithful roundtrip".into()),
        unit_cell: None,
        space_group: None,
        experimental_methods: None,
        resolution: None,
        num_bonds: 1,
        num_atoms: 2,
        num_groups: 1,
        num_chains: 1,
        num_models: 1,
        group_list: vec![Group {
            formal_charge_list: vec![0, 0],
            atom_name_list: vec!["N".into(), "CA".into()],
            element_list: Some(vec!["N".into(), "C".into()]),
            bond_atom_list: Vec::new(),
            bond_order_list: Vec::new(),
            bond_resonance_list: Vec::new(),
            name: "GLY".into(),
            single_letter_code: "G".into(),
            chem_comp_type: "L-peptide linking".into(),
        }],
        bond_atom_list: Some(encode_i32(&[0, 1]).expect("encoding")),
        bond_order_list: Some(encode_i32(&[1]).expect("encoding")),
        bond_resonance_list: Some(encode_i32(&[0]).expect("encoding")),
        x_coord_list: encode_f32(&[11.104, 12.560]).expect("encoding"),
        y_coord_list: encode_f32(&[13.207, 13.100]).expect("encoding"),
        z_coord_list: encode_f32(&[9.301, 9.250]).expect("encoding"),
        b_factor_list: Some(encode_f32(&[20.0, 21.0]).expect("encoding")),
        atom_id_list: Some(encode_i32(&[1, 2]).expect("encoding")),
        alt_loc_list: None,
        occupancy_list: Some(encode_f32(&[1.0, 0.5]).expect("encoding")),
        group_id_list: encode_i32(&[1]).expect("encoding"),
        group_type_list: encode_i32(&[0]).expect("encoding"),
        ins_code_list: None,
        sequence_index_list: Some(encode_i32(&[0]).expect("encoding")),
        chain_id_list: encode_strings(&["A".into()], 4).expect("chain encoding"),
        chain_name_list: None,
        groups_per_chain: vec![1],
        chains_per_model: vec![1],
        entity_list: Some(vec![Entity {
            chain_index_list: vec![0],
            description: "glycine".into(),
            kind: "polymer".into(),
            sequence: "G".into(),
        }]),
    };
    rmp_serde::to_vec_named(&file).expect("fixture MessagePack")
}

#[test]
fn official_reference_fixture_decodes_all_standard_columns() {
    let Some(directory) = official_mmtf_corpus() else {
        return;
    };
    let path = directory.join("3NJW.mmtf");
    if !path.exists() {
        return;
    }
    let bytes = std::fs::read(path).expect("official fixture bytes");
    let input = InputBuffer::from_bytes(bytes);
    let (structure, _) = crate::read_mmtf(&input, &ReadOptions::new()).expect("official MMTF");
    assert_eq!(structure.atom_count(), 169);
    assert_eq!(structure.residue_count(), 44);
    assert_eq!(structure.chain_count(), 2);
    assert_eq!(structure.data().entry.id.as_deref(), Some("3NJW"));
}

#[test]
fn official_reference_corpus_including_ragged_models_decodes() {
    let Some(directory) = official_mmtf_corpus() else {
        return;
    };
    if !directory.exists() {
        return;
    }
    let entries = std::fs::read_dir(directory).expect("fixture directory");
    for entry in entries {
        let path = entry.expect("fixture entry").path();
        let Some(name) = path.file_name().and_then(std::ffi::OsStr::to_str) else {
            continue;
        };
        if !path
            .extension()
            .is_some_and(|extension| extension.eq_ignore_ascii_case("mmtf"))
            || name.starts_with("empty-")
        {
            continue;
        }
        let input = InputBuffer::from_bytes(std::fs::read(&path).expect("fixture bytes"));
        if name.contains("onlyrequired") {
            assert!(crate::read_mmtf(&input, &ReadOptions::new()).is_err());
            continue;
        }
        crate::read_mmtf(&input, &ReadOptions::new())
            .unwrap_or_else(|findings| panic!("{name}: {findings:?}"));
    }
}

fn official_mmtf_corpus() -> Option<std::path::PathBuf> {
    std::env::var_os("MOLFRAME_MMTF_CORPUS").map(std::path::PathBuf::from)
}

fn repeated_glycine_mmtf() -> Vec<u8> {
    let file = File {
        mmtf_version: "1.0.0".into(),
        mmtf_producer: "fixture".into(),
        structure_id: Some("REPEAT".into()),
        title: Some("repeated chemistry".into()),
        unit_cell: None,
        space_group: None,
        experimental_methods: None,
        resolution: None,
        num_bonds: 0,
        num_atoms: 6,
        num_groups: 3,
        num_chains: 1,
        num_models: 1,
        group_list: vec![Group {
            formal_charge_list: vec![0, 0],
            atom_name_list: vec!["N".into(), "CA".into()],
            element_list: Some(vec!["N".into(), "C".into()]),
            bond_atom_list: Vec::new(),
            bond_order_list: Vec::new(),
            bond_resonance_list: Vec::new(),
            name: "GLY".into(),
            single_letter_code: "G".into(),
            chem_comp_type: "L-peptide linking".into(),
        }],
        bond_atom_list: None,
        bond_order_list: None,
        bond_resonance_list: None,
        x_coord_list: encode_f32(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).expect("encoding"),
        y_coord_list: encode_f32(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).expect("encoding"),
        z_coord_list: encode_f32(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).expect("encoding"),
        b_factor_list: None,
        atom_id_list: Some(encode_i32(&[1, 2, 3, 4, 5, 6]).expect("encoding")),
        alt_loc_list: None,
        occupancy_list: None,
        group_id_list: encode_i32(&[1, 2, 3]).expect("encoding"),
        group_type_list: encode_i32(&[0, 0, 0]).expect("encoding"),
        ins_code_list: None,
        sequence_index_list: Some(encode_i32(&[0, 1, 2]).expect("encoding")),
        chain_id_list: encode_strings(&["A".into()], 4).expect("chain encoding"),
        chain_name_list: None,
        groups_per_chain: vec![3],
        chains_per_model: vec![1],
        entity_list: Some(vec![Entity {
            chain_index_list: vec![0],
            description: "triglycine".into(),
            kind: "polymer".into(),
            sequence: "GGG".into(),
        }]),
    };
    rmp_serde::to_vec_named(&file).expect("fixture MessagePack")
}
