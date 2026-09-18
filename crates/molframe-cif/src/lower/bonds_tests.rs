use molframe_core::io::{InputBuffer, ReadOptions};

const SOURCE: &str = "data_bond\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA GLY A 1 0 0 0\nATOM 2 N N GLY A 1 1 0 0\n\
loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n\
_struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_atom_id\n\
_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n\
_struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_atom_id\n\
_struct_conn.pdbx_value_order\n1 covale A 1 GLY CA A 1 GLY N SING\n";

const NON_POLYMER_SOURCE: &str = "data_non_polymer_bond\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.auth_atom_id\n_atom_site.auth_comp_id\n\
_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
HETATM 1 ZN ZN ZN A . ZN ZN A 401 0 0 0\n\
HETATM 2 O O HOH U . O HOH W 516 1 0 0\n\
HETATM 3 O O HOH U . O HOH W 603 2 0 0\n\
loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n\
_struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_atom_id\n\
_struct_conn.ptnr1_auth_seq_id\n_struct_conn.ptnr2_label_asym_id\n\
_struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_label_comp_id\n\
_struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_auth_seq_id\n\
1 metalc A . ZN ZN 401 U . HOH O 603\n";

#[test]
fn struct_conn_resolves_label_endpoints_order_and_file_provenance() {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let (structure, _) = crate::read(&input, &ReadOptions::new()).expect("fixture reads");
    let bond = structure.data().bonds.iter().next().expect("bond exists");
    assert_eq!((bond.atom_a.get(), bond.atom_b.get()), (0, 1));
    assert_eq!(bond.order, molframe_core::BondOrder::Single);
    assert_eq!(bond.provenance, molframe_core::BondProvenance::File);
}

#[test]
fn canonical_round_trip_keeps_connectivity_order_and_provenance() {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let (structure, _) = crate::read(&input, &ReadOptions::new()).expect("fixture reads");
    let options = crate::CifWriteOptions::new()
        .with_block_id("bond")
        .with_generated_connection_ids()
        .with_connection_type_id("covale");
    let written = crate::write_canonical_with_options(&structure, &options)
        .unwrap_or_else(|error| panic!("write failed: {error}"));
    let input = InputBuffer::from_bytes(written.into_bytes());
    let (round_trip, _) = crate::read(&input, &ReadOptions::new()).expect("output reads");
    assert_eq!(
        structure.data().bonds.iter().collect::<Vec<_>>(),
        round_trip.data().bonds.iter().collect::<Vec<_>>()
    );
}

#[test]
fn struct_conn_uses_author_sequence_when_label_sequence_is_absent() {
    let input = InputBuffer::from_bytes(NON_POLYMER_SOURCE.as_bytes().to_vec());
    let (structure, diagnostics) = crate::read(&input, &ReadOptions::new()).expect("fixture reads");
    assert!(diagnostics.is_empty(), "{diagnostics:?}");
    let bond = structure.data().bonds.iter().next().expect("bond exists");
    assert_eq!((bond.atom_a.get(), bond.atom_b.get()), (0, 2));
}
