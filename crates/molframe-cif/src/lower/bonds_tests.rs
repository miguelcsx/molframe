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

/// One covalent bond, one disulfide and one metal coordination in a single table.
const MIXED_KINDS: &str = "data_mixed_kinds\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA GLY A 1 0 0 0\n\
ATOM 2 N N GLY A 1 1 0 0\n\
ATOM 3 S SG CYS A 2 2 0 0\n\
ATOM 4 S SG CYS A 3 3 0 0\n\
HETATM 5 FE FE HEM A 4 4 0 0\n\
loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n\
_struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_atom_id\n\
_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n\
_struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_atom_id\n\
_struct_conn.pdbx_value_order\n\
1 covale A 1 GLY CA A 1 GLY N SING\n\
2 disulf A 2 CYS SG A 3 CYS SG SING\n\
3 metalc A 4 HEM FE A 1 GLY CA SING\n";

/// One partner names operator 1 with a translation along b.
const IMAGE_SOURCE: &str = "data_image\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA GLY A 1 0 0 0\nATOM 2 N N GLY A 1 1 0 0\n\
loop_\n_space_group_symop.id\n_space_group_symop.operation_xyz\n1 'x,y,z'\n\
loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n\
_struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_atom_id\n\
_struct_conn.ptnr1_symmetry\n\
_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n\
_struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_atom_id\n\
_struct_conn.ptnr2_symmetry\n\
1 covale A 1 GLY CA 1_555 A 1 GLY N 1_565\n";

/// The same connection with both partners in this asymmetric unit.
const IDENTITY_SYMMETRY_SOURCE: &str = "data_identity\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA GLY A 1 0 0 0\nATOM 2 N N GLY A 1 1 0 0\n\
loop_\n_space_group_symop.id\n_space_group_symop.operation_xyz\n1 'x,y,z'\n\
loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n\
_struct_conn.ptnr1_label_comp_id\n_struct_conn.ptnr1_label_atom_id\n\
_struct_conn.ptnr1_symmetry\n\
_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n\
_struct_conn.ptnr2_label_comp_id\n_struct_conn.ptnr2_label_atom_id\n\
_struct_conn.ptnr2_symmetry\n\
1 covale A 1 GLY CA 1_555 A 1 GLY N 1_555\n";

/// The two-atom topology of [`SOURCE`], without its connectivity.
fn source_data() -> molframe_core::StructureData {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let (structure, _) = crate::read(&input, &ReadOptions::new()).expect("fixture reads");
    let mut data = structure.data().clone();
    data.bonds = molframe_core::BondTableBuilder::new().finish();
    data
}

#[test]
fn struct_conn_resolves_label_endpoints_order_and_file_provenance() {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let (structure, _) = crate::read(&input, &ReadOptions::new()).expect("fixture reads");
    let bond = structure.data().bonds.iter().next().expect("bond exists");
    assert_eq!((bond.atom_a.get(), bond.atom_b.get()), (0, 1));
    assert_eq!(bond.order, molframe_core::BondOrder::Single);
    assert_eq!(bond.provenance, molframe_core::BondProvenance::File);
    assert_eq!(
        bond.kind.get().and_then(|kind| structure.resolve(kind)),
        Some("covale"),
    );
}

#[test]
fn a_deposited_metal_coordination_keeps_its_own_class() {
    // Four 4HHB-style metal sites and one covalent bond in one table: no single
    // global connection type can describe them, and writing `covale` for the
    // metal sites would reclassify coordination as an ordinary covalent edge.
    let input = InputBuffer::from_bytes(MIXED_KINDS.as_bytes().to_vec());
    let (structure, findings) = crate::read(&input, &ReadOptions::new()).expect("fixture reads");
    assert!(findings.is_empty(), "{findings:?}");
    let kinds: Vec<_> = structure
        .data()
        .bonds
        .iter()
        .map(|bond| {
            bond.kind
                .get()
                .and_then(|kind| structure.resolve(kind))
                .unwrap_or("<unclassified>")
        })
        .collect();
    // Sorted by endpoint pair: (0,1) covalent, (0,4) metal, (2,3) disulfide.
    assert_eq!(kinds, ["covale", "metalc", "disulf"]);

    let options = crate::CifWriteOptions::new()
        .with_block_id("mixed")
        .with_generated_connection_ids();
    let written = crate::write_canonical_with_options(&structure, &options)
        .unwrap_or_else(|error| panic!("write failed: {error}"));
    for class in ["covale", "disulf", "metalc"] {
        assert!(
            written.contains(class),
            "`{class}` was flattened away:\n{written}"
        );
    }
    let input = InputBuffer::from_bytes(written.into_bytes());
    let (round_trip, _) = crate::read(&input, &ReadOptions::new()).expect("output reads");
    let after: Vec<_> = round_trip
        .data()
        .bonds
        .iter()
        .map(|bond| {
            bond.kind
                .get()
                .and_then(|kind| round_trip.resolve(kind))
                .unwrap_or("<unclassified>")
        })
        .collect();
    assert_eq!(after, kinds);
}

#[test]
fn a_bond_without_a_deposited_class_is_written_as_unknown() {
    // A caller-created edge with no classification must stay unclassified rather
    // than acquire `covale` on the way out.
    let mut data = source_data();
    let mut bonds = molframe_core::BondTableBuilder::new();
    bonds.push(molframe_core::BondRecord::unclassified(
        molframe_core::AtomIndex::new(0),
        molframe_core::AtomIndex::new(1),
        molframe_core::BondOrder::Single,
        molframe_core::BondProvenance::User,
    ));
    data.bonds = bonds.finish();
    let structure = molframe_core::Structure::new(data);
    let options = crate::CifWriteOptions::new()
        .with_block_id("plain")
        .with_generated_connection_ids();
    let written = crate::write_canonical_with_options(&structure, &options)
        .unwrap_or_else(|error| panic!("write failed: {error}"));
    let row = written
        .lines()
        .find(|line| line.starts_with("1 "))
        .unwrap_or_else(|| panic!("one connection row:\n{written}"));
    assert!(row.contains(" ? "), "{row}");
    assert!(!row.contains("covale"), "{row}");
}

#[test]
fn a_connection_to_another_crystal_image_is_skipped_with_a_warning_not_flattened() {
    // `1_565` is a translation along b, so the partner is an image, not the atom
    // that shares its identifiers. Lowering it here would invent a local bond.
    // The read still succeeds in every mode: all atoms are intact.
    let input = InputBuffer::from_bytes(IMAGE_SOURCE.as_bytes().to_vec());
    for mode in [
        molframe_core::io::ParseMode::Strict,
        molframe_core::io::ParseMode::Permissive,
        molframe_core::io::ParseMode::Recover,
    ] {
        let (structure, findings) = crate::read(&input, &ReadOptions::new().mode(mode))
            .unwrap_or_else(|failure| panic!("{mode:?} must keep the atoms: {failure:?}"));
        assert!(structure.atom_count() > 0);
        assert_eq!(structure.data().bonds.len(), 0, "no local bond is invented");
        let finding = findings
            .iter()
            .find(|finding| finding.code() == molframe_core::Code::W3013)
            .unwrap_or_else(|| panic!("{findings:?}"));
        assert_eq!(finding.severity(), molframe_core::Severity::Loose);
    }
}

#[test]
fn an_identity_symmetry_code_still_lowers_to_a_local_bond() {
    // The archive writes `1_555` for a partner already in this asymmetric unit.
    let input = InputBuffer::from_bytes(IDENTITY_SYMMETRY_SOURCE.as_bytes().to_vec());
    let (structure, findings) = crate::read(&input, &ReadOptions::new()).expect("fixture reads");
    assert!(findings.is_empty(), "{findings:?}");
    assert_eq!(structure.data().bonds.len(), 1);
}

#[test]
fn canonical_round_trip_keeps_connectivity_order_and_provenance() {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let (structure, _) = crate::read(&input, &ReadOptions::new()).expect("fixture reads");
    let options = crate::CifWriteOptions::new()
        .with_block_id("bond")
        .with_generated_connection_ids();
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
