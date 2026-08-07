use super::*;
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::index::EntityIndex;
use pdbiox_core::io::ParseMode;
use pdbiox_core::structure::{AtomRef, ResidueRef, Structure};

const DIPEPTIDE: &str = "\
data_TEST
_entry.id TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N   . GLY A 1 ? 27.340 24.430 2.614 1.00 10.00 1 A 1
ATOM 2 C CA  . GLY A 1 ? 26.266 25.413 2.842 1.00 11.00 1 A 1
ATOM 3 N N   . ASN A 2 ? 26.335 27.770 3.258 1.00 14.00 2 A 1
#
";

fn parse_structure(text: &str) -> (Structure, Vec<Diagnostic>) {
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    match read(&input, &ReadOptions::new()) {
        Ok(result) => result,
        Err(findings) => panic!("read failed: {findings:?}"),
    }
}

#[test]
fn a_dipeptide_reads_with_the_hierarchy_the_file_describes() {
    let (structure, _) = parse_structure(DIPEPTIDE);
    assert_eq!(structure.atom_count(), 3);
    assert_eq!(structure.residue_count(), 2);
    assert_eq!(structure.chain_count(), 1);
    assert_eq!(structure.data().entry.id.as_deref(), Some("TEST"));
}

#[test]
fn names_elements_and_positions_survive_the_read() {
    let (structure, _) = parse_structure(DIPEPTIDE);
    let names: Vec<_> = structure.data().atoms().filter_map(AtomRef::name).collect();
    assert_eq!(names, ["N", "CA", "N"]);

    let residues: Vec<_> = structure
        .data()
        .residues()
        .filter_map(ResidueRef::name)
        .collect();
    assert_eq!(residues, ["GLY", "ASN"]);

    let first = structure.data().atoms().next().and_then(AtomRef::position);
    assert!(first.is_some_and(|position| (position[0] - 27.340).abs() < 1e-3));
}

#[test]
fn both_namespaces_are_kept_rather_than_one_being_derived_from_the_other() {
    let (structure, _) = parse_structure(DIPEPTIDE);
    let Some(residue) = structure.data().residues().next() else {
        panic!("expected a residue")
    };
    assert_eq!(residue.label_seq_id(), Some(1));
    assert_eq!(residue.auth_seq_id(), Some(1));
}

#[test]
fn distinct_chain_namespaces_and_entity_relations_survive_lowering() {
    let text = "\
data_x
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
#
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 2
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N N GLY A 1 10 X 0 0 0
HETATM 2 C C1 LIG B . 20 Y 1 1 1
#
";
    let (structure, _) = parse_structure(text);
    let chains: Vec<_> = structure.data().chains().collect();
    assert_eq!(chains.len(), 2);
    assert_eq!(
        chains[0]
            .label_asym_id()
            .and_then(|symbol| structure.resolve(symbol)),
        Some("A")
    );
    assert_eq!(
        chains[0]
            .auth_asym_id()
            .and_then(|symbol| structure.resolve(symbol)),
        Some("X")
    );
    assert_eq!(chains[0].entity(), Some(EntityIndex::new(0)));
    assert_eq!(chains[1].entity(), Some(EntityIndex::new(1)));
}

#[test]
fn coordinate_only_reads_keep_entity_identity_because_it_is_topology() {
    let text = "\
data_x
loop_
_entity.id
_entity.type
7 polymer
#
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 N N GLY A 7 1 0 0 0
#
";
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let options = ReadOptions::new().only_atomic_coords(true);
    let (structure, _) = match read(&input, &options) {
        Ok(result) => result,
        Err(findings) => panic!("read failed: {findings:?}"),
    };
    let entity = structure
        .data()
        .chains()
        .next()
        .and_then(pdbiox_core::structure::ChainRef::entity);
    assert_eq!(entity, Some(EntityIndex::new(0)));
}

#[test]
fn strict_mode_refuses_an_ambiguous_residue_boundary_that_permissive_mode_reports() {
    let text = "\
data_x
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 N N GLY A 1 0 0 0
2 C CA GLY A 1 1 0 0
3 N N GLY A 1 2 0 0
#
";
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let strict = read(&input, &ReadOptions::new().mode(ParseMode::Strict));
    assert!(strict.is_err());

    let permissive = read(&input, &ReadOptions::new());
    assert!(permissive.is_ok());
}

#[test]
fn a_component_identity_carried_by_an_altloc_survives_structure_and_write() {
    let text = "\
data_x
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 N N A GLY A 1 0 0 0
2 C CA B ALA A 1 1 0 0
#
";
    let (structure, findings) = parse_structure(text);
    let components: Vec<_> = structure
        .data()
        .atoms()
        .filter_map(AtomRef::component_name)
        .collect();
    assert_eq!(components, ["GLY", "ALA"]);
    assert!(findings.iter().any(|finding| finding.code() == Code::W3012));

    let written = crate::write_canonical(&structure);
    let (reread, _) = parse_structure(&written);
    let reread_components: Vec<_> = reread
        .data()
        .atoms()
        .filter_map(AtomRef::component_name)
        .collect();
    assert_eq!(reread_components, ["GLY", "ALA"]);
}

#[test]
fn dense_models_keep_numbers_and_identity_changes_become_ragged() {
    let text = "\
data_x
loop_
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
1 N N GLY A 1 0 0 0 5
2 C CA GLY A 1 1 0 0 5
1 N N GLY A 1 0 1 0 9
2 C CA GLY A 1 1 1 0 9
#
";
    let (structure, _) = parse_structure(text);
    let numbers: Vec<_> = structure
        .data()
        .models()
        .filter_map(pdbiox_core::structure::ModelRef::number)
        .collect();
    assert_eq!(numbers, [5, 9]);
    assert_eq!(structure.model_count(), 2);

    let mismatched = text.replacen("2 C CA GLY A 1 1 1 0 9", "2 O O GLY A 1 1 1 0 9", 1);
    let input = InputBuffer::from_bytes(mismatched.into_bytes());
    let (ragged, _) = match read(&input, &ReadOptions::new().mode(ParseMode::Recover)) {
        Ok(result) => result,
        Err(findings) => panic!("ragged read failed: {findings:?}"),
    };
    let Some(models) = ragged.ragged_models() else {
        panic!("identity-changing models must be ragged")
    };
    assert_eq!(models.len(), 2);
    let names: Vec<Vec<_>> = models
        .iter()
        .map(|model| model.data().atoms().filter_map(AtomRef::name).collect())
        .collect();
    assert_eq!(names, [["N", "CA"], ["N", "O"]]);
    let numbers: Vec<_> = ragged
        .data()
        .models()
        .filter_map(pdbiox_core::structure::ModelRef::number)
        .collect();
    assert_eq!(numbers, [5, 9]);
}

#[test]
fn reading_only_the_first_cif_model_does_not_append_an_empty_frame() {
    let text = "\
data_x
loop_
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
1 N N GLY A 1 0 0 0 4
1 N N GLY A 1 1 1 1 8
#
";
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let options = ReadOptions::new().only_first_model(true);
    let (structure, _) = match read(&input, &options) {
        Ok(result) => result,
        Err(findings) => panic!("read failed: {findings:?}"),
    };
    assert_eq!(structure.model_count(), 1);
    let number = structure
        .data()
        .models()
        .next()
        .and_then(pdbiox_core::structure::ModelRef::number);
    assert_eq!(number, Some(4));
}

#[test]
fn the_structure_a_read_produces_satisfies_its_own_invariants() {
    let (structure, _) = parse_structure(DIPEPTIDE);
    assert!(pdbiox_core::structure::validate(structure.data()).is_empty());
}

#[test]
fn a_block_without_coordinates_is_refused_rather_than_read_as_empty() {
    let input = InputBuffer::from_bytes(b"data_x\n_entry.id X\n".to_vec());
    let refused = read(&input, &ReadOptions::new());
    let codes: Vec<_> = match refused {
        Ok(_) => Vec::new(),
        Err(findings) => findings.iter().map(Diagnostic::code).collect(),
    };
    assert!(codes.contains(&Code::E2001), "got {codes:?}");
}

#[test]
fn a_missing_element_column_is_inferred_and_the_inference_is_reported() {
    let text = "\
data_x
loop_
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 CA GLY A 1 1.0 1.0 1.0
#
";
    let (structure, findings) = parse_structure(text);
    assert_eq!(structure.atom_count(), 1);
    assert!(findings.iter().any(|finding| finding.code() == Code::W3203));
}
