use super::*;
use pdbiox_core::diagnostic::{Code, Diagnostic};
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
