use super::*;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::index::EntityIndex;
use molframe_core::io::{AmbiguousResidueBoundaryPolicy, MissingElementPolicy, ParseMode};
use molframe_core::structure::{
    AtomRef, ResidueRef, Structure, StructureDifferenceOptions, structure_difference,
};

#[path = "structure_model_tests.rs"]
mod model_tests;

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
fn direct_and_lossless_single_model_reads_are_semantically_identical() {
    let input = InputBuffer::from_bytes(DIPEPTIDE.as_bytes().to_vec());
    let options = ReadOptions::new();
    let (direct, direct_findings) = read(&input, &options).expect("direct read should succeed");
    let (_, lossless, lossless_findings) =
        read_with_document(&input, &options).expect("lossless read should succeed");
    let difference = structure_difference(
        &direct,
        &lossless,
        StructureDifferenceOptions {
            coordinate_tolerance: 0.0,
        },
    )
    .expect("zero is a valid coordinate tolerance");

    assert!(difference.is_empty(), "difference: {difference:?}");
    assert_eq!(direct_findings, lossless_findings);
}

#[test]
fn metadata_after_atom_site_is_reconciled_without_reparsing() {
    let text = "\
data_late
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
#
_entry.id LATE
loop_
_entity.id
_entity.type
7 polymer
#
loop_
_struct_asym.id
_struct_asym.entity_id
A 7
#
_cell.length_a 10
_cell.length_b 11
_cell.length_c 12
_cell.angle_alpha 90
_cell.angle_beta 90
_cell.angle_gamma 90
";
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let options = ReadOptions::new();
    let (direct, direct_findings) = read(&input, &options).expect("direct read should succeed");
    let (_, lossless, lossless_findings) =
        read_with_document(&input, &options).expect("lossless read should succeed");
    let difference = structure_difference(
        &direct,
        &lossless,
        StructureDifferenceOptions {
            coordinate_tolerance: 0.0,
        },
    )
    .expect("zero is a valid coordinate tolerance");

    assert!(difference.is_empty(), "difference: {difference:?}");
    assert_eq!(direct_findings, lossless_findings);
    assert_eq!(direct.data().entry.id.as_deref(), Some("LATE"));
    assert_eq!(
        direct.data().topology.entities.kind(EntityIndex::new(0)),
        Some(molframe_core::topology::EntityKind::Polymer)
    );
    assert_eq!(
        direct.data().cell.map(|cell| cell.lengths),
        Some([10.0, 11.0, 12.0])
    );
}

#[test]
fn atom_site_can_follow_unrelated_metadata_loops() {
    let text = "data_x\n\
loop_\n_entity.id\n_entity.type\n1 polymer\n#\n\
loop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n\
_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
1 C CA 0 0 0\n#\n";
    let (structure, _) = parse_structure(text);
    assert_eq!(structure.atom_count(), 1);
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
fn author_identifiers_are_not_copied_into_missing_label_identifiers() {
    let text = "data_x\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.auth_atom_id\n_atom_site.label_comp_id\n_atom_site.auth_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C AUTHOR_ATOM LIG AUTHOR_CHAIN 1 0 0 0\n";
    let (structure, _) = parse_structure(text);
    let Some(atom) = structure.atom(molframe_core::index::AtomIndex::new(0)) else {
        panic!("fixture atom missing")
    };
    let Some(chain) = structure.chain(molframe_core::index::ChainIndex::new(0)) else {
        panic!("fixture chain missing")
    };
    assert_ne!(atom.name(), Some("AUTHOR_ATOM"));
    assert_eq!(atom.auth_name(), Some("AUTHOR_ATOM"));
    assert_ne!(chain.label(), Some("AUTHOR_CHAIN"));
    assert_eq!(chain.auth_label(), Some("AUTHOR_CHAIN"));
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
        .and_then(molframe_core::structure::ChainRef::entity);
    assert_eq!(entity, Some(EntityIndex::new(0)));
}

#[test]
fn ambiguous_residue_boundary_requires_explicit_file_order_inference() {
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
    let implicit = read(&input, &ReadOptions::new().mode(ParseMode::Permissive));
    assert!(implicit.is_err());

    let explicit = ReadOptions::new()
        .ambiguous_residue_boundary_policy(AmbiguousResidueBoundaryPolicy::InferFromFileOrder);
    let (_, findings) = read(&input, &explicit).expect("explicit inference should read");
    assert!(findings.iter().any(|finding| finding.code() == Code::W3011));
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

    let options = crate::CifWriteOptions::new().with_block_id("reread");
    let written = crate::write_canonical_with_options(&structure, &options)
        .unwrap_or_else(|error| panic!("write failed: {error}"));
    let (reread, _) = parse_structure(&written);
    let reread_components: Vec<_> = reread
        .data()
        .atoms()
        .filter_map(AtomRef::component_name)
        .collect();
    assert_eq!(reread_components, ["GLY", "ALA"]);
}

#[test]
fn non_row_atom_site_layouts_are_rejected_without_a_document_fallback() {
    let scalar = "data_x\n\
_atom_site.id 1\n_atom_site.type_symbol C\n_atom_site.label_atom_id CA\n\
_atom_site.label_comp_id GLY\n_atom_site.label_asym_id A\n_atom_site.label_seq_id 1\n\
_atom_site.Cartn_x 0\n_atom_site.Cartn_y 0\n_atom_site.Cartn_z 0\n";
    let repeated_loop = "data_x\n\
loop_\n_atom_site.id\n_atom_site.Cartn_x\n1 0\n\
loop_\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n0 0\n";
    let duplicate = "data_x\n\
loop_\n_atom_site.id\n_atom_site.id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n\
_atom_site.Cartn_z\n1 1 0 0 0\n";

    for text in [scalar, repeated_loop, duplicate] {
        let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
        let findings = read(&input, &ReadOptions::new())
            .expect_err("ambiguous atom_site layout must be rejected");
        assert!(findings.iter().any(|finding| finding.code() == Code::E1104));
    }
}

#[test]
fn the_structure_a_read_produces_satisfies_its_own_invariants() {
    let (structure, _) = parse_structure(DIPEPTIDE);
    assert!(molframe_core::structure::validate(structure.data()).is_empty());
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
    let input = InputBuffer::from_bytes(text.as_bytes().to_vec());
    let options =
        ReadOptions::new().missing_element_policy(MissingElementPolicy::InferFromAtomName);
    let (structure, findings) = match read(&input, &options) {
        Ok(result) => result,
        Err(findings) => panic!("read failed: {findings:?}"),
    };
    assert_eq!(structure.atom_count(), 1);
    assert_eq!(
        structure.data().atoms().next().and_then(AtomRef::element),
        Some(molframe_core::Element::CALCIUM),
    );
    assert!(findings.iter().any(|finding| finding.code() == Code::W3203));
}
