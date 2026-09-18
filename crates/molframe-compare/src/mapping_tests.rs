use super::assignment::optimal_pairs;
use super::{chain_sequences, map_chains, map_sequence_to_structure};
use molframe_chem::{Component, ComponentKind, MemoryProvider};
use molframe_core::contract::DictionaryVersion;
use molframe_core::contract::Namespace;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use molframe_seq::Scoring;
use std::sync::Arc;

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn structure(body: &str) -> Structure {
    let input = InputBuffer::from_bytes(format!("{HEADER}{body}").into_bytes());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn provider() -> MemoryProvider {
    MemoryProvider::new(
        DictionaryVersion::new("test-ccd"),
        [
            component("ALA", ComponentKind::AminoAcid, Some(b'A')),
            component("GLY", ComponentKind::AminoAcid, Some(b'G')),
            component("SER", ComponentKind::AminoAcid, Some(b'S')),
            component("TRP", ComponentKind::AminoAcid, Some(b'W')),
            component("HOH", ComponentKind::Solvent, None),
        ],
    )
    .expect("component fixtures are unique")
}

fn component(id: &str, kind: ComponentKind, code: Option<u8>) -> Component {
    Component {
        id: id.into(),
        name: id.into(),
        kind,
        parent: None,
        one_letter_code: code,
        formula: None,
        atoms: Arc::from([]),
        bonds: Arc::from([]),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

// Chain A: Ala-Gly-Ser, plus a waters that must be ignored.
const CHAIN_A: &str = "\
ATOM 1 C CA ALA A 1 0 0 0\n\
ATOM 2 C CA GLY A 2 1 0 0\n\
ATOM 3 C CA SER A 3 2 0 0\n\
HETATM 4 O O HOH A 4 9 9 9\n";

#[test]
fn a_chain_reduces_to_its_one_letter_sequence() {
    let chains =
        chain_sequences(&structure(CHAIN_A), &provider(), Namespace::Label).expect("CCD sequence");
    assert_eq!(chains.len(), 1);
    assert_eq!(chains[0].label, "A");
    assert_eq!(chains[0].sequence, b"AGS");
}

#[test]
fn identical_chains_map_with_full_identity() {
    let reference = structure(CHAIN_A);
    // Same sequence under a different chain label.
    let target = structure(
        "ATOM 1 C CA ALA B 1 0 0 0\n\
ATOM 2 C CA GLY B 2 1 0 0\n\
ATOM 3 C CA SER B 3 2 0 0\n",
    );
    let mappings = map_chains(
        &reference,
        &target,
        &provider(),
        Namespace::Label,
        Scoring::simple(),
        0.5,
    )
    .expect("CCD mapping");
    assert_eq!(mappings.len(), 1);
    assert_eq!(mappings[0].reference, "A");
    assert_eq!(mappings[0].target, "B");
    assert!((mappings[0].identity - 1.0).abs() < 1e-12);
}

#[test]
fn a_query_maps_onto_the_matching_residues() {
    // The structure chain is Ala-Gly-Ser (residues 0,1,2); the query "AGS"
    // aligns onto them one to one.
    let structure = structure(
        "ATOM 1 C CA ALA A 1 0 0 0\n\
ATOM 2 C CA GLY A 2 1 0 0\n\
ATOM 3 C CA SER A 3 2 0 0\n",
    );
    let matches = map_sequence_to_structure(
        b"AGS",
        &structure,
        &provider(),
        molframe_core::contract::Namespace::Label,
        Scoring::simple(),
    )
    .expect("CCD mapping");
    let pairs: Vec<(usize, u32)> = matches
        .iter()
        .map(|m| (m.query_position, m.residue.get()))
        .collect();
    assert_eq!(pairs, vec![(0, 0), (1, 1), (2, 2)]);
}

#[test]
fn a_fragment_maps_onto_its_part_of_the_structure() {
    // Only the "GS" fragment: it should land on residues 1 and 2, not 0.
    let structure = structure(
        "ATOM 1 C CA ALA A 1 0 0 0\n\
ATOM 2 C CA GLY A 2 1 0 0\n\
ATOM 3 C CA SER A 3 2 0 0\n",
    );
    let matches = map_sequence_to_structure(
        b"GS",
        &structure,
        &provider(),
        molframe_core::contract::Namespace::Label,
        Scoring::simple(),
    )
    .expect("CCD mapping");
    let residues: Vec<u32> = matches.iter().map(|m| m.residue.get()).collect();
    assert_eq!(residues, vec![1, 2]);
}

#[test]
fn sequence_mapping_rejects_an_unqualified_explicit_namespace() {
    let result = map_sequence_to_structure(
        b"A",
        &structure("ATOM 1 C CA ALA A 1 0 0 0\n"),
        &provider(),
        Namespace::Explicit,
        Scoring::simple(),
    );
    assert!(result.is_err());
}

#[test]
fn a_dissimilar_chain_below_the_threshold_is_not_mapped() {
    let reference = structure(CHAIN_A);
    let target = structure(
        "ATOM 1 C CA TRP B 1 0 0 0\n\
ATOM 2 C CA TRP B 2 1 0 0\n\
ATOM 3 C CA TRP B 3 2 0 0\n",
    );
    let mappings = map_chains(
        &reference,
        &target,
        &provider(),
        Namespace::Label,
        Scoring::simple(),
        0.9,
    )
    .expect("CCD mapping");
    assert!(mappings.is_empty());
}

#[test]
fn homomeric_targets_each_take_their_own_reference() {
    // Two identical reference chains and two identical target chains: each of the
    // two mappings is one-to-one, not both targets onto one reference.
    let reference = structure(
        "ATOM 1 C CA ALA A 1 0 0 0\n\
ATOM 2 C CA GLY A 2 1 0 0\n\
ATOM 3 C CA ALA C 1 0 0 0\n\
ATOM 4 C CA GLY C 2 1 0 0\n",
    );
    let target = structure(
        "ATOM 1 C CA ALA X 1 0 0 0\n\
ATOM 2 C CA GLY X 2 1 0 0\n\
ATOM 3 C CA ALA Y 1 0 0 0\n\
ATOM 4 C CA GLY Y 2 1 0 0\n",
    );
    let mappings = map_chains(
        &reference,
        &target,
        &provider(),
        Namespace::Label,
        Scoring::simple(),
        0.5,
    )
    .expect("CCD mapping");
    assert_eq!(mappings.len(), 2);
    let targets: Vec<&str> = mappings.iter().map(|m| m.target.as_str()).collect();
    assert!(
        targets.contains(&"X") && targets.contains(&"Y"),
        "each target used once"
    );
}

#[test]
fn assignment_maximises_total_identity_instead_of_taking_a_greedy_pair() {
    let pairs = optimal_pairs(&[vec![0.9, 0.8], vec![0.85, 0.1]], 0.0);
    assert_eq!(pairs, vec![(0, 1), (1, 0)]);
}
