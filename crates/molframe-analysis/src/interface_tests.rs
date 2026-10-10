use molframe_core::ExecutionContext;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use molframe_spatial::SpatialBackend;
use std::fmt::Write;

// Chain A residue 1 (index 0) sits 3 Å from chain B residue 1 (index 1); chain B
// residue 2 (index 2) is far away.
const SOURCE: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA GLY A 1 0 0 0\n\
ATOM 2 C CA GLY B 1 3 0 0\n\
ATOM 3 C CA GLY B 2 100 0 0\n";

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

#[test]
fn the_interface_names_the_touching_residues_of_both_chains() {
    let Ok(residues) = chain_interface(
        &structure(),
        "A",
        "B",
        4.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    let ordinals: Vec<u32> = residues.iter().map(|residue| residue.get()).collect();
    assert_eq!(ordinals, vec![0, 1]);
}

#[test]
fn chains_that_never_approach_have_no_interface() {
    let Ok(residues) = chain_interface(
        &structure(),
        "A",
        "B",
        1.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert!(residues.is_empty());
}

#[test]
fn an_absent_chain_yields_no_interface() {
    let Ok(residues) = chain_interface(
        &structure(),
        "A",
        "Z",
        4.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert!(residues.is_empty());
}

#[test]
fn the_interface_residue_list_is_identical_at_every_worker_count() {
    // Two interleaved chains across a lattice, so the cross-selection search
    // spans several query blocks.
    let mut source = String::from(
        "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
    );
    let mut serial_id = 0u32;
    for x in 0..7i16 {
        for y in 0..7i16 {
            for z in 0..7i16 {
                serial_id += 1;
                let chain = if serial_id.is_multiple_of(2) {
                    "A"
                } else {
                    "B"
                };
                writeln!(
                    source,
                    "ATOM {serial_id} C CA GLY {chain} {serial_id} {:.3} {:.3} {:.3}",
                    f32::from(x) * 2.5,
                    f32::from(y) * 2.5,
                    f32::from(z) * 2.5
                )
                .expect("fixture row");
            }
        }
    }

    let input = InputBuffer::from_bytes(source.into_bytes());
    let structure = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let Ok(serial) = chain_interface(
        &structure,
        "A",
        "B",
        4.0,
        SpatialBackend::CellList,
        &ExecutionContext::default(),
    ) else {
        panic!("valid lattice");
    };
    assert!(!serial.is_empty(), "the fixture must produce an interface");

    for workers in [2, 4, 8, 16] {
        let context = match ExecutionContext::builder().worker_budget(workers).build() {
            Ok(context) => context,
            Err(error) => panic!("valid execution context: {error}"),
        };
        let Ok(parallel) = chain_interface(
            &structure,
            "A",
            "B",
            4.0,
            SpatialBackend::CellList,
            &context,
        ) else {
            panic!("valid lattice");
        };
        assert_eq!(
            parallel, serial,
            "worker count {workers} changed the result"
        );
    }
}

use super::{
    ChainSelector, InterfaceError, InterfaceOptions, chain_interface, chain_interface_with_options,
};

fn read(source: &str) -> Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

// Label A is author B and label B is author A, so the bare name "A" means two
// different chains depending on the namespace.
const SWAPPED: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.auth_asym_id\n_atom_site.label_seq_id\n\
_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA GLY A B 1 0 0 0\n\
ATOM 2 C CA GLY B A 1 3 0 0\n\
ATOM 3 C CA GLY B A 2 100 0 0\n";

fn options<'a>(first: ChainSelector<'a>, second: ChainSelector<'a>) -> InterfaceOptions<'a> {
    InterfaceOptions {
        first,
        second,
        cutoff: 4.0,
        backend: SpatialBackend::BruteForce,
        periodic: false,
    }
}

fn run(structure: &Structure, options: &InterfaceOptions<'_>) -> Result<Vec<u32>, InterfaceError> {
    chain_interface_with_options(structure, options, &ExecutionContext::default())
        .map(|residues| residues.iter().map(|residue| residue.get()).collect())
}

#[test]
fn a_name_that_means_different_chains_per_namespace_is_rejected() {
    let structure = read(SWAPPED);
    let result = chain_interface(
        &structure,
        "A",
        "B",
        4.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    );
    assert!(matches!(result, Err(InterfaceError::AmbiguousChain(_))));
}

#[test]
fn an_explicit_namespace_resolves_swapped_identifiers() {
    let structure = read(SWAPPED);
    let by_label = run(
        &structure,
        &options(ChainSelector::label("A"), ChainSelector::label("B")),
    );
    assert!(matches!(by_label, Ok(ref found) if found == &[0, 1]));
    let by_auth = run(
        &structure,
        &options(ChainSelector::auth("A"), ChainSelector::auth("B")),
    );
    assert!(matches!(by_auth, Ok(ref found) if found == &[0, 1]));
    // Label A is the single-residue chain; only its own partner is returned.
    let mixed = run(
        &structure,
        &options(ChainSelector::label("A"), ChainSelector::auth("A")),
    );
    assert!(matches!(mixed, Ok(ref found) if found == &[0, 1]));
}

#[test]
fn an_interface_of_a_chain_with_itself_is_rejected() {
    let structure = read(SWAPPED);
    let same_name = run(
        &structure,
        &options(ChainSelector::label("A"), ChainSelector::label("A")),
    );
    assert!(matches!(same_name, Err(InterfaceError::SameChain)));
    // Label A and author B are the same chain under two names.
    let same_chain = run(
        &structure,
        &options(ChainSelector::label("A"), ChainSelector::auth("B")),
    );
    assert!(matches!(same_chain, Err(InterfaceError::SameChain)));
}

fn cell_source(lengths: f32) -> String {
    format!(
        "data_s\n_cell.length_a {lengths}\n_cell.length_b {lengths}\n_cell.length_c {lengths}\n\
_cell.angle_alpha 90\n_cell.angle_beta 90\n_cell.angle_gamma 90\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA GLY A 1 0.5 0 0\n\
ATOM 2 C CA GLY B 1 9.5 0 0\n"
    )
}

#[test]
fn periodic_interface_finds_contacts_across_the_cell_boundary() {
    let structure = read(&cell_source(10.0));
    let mut request = options(ChainSelector::label("A"), ChainSelector::label("B"));
    request.cutoff = 2.0;
    assert!(matches!(run(&structure, &request), Ok(ref found) if found.is_empty()));
    request.periodic = true;
    assert!(matches!(run(&structure, &request), Ok(ref found) if found == &[0, 1]));
}

#[test]
fn periodic_interface_rejects_a_placeholder_cell() {
    let structure = read(&cell_source(1.0));
    let mut request = options(ChainSelector::label("A"), ChainSelector::label("B"));
    request.periodic = true;
    assert!(matches!(
        run(&structure, &request),
        Err(InterfaceError::PlaceholderCell)
    ));
}

#[test]
fn periodic_interface_requires_a_cell() {
    let structure = read(SOURCE);
    let mut request = options(ChainSelector::label("A"), ChainSelector::label("B"));
    request.periodic = true;
    assert!(matches!(
        run(&structure, &request),
        Err(InterfaceError::MissingCell)
    ));
}
