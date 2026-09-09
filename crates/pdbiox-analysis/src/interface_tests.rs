use super::chain_interface;
use pdbiox_core::ExecutionContext;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::structure::Structure;
use pdbiox_spatial::SpatialBackend;
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
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
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
    let structure = match pdbiox_cif::read(&input, &ReadOptions::new()) {
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
