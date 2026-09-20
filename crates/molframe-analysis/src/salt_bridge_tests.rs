use super::salt_bridges;
use molframe_core::ExecutionContext;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use molframe_spatial::SpatialBackend;
use std::fmt::Write;

fn structure(source: &str) -> Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

#[test]
fn an_aspartate_oxygen_near_a_lysine_nitrogen_is_a_bridge() {
    let source = format!(
        "{HEADER}\
ATOM 1 O OD1 ASP A 1 0 0 0\n\
ATOM 2 N NZ LYS A 2 3.5 0 0\n"
    );
    let structure = crate::chemistry_test_support::charges(&structure(&source), &[(0, -1), (1, 1)]);
    let Ok(bridges) = salt_bridges(
        &structure,
        4.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert_eq!(bridges.len(), 1);
    let bridge = bridges.row(0).expect("one bridge");
    assert_eq!(bridge.anion.get(), 0);
    assert_eq!(bridge.cation.get(), 1);
    assert!((bridge.distance - 3.5).abs() < 1e-5);
}

#[test]
fn a_pair_beyond_the_cutoff_is_not_a_bridge() {
    let source = format!(
        "{HEADER}\
ATOM 1 O OD1 ASP A 1 0 0 0\n\
ATOM 2 N NZ LYS A 2 5 0 0\n"
    );
    let structure = crate::chemistry_test_support::charges(&structure(&source), &[(0, -1), (1, 1)]);
    let Ok(bridges) = salt_bridges(
        &structure,
        4.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert!(bridges.is_empty());
}

#[test]
fn two_anions_do_not_bridge_each_other() {
    let source = format!(
        "{HEADER}\
ATOM 1 O OD1 ASP A 1 0 0 0\n\
ATOM 2 O OE1 GLU A 2 3 0 0\n"
    );
    let structure =
        crate::chemistry_test_support::charges(&structure(&source), &[(0, -1), (1, -1)]);
    let Ok(bridges) = salt_bridges(
        &structure,
        4.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert!(bridges.is_empty(), "like charges do not bridge");
}

#[test]
fn a_non_charged_side_chain_atom_is_ignored() {
    // Aspartate's CB carbon is not part of the carboxylate.
    let source = format!(
        "{HEADER}\
ATOM 1 C CB ASP A 1 0 0 0\n\
ATOM 2 N NZ LYS A 2 3 0 0\n"
    );
    let structure = crate::chemistry_test_support::charges(&structure(&source), &[(1, 1)]);
    let Ok(bridges) = salt_bridges(
        &structure,
        4.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert!(bridges.is_empty());
}

#[test]
fn the_bridge_list_is_identical_at_every_worker_count() {
    // Alternating charged atoms across a lattice, dense enough that the search
    // spans several query blocks.
    let mut source = String::from(HEADER);
    let mut charges = Vec::new();
    let mut serial_id = 0u32;
    for x in 0..7i16 {
        for y in 0..7i16 {
            for z in 0..7i16 {
                let anion = serial_id.is_multiple_of(2);
                let (name, element) = if anion { ("OD1", "O") } else { ("NZ", "N") };
                serial_id += 1;
                writeln!(
                    source,
                    "ATOM {serial_id} {element} {name} LIG A {serial_id} {:.3} {:.3} {:.3}",
                    f32::from(x) * 3.0,
                    f32::from(y) * 3.0,
                    f32::from(z) * 3.0
                )
                .expect("fixture row");
                charges.push((serial_id - 1, if anion { -1 } else { 1 }));
            }
        }
    }

    let structure = crate::chemistry_test_support::charges(&structure(&source), &charges);
    let Ok(serial) = salt_bridges(
        &structure,
        4.0,
        SpatialBackend::CellList,
        &ExecutionContext::default(),
    ) else {
        panic!("valid lattice");
    };
    assert!(!serial.is_empty(), "the fixture must produce bridges");

    for workers in [2, 4, 8, 16] {
        let context = match ExecutionContext::builder().worker_budget(workers).build() {
            Ok(context) => context,
            Err(error) => panic!("valid execution context: {error}"),
        };
        let Ok(parallel) = salt_bridges(&structure, 4.0, SpatialBackend::CellList, &context) else {
            panic!("valid lattice");
        };
        assert_eq!(
            parallel, serial,
            "worker count {workers} changed the result"
        );
    }
}
