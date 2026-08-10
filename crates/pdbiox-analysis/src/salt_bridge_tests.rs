use super::salt_bridges;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::structure::Structure;
use pdbiox_spatial::SpatialBackend;

fn structure(source: &str) -> Structure {
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
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
    let Ok(bridges) = salt_bridges(&structure, 4.0, SpatialBackend::BruteForce) else {
        panic!("valid");
    };
    assert_eq!(bridges.len(), 1);
    assert_eq!(bridges[0].anion.get(), 0);
    assert_eq!(bridges[0].cation.get(), 1);
    assert!((bridges[0].distance - 3.5).abs() < 1e-5);
}

#[test]
fn a_pair_beyond_the_cutoff_is_not_a_bridge() {
    let source = format!(
        "{HEADER}\
ATOM 1 O OD1 ASP A 1 0 0 0\n\
ATOM 2 N NZ LYS A 2 5 0 0\n"
    );
    let structure = crate::chemistry_test_support::charges(&structure(&source), &[(0, -1), (1, 1)]);
    let Ok(bridges) = salt_bridges(&structure, 4.0, SpatialBackend::BruteForce) else {
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
    let Ok(bridges) = salt_bridges(&structure, 4.0, SpatialBackend::BruteForce) else {
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
    let Ok(bridges) = salt_bridges(&structure, 4.0, SpatialBackend::BruteForce) else {
        panic!("valid");
    };
    assert!(bridges.is_empty());
}
