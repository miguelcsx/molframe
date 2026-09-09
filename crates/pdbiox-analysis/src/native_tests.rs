use super::native_contact_fraction;
use pdbiox_core::ExecutionContext;
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

const COMPACT: &str = "\
ATOM 1 C C1 LIG A 1 0 0 0\n\
ATOM 2 C C2 LIG A 2 1 0 0\n\
ATOM 3 C C3 LIG A 3 2 0 0\n";

#[test]
fn a_structure_keeps_all_of_its_own_contacts() {
    let source = format!("{HEADER}{COMPACT}");
    let reference = structure(&source);
    let target = structure(&source);
    let Ok(result) = native_contact_fraction(
        &reference,
        &target,
        1.5,
        1.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert!(result.native > 0);
    assert_eq!(result.kept, result.native);
    assert!((result.fraction - 1.0).abs() < 1e-12);
}

#[test]
fn a_pulled_apart_target_keeps_none() {
    let reference = structure(&format!("{HEADER}{COMPACT}"));
    let spread = "\
ATOM 1 C C1 LIG A 1 0 0 0\n\
ATOM 2 C C2 LIG A 2 50 0 0\n\
ATOM 3 C C3 LIG A 3 100 0 0\n";
    let target = structure(&format!("{HEADER}{spread}"));
    let Ok(result) = native_contact_fraction(
        &reference,
        &target,
        1.5,
        1.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert!(result.native > 0);
    assert_eq!(result.kept, 0);
    assert!(result.fraction.abs() < 1e-12);
}

#[test]
fn structures_of_different_size_cannot_be_compared() {
    let reference = structure(&format!("{HEADER}{COMPACT}"));
    let smaller = "ATOM 1 C C1 LIG A 1 0 0 0\n";
    let target = structure(&format!("{HEADER}{smaller}"));
    assert!(
        native_contact_fraction(
            &reference,
            &target,
            1.5,
            1.0,
            SpatialBackend::BruteForce,
            &ExecutionContext::default(),
        )
        .is_err()
    );
}

#[test]
fn worker_count_does_not_change_q() {
    let source = format!("{HEADER}{COMPACT}");
    let reference = structure(&source);
    let target = structure(&source);
    let serial = match native_contact_fraction(
        &reference,
        &target,
        1.5,
        1.0,
        SpatialBackend::CellList,
        &ExecutionContext::default(),
    ) {
        Ok(result) => result,
        Err(error) => panic!("serial Q failed: {error}"),
    };
    for workers in [1, 2, 4, 8] {
        let context = match ExecutionContext::builder().worker_budget(workers).build() {
            Ok(context) => context,
            Err(error) => panic!("valid execution context: {error}"),
        };
        let parallel = match native_contact_fraction(
            &reference,
            &target,
            1.5,
            1.0,
            SpatialBackend::CellList,
            &context,
        ) {
            Ok(result) => result,
            Err(error) => panic!("parallel Q failed: {error}"),
        };
        assert_eq!(serial, parallel, "worker count {workers} changed Q");
    }
}
