use super::{NativeError, native_contact_fraction, native_contact_fraction_matched};
use molframe_core::ExecutionContext;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use molframe_spatial::SpatialBackend;

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
        0.0,
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
        0.0,
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
            0.0,
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
        0.0,
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
            0.0,
            SpatialBackend::CellList,
            &context,
        ) {
            Ok(result) => result,
            Err(error) => panic!("parallel Q failed: {error}"),
        };
        assert_eq!(serial, parallel, "worker count {workers} changed Q");
    }
}

fn q(reference: &Structure, target: &Structure, cutoff: f32, tolerance: f32) -> f64 {
    match native_contact_fraction(
        reference,
        target,
        cutoff,
        tolerance,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) {
        Ok(result) => result.fraction,
        Err(error) => panic!("Q failed: {error}"),
    }
}

#[test]
fn a_contact_exactly_at_the_cutoff_survives_comparison_with_itself() {
    // Distances 1.1 and 2.3 are not representable in binary, so the pair sits
    // on the cutoff only up to rounding; the comparison must still keep it.
    let source =
        format!("{HEADER}ATOM 1 C C1 LIG A 1 0.1 0.1 0.1\nATOM 2 C C2 LIG A 2 2.4 0.1 0.1\n");
    let reference = structure(&source);
    let target = structure(&source);
    assert!((q(&reference, &target, 2.3, 0.0) - 1.0).abs() < 1e-12);
}

#[test]
fn the_tolerance_is_an_additive_distance_in_angstrom() {
    // Native pair at 1.0 A (cutoff 1.5). The target stretches it to 1.7 A, so it
    // is kept only once cutoff + tolerance reaches 1.7, i.e. tolerance >= 0.2 A.
    let reference = structure(&format!(
        "{HEADER}ATOM 1 C C1 LIG A 1 0 0 0\nATOM 2 C C2 LIG A 2 1 0 0\n"
    ));
    let target = structure(&format!(
        "{HEADER}ATOM 1 C C1 LIG A 1 0 0 0\nATOM 2 C C2 LIG A 2 1.7 0 0\n"
    ));
    assert!(q(&reference, &target, 1.5, 0.0).abs() < 1e-12);
    assert!(q(&reference, &target, 1.5, 0.19).abs() < 1e-12);
    assert!((q(&reference, &target, 1.5, 0.21) - 1.0).abs() < 1e-12);
}

#[test]
fn matched_atoms_are_compared_through_the_correspondence() {
    // Reference chain 0-1-2 along x at unit spacing: contacts (0,1) and (1,2)
    // at cutoff 1.5. The target lists the atoms in reverse order with an extra
    // far atom first; atom 1 of the target is therefore reference atom 2.
    let reference = structure(&format!("{HEADER}{COMPACT}"));
    let positions = [
        [90.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [5.0, 0.0, 0.0],
    ];
    let map = [Some(1), Some(2), Some(3)];
    let Ok(result) = native_contact_fraction_matched(
        &reference,
        &positions,
        &map,
        1.5,
        0.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    // (0,1): 2->1 gives 1.0, kept. (1,2): 1->5 gives 4.0, broken.
    assert_eq!(
        (result.native, result.kept, result.matched_atoms),
        (2, 1, 3)
    );
    assert!((result.fraction - 0.5).abs() < 1e-12);
}

#[test]
fn unmatched_reference_atoms_take_no_part() {
    let reference = structure(&format!("{HEADER}{COMPACT}"));
    let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let map = [Some(0), Some(1), None];
    let Ok(result) = native_contact_fraction_matched(
        &reference,
        &positions,
        &map,
        1.5,
        0.0,
        SpatialBackend::BruteForce,
        &ExecutionContext::default(),
    ) else {
        panic!("valid");
    };
    assert_eq!(
        (result.native, result.kept, result.matched_atoms),
        (1, 1, 2)
    );
}

#[test]
fn a_correspondence_must_fit_both_structures() {
    let reference = structure(&format!("{HEADER}{COMPACT}"));
    let positions = [[0.0, 0.0, 0.0]];
    let run = |map: &[Option<u32>]| {
        native_contact_fraction_matched(
            &reference,
            &positions,
            map,
            1.5,
            0.0,
            SpatialBackend::BruteForce,
            &ExecutionContext::default(),
        )
    };
    assert!(matches!(
        run(&[Some(0)]),
        Err(NativeError::CorrespondenceLength {
            given: 1,
            reference: 3
        })
    ));
    assert!(matches!(
        run(&[Some(0), Some(7), None]),
        Err(NativeError::TargetAtomAbsent {
            index: 7,
            target: 1
        })
    ));
}
