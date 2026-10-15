use super::super::tests::{SOURCE, annotated_structure, water_lattice};
use super::*;

fn options() -> HydrogenBondCandidateOptions {
    HydrogenBondCandidateOptions {
        maximum_donor_acceptor_distance: 3.5,
        maximum_candidates: 100,
        backend: SpatialBackend::CellList,
    }
}

#[test]
fn motion_candidates_include_a_contact_that_exists_only_inside_the_segment() {
    let structure = annotated_structure(SOURCE);
    let start = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [6.0, 0.0, 0.0]];
    let end = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-6.0, 0.0, 0.0]];
    for backend in [
        SpatialBackend::CellList,
        SpatialBackend::KdTree,
        SpatialBackend::BruteForce,
    ] {
        let rows = hydrogen_bond_candidates_over_segment(
            &structure,
            &start,
            &end,
            HydrogenBondCandidateOptions {
                backend,
                ..options()
            },
            &ExecutionContext::default(),
        )
        .unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(
            rows.row(0).unwrap(),
            HydrogenBondCandidate {
                donor: AtomIndex::new(0),
                hydrogen: AtomIndex::new(1),
                acceptor: AtomIndex::new(2)
            }
        );
    }
}

#[test]
fn common_translation_and_hydrogen_motion_do_not_expand_the_heavy_atom_query() {
    let structure = annotated_structure(SOURCE);
    let start = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [6.0, 0.0, 0.0]];
    let end = [[1000.0, 0.0, 0.0], [-10000.0, 0.0, 0.0], [1006.0, 0.0, 0.0]];
    let roles = role_selection(&structure, molframe_core::HBOND_DONOR_ANNOTATION).unwrap();
    let acceptors = role_selection(&structure, molframe_core::HBOND_ACCEPTOR_ANNOTATION).unwrap();
    assert!(envelope_cutoff(&start, &end, &roles, &acceptors, 3.5).unwrap() < 3.50001);
    assert!(
        hydrogen_bond_candidates_over_segment(
            &structure,
            &start,
            &end,
            options(),
            &ExecutionContext::default()
        )
        .unwrap()
        .is_empty()
    );
}

#[test]
fn motion_candidates_reject_invalid_coordinate_columns_and_missing_chemistry() {
    let input = molframe_core::io::InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let unannotated = molframe_cif::read(&input, &molframe_core::io::ReadOptions::new())
        .expect("unannotated source")
        .0;
    assert!(matches!(
        hydrogen_bond_candidates_over_segment(
            &unannotated,
            unannotated.positions(),
            unannotated.positions(),
            options(),
            &ExecutionContext::default(),
        ),
        Err(HydrogenBondError::MissingChemistry)
    ));
    let structure = annotated_structure(SOURCE);
    let start = structure.positions();
    let mut end = start.to_vec();
    end[2][0] = f32::NAN;
    assert!(matches!(
        hydrogen_bond_candidates_over_segment(
            &structure,
            start,
            &end,
            options(),
            &ExecutionContext::default()
        ),
        Err(HydrogenBondError::InvalidOptions)
    ));
    assert!(matches!(
        hydrogen_bond_candidates_over_segment(
            &structure,
            &start[..2],
            start,
            options(),
            &ExecutionContext::default()
        ),
        Err(HydrogenBondError::InvalidOptions)
    ));
    assert!(matches!(
        hydrogen_bond_candidates_over_segment(
            &structure,
            start,
            start,
            HydrogenBondCandidateOptions {
                maximum_candidates: 0,
                ..options()
            },
            &ExecutionContext::default()
        ),
        Err(HydrogenBondError::InvalidOptions)
    ));
}

#[test]
fn a_candidate_budget_rejects_the_whole_result_and_parallel_order_is_stable() {
    let structure = water_lattice(3);
    let positions = structure.positions();
    let policy = HydrogenBondCandidateOptions {
        maximum_candidates: 1000,
        ..options()
    };
    let serial = hydrogen_bond_candidates_over_segment(
        &structure,
        positions,
        positions,
        policy,
        &ExecutionContext::builder()
            .worker_budget(1)
            .build()
            .unwrap(),
    )
    .unwrap();
    assert!(serial.len() > 1);
    let parallel = hydrogen_bond_candidates_over_segment(
        &structure,
        positions,
        positions,
        policy,
        &ExecutionContext::builder()
            .worker_budget(4)
            .build()
            .unwrap(),
    )
    .unwrap();
    assert_eq!(serial, parallel);
    assert!(matches!(
        hydrogen_bond_candidates_over_segment(
            &structure,
            positions,
            positions,
            HydrogenBondCandidateOptions {
                maximum_candidates: 1,
                ..policy
            },
            &ExecutionContext::default()
        ),
        Err(HydrogenBondError::CandidateBudget { maximum: 1 })
    ));
}
