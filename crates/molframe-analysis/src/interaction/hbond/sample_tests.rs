use super::*;
use crate::hbond::tests::{SOURCE, annotated_structure};
use crate::{
    HydrogenBondCandidate, HydrogenBondCandidateOptions, hydrogen_bond_candidates_over_segment,
};
use molframe_core::{AtomIndex, MemoryBudget};
use molframe_spatial::SpatialBackend;

fn options() -> HydrogenBondOptions {
    HydrogenBondOptions {
        maximum_donor_acceptor_distance: 3.5,
        minimum_angle_degrees: 150.0,
        backend: SpatialBackend::CellList,
        periodic: false,
    }
}

#[test]
fn sampled_snapshots_find_contacts_inside_an_interval_and_share_static_geometry_rules() {
    let structure = annotated_structure(SOURCE);
    let start = [[0.0; 3], [1.0, 0.0, 0.0], [6.0, 0.0, 0.0]];
    let end = [[0.0; 3], [1.0, 0.0, 0.0], [-6.0, 0.0, 0.0]];
    let context = ExecutionContext::default();
    let candidates = hydrogen_bond_candidates_over_segment(
        &structure,
        &start,
        &end,
        HydrogenBondCandidateOptions {
            maximum_donor_acceptor_distance: 3.5,
            maximum_candidates: 100,
            backend: SpatialBackend::CellList,
        },
        &context,
    )
    .expect("conservative interval");
    for (weight, count) in [(0.0, 0), (0.25, 1), (0.5, 0), (0.75, 0), (1.0, 0)] {
        let table = hydrogen_bonds_at_sample(
            &structure,
            CoordinateSample::new(&start, &end, weight).expect("sample"),
            &candidates,
            options(),
            HydrogenBondPolicy::default(),
            &context,
        )
        .expect("snapshot");
        assert_eq!(table.len(), count);
        if let Some(bond) = table.row(0) {
            assert!((bond.donor_acceptor_distance - 3.0).abs() < 1e-6);
            assert!((bond.hydrogen_acceptor_distance - 2.0).abs() < 1e-6);
            assert!((bond.angle_degrees - 180.0).abs() < 1e-6);
        }
    }
    let sample = CoordinateSample::new(structure.positions(), structure.positions(), 0.5)
        .expect("static sample");
    let sampled = hydrogen_bonds_at_sample(
        &structure,
        sample,
        &candidates,
        options(),
        HydrogenBondPolicy::default(),
        &context,
    )
    .expect("same source");
    let direct = crate::hydrogen_bonds(&structure, options(), &context).expect("direct detection");
    assert_eq!(
        sampled.iter().collect::<Vec<_>>(),
        direct.iter().collect::<Vec<_>>()
    );
    let strict = hydrogen_bonds_at_sample(
        &structure,
        sample,
        &candidates,
        options(),
        HydrogenBondPolicy {
            maximum_hydrogen_acceptor_distance: Some(1.0),
            ..HydrogenBondPolicy::default()
        },
        &context,
    )
    .expect("hydrogen distance");
    assert!(strict.is_empty());
}

#[test]
fn malformed_samples_and_forged_chemical_triples_never_publish_partial_tables() {
    let structure = annotated_structure(SOURCE);
    let coordinates = structure.positions();
    let context = ExecutionContext::default();
    for [donor, hydrogen, acceptor] in [[0, 1, 99], [1, 0, 2], [0, 2, 1]] {
        let candidates = [HydrogenBondCandidate {
            donor: AtomIndex::new(donor),
            hydrogen: AtomIndex::new(hydrogen),
            acceptor: AtomIndex::new(acceptor),
        }]
        .into_iter()
        .collect();
        assert!(matches!(
            hydrogen_bonds_at_sample(
                &structure,
                CoordinateSample::new(coordinates, coordinates, 0.5).expect("sample"),
                &candidates,
                options(),
                HydrogenBondPolicy::default(),
                &context,
            ),
            Err(HydrogenBondError::InvalidSample)
        ));
    }
    let candidates = HydrogenBondCandidateTable::default();
    for coordinates in [vec![[0.0; 3]], vec![[f32::NAN; 3]; 3]] {
        assert!(matches!(
            hydrogen_bonds_at_sample(
                &structure,
                CoordinateSample::new(&coordinates, &coordinates, 0.5).expect("metadata"),
                &candidates,
                options(),
                HydrogenBondPolicy::default(),
                &context,
            ),
            Err(HydrogenBondError::InvalidSample)
        ));
    }
}

#[test]
fn snapshots_honor_cancellation_and_reject_the_whole_output_when_memory_is_exhausted() {
    let structure = annotated_structure(SOURCE);
    let sample =
        CoordinateSample::new(structure.positions(), structure.positions(), 0.5).expect("sample");
    let candidates = [HydrogenBondCandidate {
        donor: AtomIndex::new(0),
        hydrogen: AtomIndex::new(1),
        acceptor: AtomIndex::new(2),
    }]
    .into_iter()
    .collect();
    let cancelled = ExecutionContext::default();
    cancelled.cancellation().cancel();
    assert!(matches!(
        hydrogen_bonds_at_sample(
            &structure,
            sample,
            &candidates,
            options(),
            HydrogenBondPolicy::default(),
            &cancelled
        ),
        Err(HydrogenBondError::Spatial(SpatialError::Cancelled))
    ));
    let bounded = ExecutionContext::builder()
        .memory_budget(MemoryBudget::new(16_000_000).expect("budget"))
        .build()
        .expect("context");
    let _held = bounded
        .try_reserve(16_000_000)
        .expect("other retained work");
    assert!(matches!(
        hydrogen_bonds_at_sample(
            &structure,
            sample,
            &candidates,
            options(),
            HydrogenBondPolicy::default(),
            &bounded
        ),
        Err(HydrogenBondError::Spatial(SpatialError::Memory(_)))
    ));
}
