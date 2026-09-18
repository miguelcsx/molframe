use super::*;
use crate::SpatialBackend;

fn lattice(side: u16) -> Vec<[f32; 3]> {
    let mut positions = Vec::new();
    for x in 0..side {
        for y in 0..side {
            for z in 0..side {
                positions.push([f32::from(x), f32::from(y), f32::from(z)]);
            }
        }
    }
    positions
}

fn counted(positions: &[[f32; 3]], workers: usize, backend: SpatialBackend) -> usize {
    let context = ExecutionContext::builder()
        .worker_budget(workers)
        .build()
        .expect("worker context is valid");
    let all = AtomSelection::All(u32::try_from(positions.len()).expect("lattice fits u32"));
    let query = PairQuery {
        positions,
        left: &all,
        right: &all,
        cutoff: 1.5,
        options: SpatialSearchOptions::with_backend(backend),
        periodic: None,
        context: &context,
    };
    let parts = reduce_pairs_within_unsorted(&query, || 0usize, |count, _| *count += 1)
        .expect("reduction succeeds");
    parts.into_iter().sum()
}

#[test]
fn the_pair_count_is_identical_at_every_worker_count() {
    let positions = lattice(8);
    let serial = counted(&positions, 1, SpatialBackend::CellList);

    for workers in [2, 3, 4, 8, 16] {
        assert_eq!(
            counted(&positions, workers, SpatialBackend::CellList),
            serial,
            "worker count {workers} changed the pair count"
        );
    }
}

#[test]
fn the_reduction_agrees_with_the_materialising_query() {
    let positions = lattice(6);
    let all = AtomSelection::All(u32::try_from(positions.len()).expect("lattice fits u32"));
    let expected = crate::pairs_within(
        &positions,
        &all,
        &all,
        1.5,
        SpatialBackend::CellList,
        None,
        &ExecutionContext::default(),
    )
    .expect("query succeeds")
    .len();

    assert_eq!(counted(&positions, 4, SpatialBackend::CellList), expected);
}

#[test]
fn a_backend_without_a_block_decomposition_still_reduces() {
    let positions = lattice(4);

    assert_eq!(
        counted(&positions, 4, SpatialBackend::BruteForce),
        counted(&positions, 1, SpatialBackend::CellList)
    );
}

fn cross_count(positions: &[[f32; 3]], workers: usize, backend: SpatialBackend) -> usize {
    let context = ExecutionContext::builder()
        .worker_budget(workers)
        .build()
        .expect("worker context is valid");
    let left = AtomSelection::from_sorted((0..96).collect());
    let right = AtomSelection::from_sorted((96..256).collect());
    let query = PairQuery {
        positions,
        left: &left,
        right: &right,
        cutoff: 1.5,
        options: SpatialSearchOptions::with_backend(backend),
        periodic: None,
        context: &context,
    };
    let parts = reduce_pairs_within_unsorted(&query, || 0usize, |count, _| *count += 1)
        .expect("reduction succeeds");
    parts.into_iter().sum()
}

#[test]
fn a_cross_selection_query_counts_the_same_at_every_worker_count() {
    let positions = lattice(8);
    let serial = cross_count(&positions, 1, SpatialBackend::CellList);
    assert!(serial > 0, "the fixture must produce cross pairs");

    for workers in [2, 3, 4, 8, 16] {
        assert_eq!(
            cross_count(&positions, workers, SpatialBackend::CellList),
            serial,
            "worker count {workers} changed the cross-selection pair count"
        );
    }
}

#[test]
fn a_blocked_cross_selection_agrees_with_the_materialising_query() {
    let positions = lattice(8);
    let left = AtomSelection::from_sorted((0..96).collect());
    let right = AtomSelection::from_sorted((96..256).collect());
    let expected = crate::pairs_within(
        &positions,
        &left,
        &right,
        1.5,
        SpatialBackend::CellList,
        None,
        &ExecutionContext::default(),
    )
    .expect("query succeeds")
    .len();

    assert_eq!(
        cross_count(&positions, 8, SpatialBackend::CellList),
        expected
    );
}
