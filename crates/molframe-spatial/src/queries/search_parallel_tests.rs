use super::*;
use molframe_core::ExecutionContext;

fn context() -> ExecutionContext {
    ExecutionContext::default()
}

fn pair_indices(pairs: &[NeighborPair]) -> Vec<(u32, u32)> {
    pairs.iter().map(|pair| (pair.first, pair.second)).collect()
}

// A dense cubic lattice forces the cell-list backend and many non-empty cells,
// so the parallel partition spans multiple cells per worker.
fn lattice(side: u8) -> Vec<[f32; 3]> {
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

#[test]
fn parallel_pairs_match_serial_for_every_worker_count() {
    let positions = lattice(16); // 4096 atoms -> cell-list backend
    let all = AtomSelection::All(4096);
    let options = SpatialSearchOptions::with_backend(SpatialBackend::CellList);
    let serial = pairs_within_with_options(&positions, &all, &all, 1.5, options, None, &context())
        .expect("serial cell-list search is valid");
    for workers in [1, 2, 3, 4, 8, 16] {
        let worker_context = ExecutionContext::builder()
            .worker_budget(workers)
            .build()
            .expect("worker context is valid");
        let parallel =
            pairs_within_with_options(&positions, &all, &all, 1.5, options, None, &worker_context)
                .expect("parallel cell-list search is valid");
        assert_eq!(
            pair_indices(&serial),
            pair_indices(&parallel),
            "worker count {workers} changed the pair set"
        );
        // Distances must match too, not only the index pairs.
        assert_eq!(
            serial
                .iter()
                .map(|p| p.distance_squared.to_bits())
                .collect::<Vec<_>>(),
            parallel
                .iter()
                .map(|p| p.distance_squared.to_bits())
                .collect::<Vec<_>>(),
            "worker count {workers} changed a squared distance"
        );
    }
}

fn budgeted(workers: usize, bytes: usize) -> ExecutionContext {
    ExecutionContext::builder()
        .worker_budget(workers)
        .memory_budget(molframe_core::MemoryBudget::new(bytes).expect("a positive budget"))
        .scratch_policy(molframe_core::ScratchPolicy::new(0))
        .build()
        .expect("a valid context")
}

#[test]
fn a_bounded_budget_streams_the_pair_blocks_and_matches_the_unbounded_run() {
    let positions = lattice(16);
    let all = AtomSelection::All(4096);
    let options = SpatialSearchOptions::with_backend(SpatialBackend::CellList);
    let reference = pairs_within_with_options(
        &positions,
        &all,
        &all,
        1.5,
        options,
        None,
        &budgeted(1, 512 * 1024 * 1024),
    )
    .expect("serial run");
    for workers in [2, 4, 8] {
        let context = budgeted(workers, 64 * 1024 * 1024);
        let streamed =
            pairs_within_with_options(&positions, &all, &all, 1.5, options, None, &context)
                .expect("a bounded budget suffices");
        assert_eq!(
            pair_indices(&streamed),
            pair_indices(&reference),
            "{workers} workers"
        );
        assert_eq!(context.reserved_bytes(), 0);
    }
}

#[test]
fn a_budget_below_one_pair_block_is_refused() {
    let positions = lattice(16);
    let all = AtomSelection::All(4096);
    let options = SpatialSearchOptions::with_backend(SpatialBackend::CellList);
    let tiny = budgeted(4, 1024);
    assert!(pairs_within_with_options(&positions, &all, &all, 1.5, options, None, &tiny).is_err());
}

#[test]
fn sparse_index_validation_borrows_the_existing_allocation() {
    let selection = AtomSelection::Sparse(vec![1, 5, 9]);
    let indices = super::checked_indices(&selection, 10).expect("valid indices");
    let AtomSelection::Sparse(original) = &selection else {
        panic!("sparse input")
    };
    assert!(matches!(indices, std::borrow::Cow::Borrowed(_)));
    assert_eq!(indices.as_ptr(), original.as_ptr());
    assert_eq!(
        super::checked_indices(&selection, 9),
        Err(SpatialError::AtomOutOfBounds(9))
    );
}
